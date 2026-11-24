"""
Immutable 100 kW_th Industrial High-Temperature Heat Pump (HTHP) Skid Benchmark
across 3 standardized <150 C decarbonization duties:
  - DAIRY_90C     : Sink 35 -> 90 C (55 K), Waste-heat source 45 -> 20 C (25 K)
  - DISTRICT_110C : Sink 40 -> 110 C (70 K), Waste-heat source 55 -> 25 C (30 K)
  - CHEMICAL_135C : Sink 60 -> 135 C (75 K), Waste-heat source 75 -> 40 C (35 K)

Computes exact finite-volume BPHE heat transfer (with mandatory Silver-Bell-Ghaly
and Han/Yan two-phase correlations), water pumping power, volumetric heating
capacity (VHC), Gouy-Stodola exergy destruction, and Levelized Cost of Heat (LCOH).

Hard industrial safety gates:
  - Fixed total BPHE area <= 18.0 m2 (condenser + evaporator)
  - No sub-atmospheric vacuum suction: P_evap_out >= 1.05 bar (ATEX safety)
  - Discharge pressure <= 38.0 bar, T_dis <= 155 C, dry compression (>= 4 K SH)
  - Zero PFAS & GWP <= 5 for natural tracks
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from brazed_plate_hx import BPHEGeometry, BPHEResult, _water_props, solve_counterflow_bphe_area
from thermo_backend import FLUID_METADATA, MixtureBackend

# Locked industrial skid limits (100 kW_th nominal heating output)
Q_DUTY_W = 100_000.0
MAX_TOTAL_BPHE_AREA_M2 = 18.0
MIN_SUCTION_PRESSURE_BAR = 1.05   # Strict above-atmospheric ATEX suction gate
MAX_COND_PRESSURE_BAR = 38.0      # Standard industrial semi-hermetic class
MAX_DISCHARGE_TEMP_C = 155.0      # Synthetic POE/PAG industrial oil limit
MIN_DISCHARGE_SH_K = 4.0          # Strict dry-compression guard

# Techno-economic assumptions (transparent, realistic industrial Europe/US)
ELEC_PRICE_USD_MWH = 85.0
OPER_HOURS_YR = 6000.0
CRF_15YR_7PCT = 0.07 * (1.07**15) / (1.07**15 - 1.0)  # 0.10979


@dataclass(frozen=True)
class IndustrialDuty:
    name: str
    T_sink_in_C: float
    T_sink_out_C: float
    T_src_in_C: float
    T_src_out_C: float
    Q_cond_W: float = Q_DUTY_W


INDUSTRIAL_DUTIES: Dict[str, IndustrialDuty] = {
    "DAIRY_90C": IndustrialDuty("DAIRY_90C", 35.0, 90.0, 45.0, 20.0),
    "DISTRICT_110C": IndustrialDuty("DISTRICT_110C", 40.0, 110.0, 55.0, 25.0),
    "CHEMICAL_135C": IndustrialDuty("CHEMICAL_135C", 60.0, 135.0, 75.0, 40.0),
}


@dataclass(frozen=True)
class HTHPDesign:
    label: str
    duty_name: str
    # Layer 1: Discrete flowsheet architecture
    topology: str                  # BASIC_1STAGE | IHX_SUBCOOL | FTVI_ECON | FTVI_IHX
    ihx_eff: float                 # 0.0..0.82
    ft_pr_bias: float              # 0.78..1.28 multiplier on sqrt(Pe*Pc)
    # Layer 2: Natural working fluid & mole fractions
    components: Tuple[str, ...]
    mole_fractions: Tuple[float, ...]
    # Layer 3: Asymmetric BPHE allocation & chevron geometry
    cond_area_frac: float          # 0.30..0.75 of 18.0 m2 total area
    cond_hx: BPHEGeometry
    evap_hx: BPHEGeometry
    # Layer 4: Cycle control setpoints
    superheat_K: float             # >= 5.0 K (higher for retrograde pentanes/hexanes)
    subcool_K: float               # >= 2.0 K


@dataclass
class HTHPEvaluation:
    label: str
    duty_name: str
    feasible: bool
    violations: List[str]
    COP_sys: float                 # Including condenser + evaporator water pump power!
    COP_lorenz_ideal: float        # Exact ideal Lorenz COP between sink & source glides
    eta_lorenz: float              # Second-law Lorenz efficiency = COP_sys / COP_lorenz_ideal
    LCOH_usd_MWh: float            # Levelized Cost of Heat ($/MWh_th)
    CAPEX_usd: float
    VHC_kJ_m3: float               # Volumetric Heating Capacity (kJ/m3 suction)
    V_suc_m3_h: float              # Required compressor suction volume flow (m3/h)
    W_comp_kW: float
    W_pumps_kW: float
    P_evap_bar: float
    P_evap_out_bar: float
    P_cond_bar: float
    PR: float
    T_dis_C: float
    dis_SH_K: float
    glide_cond_K: float
    glide_evap_K: float
    sbg_cond: float
    sbg_evap: float
    dP_evap_kPa: float
    dP_cond_kPa: float
    pinch_cond_K: float
    pinch_evap_K: float
    gwp: float
    is_pfas: bool
    # Gouy-Stodola Exergy Destruction (kW)
    I_comp_kW: float
    I_cond_kW: float
    I_evap_kW: float
    I_exp_kW: float
    I_ihx_ft_kW: float
    I_pumps_kW: float
    # Auxiliary VLE info for Dual-EOS audit
    z_evap: Tuple[float, ...] = field(default_factory=tuple)


def _industrial_comp_eff(PR: float) -> Tuple[float, float, float]:
    """Industrial semi-hermetic screw/piston compressor efficiency map (Bitzer/Mycom)."""
    PR_c = max(1.3, min(14.0, PR))
    eta_is = max(0.50, min(0.81, 0.80 - 0.018 * ((PR_c - 3.6) ** 2) / PR_c))
    eta_v = max(0.62, min(0.95, 0.96 - 0.024 * ((PR_c ** 0.82) - 1.0)))
    eta_em = 0.93  # Industrial IE4 motor + VFD
    return eta_is, eta_v, eta_em


def evaluate_hthp(
    design: HTHPDesign,
    enforce_natural_zero_pfas: bool = True,
    enforce_above_atm: bool = True,
    n_seg: int = 24,
    htc_mult: float = 1.0,
    sbg_mult: float = 1.0,
    dp_mult: float = 1.0,
) -> HTHPEvaluation:
    duty = INDUSTRIAL_DUTIES[design.duty_name]
    violations: List[str] = []

    tot_area = design.cond_hx.area_m2 + design.evap_hx.area_m2
    if tot_area > MAX_TOTAL_BPHE_AREA_M2 + 1e-4:
        violations.append(f"Total BPHE area {tot_area:.2f}m2 > {MAX_TOTAL_BPHE_AREA_M2:.1f}m2")

    is_pfas = any("PFAS" in str(FLUID_METADATA[c]["safety"]) for c in design.components)
    try:
        bk_cond = MixtureBackend(design.components, design.mole_fractions)
    except Exception as e:
        return _infeasible_hthp(design, [f"EOS init failed: {e}"], is_pfas, 9999.0)

    if enforce_natural_zero_pfas:
        if bk_cond.gwp > 5.0:
            violations.append(f"GWP {bk_cond.gwp:.1f} > 5.0")
        if is_pfas:
            violations.append("Contains restricted PFAS/HFO refrigerant")

    # Exact ideal Lorenz COP between log-mean sink and source temperatures
    T_sk_in, T_sk_out = duty.T_sink_in_C + 273.15, duty.T_sink_out_C + 273.15
    T_sr_in, T_sr_out = duty.T_src_in_C + 273.15, duty.T_src_out_C + 273.15
    T_lm_sink = (T_sk_out - T_sk_in) / math.log(T_sk_out / T_sk_in)
    T_lm_src = (T_sr_in - T_sr_out) / math.log(T_sr_in / T_sr_out)
    cop_lorenz = T_lm_sink / max(1.0, T_lm_sink - T_lm_src)
    T_dead = T_sr_out

    # Fixed process sink water mass flow rate to hit exact (T_sink_in -> T_sink_out) at 100 kW
    _, cp_sk, _, _ = _water_props(0.5 * (T_sk_in + T_sk_out))
    m_dot_w_sink = duty.Q_cond_W / (cp_sk * (T_sk_out - T_sk_in))
    _, cp_sr, _, _ = _water_props(0.5 * (T_sr_in + T_sr_out))

    sh_K = max(5.0, design.superheat_K)
    sc_K = max(2.0, design.subcool_K)
    use_ftvi = design.topology in ("FTVI_ECON", "FTVI_IHX")
    use_ihx = design.topology in ("IHX_SUBCOOL", "FTVI_IHX")
    ihx_eff = max(0.0, min(0.82, design.ihx_eff)) if use_ihx else 0.0

    # Robust initial pressure brackets guaranteed ABOVE sink pinch and BELOW source pinch!
    try:
        # Pure/low-glide fluids need P_cond near T_sk_out and P_evap near T_sr_out;
        # high-glide mixtures sit closer to midpoint temperatures.
        P_evap = bk_cond.p_for_t_mid(0.75 * T_sr_out + 0.25 * T_sr_in - 3.0)
        P_cond = bk_cond.p_for_t_mid(0.75 * T_sk_out + 0.25 * T_sk_in + 4.0)
    except Exception as e:
        return _infeasible_hthp(design, [f"Initial P bracket failed: {e}"], is_pfas, bk_cond.gwp)

    ln_pe_lo = math.log(0.12e5)
    ln_pe_hi = math.log(24e5)
    ln_pc_hi = math.log(min(MAX_COND_PRESSURE_BAR * 1e5 * 1.15, 0.94 * bk_cond.Pc))

    P_evap = max(0.15e5, min(20e5, P_evap))
    P_cond = max(P_evap * 1.45, min(math.exp(ln_pc_hi) * 0.96, P_cond))

    cache = {"state": None}

    def eval_res(x_vec: np.ndarray) -> np.ndarray:
        ln_pe = float(np.clip(x_vec[0], ln_pe_lo, ln_pe_hi))
        ln_pc = float(np.clip(x_vec[1], ln_pe + math.log(1.30), ln_pc_hi))
        Pe = math.exp(ln_pe)
        Pc = math.exp(ln_pc)

        sl_cond = bk_cond.get_isobar(Pc, nq=17)
        T_cond_out = sl_cond.Tb - sc_K
        h_cond_out, s_cond_out, _ = bk_cond.liquid_state_PT(Pc, T_cond_out)

        if use_ftvi:
            P_mid = design.ft_pr_bias * math.sqrt(Pe * Pc)
            P_mid = max(Pe * 1.15, min(Pc * 0.85, P_mid))
            T_ft, q_ft, x_l_ft, y_v_ft, h_ft_l, h_ft_v = bk_cond.flash_phase_split(P_mid, h_cond_out)
            q_ft = max(0.02, min(0.55, q_ft))
            bk_evap = MixtureBackend(design.components, tuple(float(v) for v in x_l_ft)) if len(design.components) > 1 else bk_cond
            h_pre_eev2 = h_ft_l
            T_pre_eev2 = T_ft
            z_ev = tuple(float(v) for v in x_l_ft)
        else:
            P_mid = Pe
            q_ft = 0.0
            bk_evap = bk_cond
            h_pre_eev2 = h_cond_out
            T_pre_eev2 = T_cond_out
            h_ft_v = 0.0
            y_v_ft = np.array(design.mole_fractions)
            z_ev = design.mole_fractions

        sl_evap = bk_evap.get_isobar(Pe, nq=17)
        T_evap_out = sl_evap.Td + sh_K
        h_evap_out, s_evap_out, rho_evap_out = bk_evap.vapor_state_PT(Pe, T_evap_out)

        if use_ihx and ihx_eff > 1e-3 and T_pre_eev2 > T_evap_out + 2.0:
            C_min = min(sl_evap.cp_v, sl_cond.cp_l)
            dh_ihx = ihx_eff * C_min * (T_pre_eev2 - T_evap_out)
            dh_ihx = min(dh_ihx, 0.88 * sl_cond.cp_l * max(0.0, T_pre_eev2 - T_evap_out - 2.0))
            h_suc = h_evap_out + dh_ihx
            T_suc = T_evap_out + dh_ihx / max(sl_evap.cp_v, 500.0)
            h_evap_in = h_pre_eev2 - dh_ihx
            _, s_suc, rho_suc = bk_evap.vapor_state_PT(Pe, T_suc)
        else:
            dh_ihx = 0.0
            h_suc = h_evap_out
            T_suc = T_evap_out
            s_suc = s_evap_out
            rho_suc = rho_evap_out
            h_evap_in = h_pre_eev2

        if not use_ftvi:
            eta_is, eta_v, eta_em = _industrial_comp_eff(Pc / Pe)
            h_dis_s, _ = bk_cond.state_PS(Pc, s_suc)
            w_shaft = max(1e3, (h_dis_s - h_suc) / eta_is)
            h_cond_in = h_suc + w_shaft
            T_dis, s_dis, _ = bk_cond.state_PH(Pc, h_cond_in)
            dh_cond = max(25e3, h_cond_in - h_cond_out)
            m_dot_cond = duty.Q_cond_W / dh_cond
            m_dot_evap = m_dot_cond
            W_comp_elec = m_dot_evap * w_shaft / eta_em
            s_gen_comp = m_dot_evap * max(0.0, s_dis - s_suc)
            s_gen_ihx_ft = 0.0
            if dh_ihx > 10.0:
                T_l_avg = T_pre_eev2 - 0.5 * dh_ihx / max(sl_cond.cp_l, 500.0)
                T_v_avg = 0.5 * (T_evap_out + T_suc)
                s_gen_ihx_ft += m_dot_evap * dh_ihx * max(0.0, 1.0 / T_v_avg - 1.0 / T_l_avg)
            Q_evap_req = m_dot_evap * max(10e3, h_evap_out - h_evap_in)
        else:
            PR1, PR2 = P_mid / Pe, Pc / P_mid
            eta_is1, eta_v, eta_em = _industrial_comp_eff(PR1)
            eta_is2, _, _ = _industrial_comp_eff(PR2)
            h_mid_s, _ = bk_evap.state_PS(P_mid, s_suc)
            w1 = max(500.0, (h_mid_s - h_suc) / eta_is1)
            h_mid_1 = h_suc + w1
            _, s_mid_1, _ = bk_evap.state_PH(P_mid, h_mid_1)

            inj_r = q_ft / max(1.0 - q_ft, 0.2)
            h_mid_mix = (h_mid_1 + inj_r * h_ft_v) / (1.0 + inj_r)
            _, s_mid_mix, _ = bk_cond.state_PH(P_mid, h_mid_mix)
            try:
                bk_v = MixtureBackend(design.components, tuple(float(v) for v in y_v_ft)) if len(design.components) > 1 else bk_cond
                _, s_ft_v, _ = bk_v.state_PH(P_mid, h_ft_v)
            except Exception:
                s_ft_v = s_mid_mix

            h_dis_s, _ = bk_cond.state_PS(Pc, s_mid_mix)
            w2 = max(500.0, (h_dis_s - h_mid_mix) / eta_is2)
            h_cond_in = h_mid_mix + w2
            T_dis, s_dis, _ = bk_cond.state_PH(Pc, h_cond_in)

            dh_cond = max(25e3, h_cond_in - h_cond_out)
            m_dot_cond = duty.Q_cond_W / dh_cond
            m_dot_evap = m_dot_cond / (1.0 + inj_r)
            W_comp_elec = (m_dot_evap * w1 + m_dot_cond * w2) / eta_em
            s_gen_comp = m_dot_evap * max(0.0, s_mid_1 - s_suc) + m_dot_cond * max(0.0, s_dis - s_mid_mix)
            s_gen_ihx_ft = max(0.0, m_dot_cond * s_mid_mix - m_dot_evap * s_mid_1 - (m_dot_cond - m_dot_evap) * s_ft_v)
            if dh_ihx > 10.0:
                T_l_avg = T_pre_eev2 - 0.5 * dh_ihx / max(sl_cond.cp_l, 500.0)
                T_v_avg = 0.5 * (T_evap_out + T_suc)
                s_gen_ihx_ft += m_dot_evap * dh_ihx * max(0.0, 1.0 / T_v_avg - 1.0 / T_l_avg)
            Q_evap_req = m_dot_evap * max(10e3, h_evap_out - h_evap_in)

        m_dot_w_src = Q_evap_req / (cp_sr * (T_sr_in - T_sr_out))

        res_c = solve_counterflow_bphe_area(
            design.cond_hx, bk_cond, True, Pc, h_cond_in, h_cond_out, m_dot_cond,
            T_sk_in, T_sk_out, m_dot_w_sink, n_seg=n_seg, htc_mult=htc_mult, sbg_mult=sbg_mult, dp_mult=dp_mult,
        )
        res_e = solve_counterflow_bphe_area(
            design.evap_hx, bk_evap, False, Pe, h_evap_in, h_evap_out, m_dot_evap,
            T_sr_in, T_sr_out, m_dot_w_src, n_seg=n_seg, htc_mult=htc_mult, sbg_mult=sbg_mult, dp_mult=dp_mult,
        )

        # Strictly monotone log-area residual:
        # A_req > A_avail -> err > 0 -> condenser needs HIGHER P_cond; evaporator needs LOWER P_evap!
        err_e = math.log(max(0.05, res_e.A_req_m2) / design.evap_hx.area_m2)
        err_c = math.log(max(0.05, res_c.A_req_m2) / design.cond_hx.area_m2)

        cache["state"] = (
            sl_cond, sl_evap, res_c, res_e, W_comp_elec, T_dis,
            m_dot_evap, m_dot_cond, rho_suc, eta_v, s_gen_comp, s_gen_ihx_ft,
            s_cond_out, h_evap_in, h_cond_out, bk_evap, Pe, Pc, z_ev, dh_ihx,
        )
        return np.array([err_e, err_c], dtype=float)

    try:
        # Rock-solid alternating monotone 1D bracketed bisection on ln(Pe) and ln(Pc):
        # - err_e(Pe) = ln(A_req_e / A_e) is STRICTLY INCREASING in Pe (higher Pe -> smaller dT -> larger A_req_e)
        # - err_c(Pc) = ln(A_req_c / A_c) is STRICTLY DECREASING in Pc (higher Pc -> larger dT -> smaller A_req_c)
        from scipy.optimize import brentq

        ln_pe_cur = math.log(P_evap)
        ln_pc_cur = math.log(P_cond)

        def solve_pc_given_pe(ln_pe_fixed: float) -> float:
            lo = ln_pe_fixed + math.log(1.30)
            hi = ln_pc_hi
            f_hi = float(eval_res(np.array([ln_pe_fixed, hi]))[1])
            if f_hi > 0.0:
                # Even at max subcritical Pc, condenser area is insufficient (fluid Tc too low!)
                return hi
            # Walk down from hi to bracket root tightly
            f_lo = float(eval_res(np.array([ln_pe_fixed, lo]))[1])
            if f_lo <= 0.0:
                return lo
            return float(brentq(lambda lpc: float(eval_res(np.array([ln_pe_fixed, lpc]))[1]), lo, hi, xtol=1e-3, maxiter=22))

        def solve_pe_given_pc(ln_pc_fixed: float) -> float:
            lo = ln_pe_lo
            hi = min(ln_pe_hi, ln_pc_fixed - math.log(1.30))
            if hi <= lo + 0.02:
                return lo
            f_lo = float(eval_res(np.array([lo, ln_pc_fixed]))[0])
            if f_lo > 0.0:
                return lo
            f_hi = float(eval_res(np.array([hi, ln_pc_fixed]))[0])
            if f_hi <= 0.0:
                return hi
            return float(brentq(lambda lpe: float(eval_res(np.array([lpe, ln_pc_fixed]))[0]), lo, hi, xtol=1e-3, maxiter=22))

        converged = False
        for outer_k in range(5):
            ln_pc_new = solve_pc_given_pe(ln_pe_cur)
            ln_pe_new = solve_pe_given_pc(ln_pc_new)
            r_chk = eval_res(np.array([ln_pe_new, ln_pc_new]))
            ln_pe_cur, ln_pc_cur = ln_pe_new, ln_pc_new
            if max(abs(r_chk[0]), abs(r_chk[1])) < 0.012:
                converged = True
                break
        if not converged:
            r_chk = eval_res(np.array([ln_pe_cur, ln_pc_cur]))
            converged = max(abs(r_chk[0]), abs(r_chk[1])) < 0.025
    except Exception as e:
        return _infeasible_hthp(design, [f"Solver exception: {e}"], is_pfas, bk_cond.gwp)

    if cache["state"] is None or not converged:
        return _infeasible_hthp(design, ["BPHE-cycle balance did not converge"], is_pfas, bk_cond.gwp)

    (
        sl_cond, sl_evap, res_c, res_e, W_comp_elec, T_dis,
        m_dot_evap, m_dot_cond, rho_suc, eta_v, s_gen_comp, s_gen_ihx_ft,
        s_cond_out, h_evap_in, h_cond_out, bk_evap, Pe, Pc, z_ev, dh_ihx,
    ) = cache["state"]

    # Account for evaporator suction pressure drop penalty on compressor power & suction density!
    P_suc_out = max(0.10e5, Pe - res_e.dP_ref_Pa)
    dp_comp_fac = math.log(Pc / P_suc_out) / max(0.15, math.log(Pc / Pe))
    W_comp_elec *= max(1.0, min(1.40, dp_comp_fac))
    rho_suc_actual = rho_suc * (P_suc_out / Pe)

    # Volumetric Heating Capacity & required industrial compressor swept flow (m3/h)
    V_suc_m3_s = m_dot_evap / max(rho_suc_actual * eta_v, 1e-6)
    V_suc_m3_h = V_suc_m3_s * 3600.0
    VHC_kJ_m3 = (duty.Q_cond_W / 1000.0) / max(V_suc_m3_s, 1e-6)

    W_pumps = res_c.W_pump_W + res_e.W_pump_W
    W_tot = W_comp_elec + W_pumps
    COP_sys = duty.Q_cond_W / max(W_tot, 100.0)
    eta_lorenz = COP_sys / cop_lorenz

    # Hard industrial safety & dry-compression checks
    dis_SH_K = T_dis - sl_cond.Td
    if enforce_above_atm and (P_suc_out / 1e5) < MIN_SUCTION_PRESSURE_BAR:
        violations.append(f"Vacuum suction {P_suc_out/1e5:.2f}bar < {MIN_SUCTION_PRESSURE_BAR:.2f}bar")
    if (Pc / 1e5) > MAX_COND_PRESSURE_BAR:
        violations.append(f"P_cond {Pc/1e5:.1f}bar > {MAX_COND_PRESSURE_BAR:.1f}bar")
    if (T_dis - 273.15) > MAX_DISCHARGE_TEMP_C:
        violations.append(f"T_dis {T_dis-273.15:.1f}C > {MAX_DISCHARGE_TEMP_C:.1f}C")
    if dis_SH_K < MIN_DISCHARGE_SH_K:
        violations.append(f"Wet compression risk (dis SH {dis_SH_K:.1f}K < {MIN_DISCHARGE_SH_K:.1f}K)")
    if res_c.pinch_min_K < 0.015 or res_e.pinch_min_K < 0.015:
        violations.append(f"Pinch cross ({res_c.pinch_min_K:.2f}/{res_e.pinch_min_K:.2f}K)")

    # Exact expansion entropy generation (post-IHX!)
    try:
        if not use_ftvi:
            _, s_pre, _ = bk_cond.state_PH(Pc, h_evap_in)
            _, s_out, _ = bk_evap.state_PH(Pe, h_evap_in)
            s_gen_exp = m_dot_evap * max(0.0, s_out - s_pre)
        else:
            P_mid_f = max(Pe * 1.15, min(Pc * 0.85, design.ft_pr_bias * math.sqrt(Pe * Pc)))
            _, s_ft_m, _ = bk_cond.state_PH(P_mid_f, h_cond_out)
            _, s_ft_l, _ = bk_evap.state_PH(P_mid_f, h_evap_in)
            _, s_e2, _ = bk_evap.state_PH(Pe, h_evap_in)
            s_gen_exp = m_dot_cond * max(0.0, s_ft_m - s_cond_out) + m_dot_evap * max(0.0, s_e2 - s_ft_l)
    except Exception:
        s_gen_exp = 0.0

    # Real industrial Levelized Cost of Heat (LCOH, $/MWh_th)
    capex_comp = 9500.0 + 285.0 * (V_suc_m3_h ** 0.82)
    a_ihx_m2 = (0.4 + 2.2 * ihx_eff) if use_ihx else 0.0
    capex_hx = 450.0 * (tot_area + a_ihx_m2)
    capex_ft = 3500.0 if use_ftvi else 0.0
    capex_tot = 1.35 * (18000.0 + capex_comp + capex_hx + capex_ft)

    annual_capex_om = capex_tot * (CRF_15YR_7PCT + 0.025)
    annual_MWh_th = (duty.Q_cond_W / 1e6) * OPER_HOURS_YR
    annual_elec_usd = (W_tot / 1e6) * OPER_HOURS_YR * ELEC_PRICE_USD_MWH
    LCOH = (annual_capex_om + annual_elec_usd) / annual_MWh_th

    return HTHPEvaluation(
        label=design.label,
        duty_name=design.duty_name,
        feasible=(len(violations) == 0),
        violations=violations,
        COP_sys=COP_sys,
        COP_lorenz_ideal=cop_lorenz,
        eta_lorenz=eta_lorenz,
        LCOH_usd_MWh=LCOH,
        CAPEX_usd=capex_tot,
        VHC_kJ_m3=VHC_kJ_m3,
        V_suc_m3_h=V_suc_m3_h,
        W_comp_kW=W_comp_elec / 1e3,
        W_pumps_kW=W_pumps / 1e3,
        P_evap_bar=Pe / 1e5,
        P_evap_out_bar=P_suc_out / 1e5,
        P_cond_bar=Pc / 1e5,
        PR=Pc / P_suc_out,
        T_dis_C=T_dis - 273.15,
        dis_SH_K=dis_SH_K,
        glide_cond_K=sl_cond.glide,
        glide_evap_K=sl_evap.glide,
        sbg_cond=res_c.sbg_penalty_avg,
        sbg_evap=res_e.sbg_penalty_avg,
        dP_evap_kPa=res_e.dP_ref_Pa / 1e3,
        dP_cond_kPa=res_c.dP_ref_Pa / 1e3,
        pinch_cond_K=res_c.pinch_min_K,
        pinch_evap_K=res_e.pinch_min_K,
        gwp=bk_cond.gwp,
        is_pfas=is_pfas,
        I_comp_kW=(T_dead * s_gen_comp + 0.07 * W_comp_elec) / 1e3,
        I_cond_kW=T_dead * res_c.s_gen_W_K / 1e3,
        I_evap_kW=T_dead * res_e.s_gen_W_K / 1e3,
        I_exp_kW=T_dead * s_gen_exp / 1e3,
        I_ihx_ft_kW=T_dead * s_gen_ihx_ft / 1e3,
        I_pumps_kW=W_pumps / 1e3,
        z_evap=z_ev,
    )


def _infeasible_hthp(design: HTHPDesign, violations: List[str], is_pfas: bool, gwp: float) -> HTHPEvaluation:
    return HTHPEvaluation(
        label=design.label, duty_name=design.duty_name, feasible=False, violations=violations,
        COP_sys=0.0, COP_lorenz_ideal=1.0, eta_lorenz=0.0, LCOH_usd_MWh=999.0, CAPEX_usd=1e6,
        VHC_kJ_m3=0.0, V_suc_m3_h=9999.0, W_comp_kW=999.0, W_pumps_kW=0.0,
        P_evap_bar=0.0, P_evap_out_bar=0.0, P_cond_bar=0.0, PR=99.0, T_dis_C=999.0, dis_SH_K=-99.0,
        glide_cond_K=0.0, glide_evap_K=0.0, sbg_cond=1.0, sbg_evap=1.0,
        dP_evap_kPa=0.0, dP_cond_kPa=0.0, pinch_cond_K=0.0, pinch_evap_K=0.0,
        gwp=gwp, is_pfas=is_pfas,
        I_comp_kW=0.0, I_cond_kW=0.0, I_evap_kW=0.0, I_exp_kW=0.0, I_ihx_ft_kW=0.0, I_pumps_kW=0.0,
    )
