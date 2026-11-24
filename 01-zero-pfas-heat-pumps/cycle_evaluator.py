"""
Immutable, cheat-proof system-level reversible heat pump & tropical cooling
benchmark evaluator (6-bin dual-climate EN 14825 + ISO 5151 T1/T3 standard).

Supports four discrete thermodynamic cycle architectures:
  - TXV_BASIC   : Standard 1-stage vapor compression (+$0 BOM)
  - IHX_SUBCOOL : Suction-line internal heat exchanger (+$18 BOM)
  - FTVI_ECON   : Flash-tank vapor injection with exact zeotropic VLE split (+$55 BOM)
  - FTVI_IHX    : Combined FTVI + subcooling IHX (+$70 BOM)

Plus optional active/passive composition-shifting receiver (+$28 BOM) and
counterflow-preserving coil bridge headers (+$12/coil BOM).

All physical constraints are hard-gated:
  - Fixed coil face areas & aluminum mass caps (<= 3.2 kg indoor, <= 7.5 kg outdoor)
  - Explicit Kim-Bullard air-side fan power deducted from system COP
  - Mandatory Silver-Bell-Ghaly zeotropic HTC degradation & Muller-Steinhagen-Heck dP
  - Fixed compressor swept volume across all 6 bins (20..110 Hz); EN 14825 backup heat
  - Discharge temperature <= 115 C, High pressure <= 42 bar
  - IEC 60335-2-40 flammability charge cap (<= 500 g total flammable refrigerant)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from microchannel_hx import CoilGeometry, HXResult, solve_microchannel_hx
from thermo_backend import (
    FLUID_METADATA,
    MixtureBackend,
    compute_composition_shift_via_receiver,
)

# Locked physical coil envelope limits (standard 7 kW residential/light-commercial split)
INDOOR_W, INDOOR_H, INDOOR_MAX_DEPTH, INDOOR_MAX_AL_KG = 0.78, 0.30, 0.044, 3.20
OUTDOOR_W, OUTDOOR_H, OUTDOOR_MAX_DEPTH, OUTDOOR_MAX_AL_KG = 0.85, 0.62, 0.052, 7.50

# Locked safety & mechanical gates
MAX_T_DISCHARGE_C = 115.0
MAX_P_HIGH_BAR = 42.0
MAX_FLAMMABLE_CHARGE_G = 500.0   # Conservative safety margin below IEC 60335-2-40 Ed.7 988g limit
MIN_SUPERHEAT_K = 3.5
MIN_SUBCOOL_K = 2.0
COMP_F_MIN_HZ = 20.0
COMP_F_MAX_HZ = 110.0

# 6 Standardized Dual-Climate Operating Bins (ISO 5151 T1/T3 + EN 14825 Average Climate)
@dataclass(frozen=True)
class OperatingBin:
    name: str
    is_heating: bool
    T_out_C: float      # Outdoor air dry-bulb (C)
    T_in_C: float       # Indoor air dry-bulb (C)
    RH_in: float        # Indoor relative humidity (fraction)
    Q_load_W: float     # Required useful thermal load (W)
    weight: float       # Annual operating hours weight (sum = 1.0)


STANDARD_BINS: Tuple[OperatingBin, ...] = (
    # Tropical & Temperate Cooling (55% annual weight)
    OperatingBin("COOL_46C_T3",  False, 46.0, 29.0, 0.45, 6500.0, 0.14),  # Extreme heatwave survival bin
    OperatingBin("COOL_35C_T1",  False, 35.0, 27.0, 0.50, 7000.0, 0.18),  # ISO 5151 T1 rated cooling
    OperatingBin("COOL_29C_PL",  False, 29.0, 27.0, 0.50, 3850.0, 0.23),  # Dominant part-load cooling
    # EN 14825 Reversible Heating (45% annual weight)
    OperatingBin("HEAT_12C_PL",  True,  12.0, 20.0, 0.45, 2800.0, 0.15),  # Mild shoulder heating
    OperatingBin("HEAT_7C_RATED", True,  7.0, 20.0, 0.45, 4550.0, 0.17),  # Rated EN 14825 heating
    OperatingBin("HEAT_-7C_COLD", True, -7.0, 20.0, 0.45, 6200.0, 0.13),  # Cold-climate design heating
)


@dataclass(frozen=True)
class SystemDesign:
    """Complete discrete + continuous system specification evolved by autoresearch."""
    label: str
    # Layer 1: Discrete cycle architecture
    topology: str                        # TXV_BASIC | IHX_SUBCOOL | FTVI_ECON | FTVI_IHX
    ihx_eff: float                       # 0.0..0.75 (used in IHX_SUBCOOL, FTVI_IHX)
    ft_pr_bias: float                    # 0.80..1.25 multiplier on sqrt(P_evap*P_cond)
    comp_shift_enabled: bool             # Active composition-shifting receiver
    # Layer 2: Refrigerant fluid & composition
    components: Tuple[str, ...]
    mole_fractions: Tuple[float, ...]
    # Layer 3: Reversible microchannel HX geometry & circuitry
    indoor_coil: CoilGeometry
    outdoor_coil: CoilGeometry
    comp_disp_cc: float                  # Compressor swept displacement (cm3/rev)
    # Layer 4: Per-bin operating policy (6 bins)
    v_air_in_m3s: Tuple[float, ...]      # Indoor fan airflow per bin (m3/s)
    v_air_out_m3s: Tuple[float, ...]     # Outdoor fan airflow per bin (m3/s)
    superheat_K: Tuple[float, ...]       # Target evaporator superheat per bin (K)
    subcool_K: Tuple[float, ...]         # Target condenser subcooling per bin (K)
    liq_hold_frac: Tuple[float, ...]     # Receiver liquid hold fraction per bin (0..0.38)

    def incremental_bom_usd(self) -> float:
        """Transparent incremental hardware BOM vs standard 1-stage TXV split AC ($)."""
        bom = 0.0
        if self.topology == "IHX_SUBCOOL":
            bom += 18.0
        elif self.topology == "FTVI_ECON":
            bom += 55.0
        elif self.topology == "FTVI_IHX":
            bom += 70.0
        if self.comp_shift_enabled:
            bom += 28.0
        if self.indoor_coil.counterflow_bridge:
            bom += 12.0
        if self.outdoor_coil.counterflow_bridge:
            bom += 12.0
        # Extra aluminum cost delta above baseline ~6.5 kg total ($4/kg brazed Al)
        al_tot = self.indoor_coil.aluminum_mass + self.outdoor_coil.aluminum_mass
        bom += max(0.0, (al_tot - 6.5) * 4.5)
        return bom


@dataclass
class BinEvaluation:
    bin_name: str
    feasible: bool
    reason: str
    Q_delivered_W: float
    Q_load_W: float
    W_comp_W: float
    W_fan_in_W: float
    W_fan_out_W: float
    W_backup_W: float
    COP_sys: float
    COP_carnot_ref: float
    P_evap_bar: float
    P_cond_bar: float
    T_dis_C: float
    f_comp_Hz: float
    charge_g: float
    glide_evap_K: float
    glide_cond_K: float
    sbg_evap: float
    sbg_cond: float
    r_air_evap: float
    r_air_cond: float
    dP_evap_kPa: float
    dP_cond_kPa: float
    z_circ: Tuple[float, ...]
    # Gouy-Stodola Exergy Destruction Breakdown (W)
    I_comp_W: float
    I_cond_W: float
    I_evap_W: float
    I_exp_W: float
    I_ihx_ft_W: float
    I_fans_W: float


@dataclass
class SystemEvaluation:
    label: str
    feasible: bool
    violations: List[str]
    ASPF: float                  # Blended Annual Seasonal Performance Factor (primary metric)
    SEER_cool: float             # Cooling-only weighted COP
    SCOP_heat: float             # Heating-only weighted COP
    COP_46C_T3: float            # Extreme 46 C heatwave COP
    COP_m7C_cold: float          # Cold -7 C heating COP (incl. backup heat)
    max_charge_g: float          # Worst-case total charge across bins (g)
    flammable_charge_g: float    # Flammable refrigerant mass (g)
    max_T_dis_C: float
    max_P_high_bar: float
    max_glide_K: float
    gwp: float
    al_mass_in_kg: float
    al_mass_out_kg: float
    inc_bom_usd: float
    bins: List[BinEvaluation] = field(default_factory=list)


def _compressor_efficiencies(PR: float, f_hz: float) -> Tuple[float, float, float]:
    """
    Realistic variable-speed inverter scroll/twin-rotary compressor efficiency model
    calibrated to published R-290/R-32 OEM maps (GMCC / Copeland).
    Returns (eta_isentropic_mech, eta_volumetric, eta_motor_inverter).
    """
    PR_c = max(1.25, min(14.0, PR))
    # Peak isentropic efficiency ~0.77 around PR=2.8-3.4, dropping at extreme PR or speed
    pr_pen = 0.024 * ((PR_c - 3.0) ** 2) / PR_c
    f_norm = f_hz / 65.0
    f_pen = 0.045 * ((f_norm - 0.85) ** 2)
    eta_is = max(0.42, min(0.79, 0.78 - pr_pen - f_pen))
    # Clearance volumetric efficiency
    eta_v = max(0.55, min(0.95, 0.96 - 0.030 * ((PR_c ** 0.85) - 1.0)))
    # Permanent-magnet / inverter drive efficiency (drops lightly at very low speed)
    eta_em = max(0.84, min(0.93, 0.92 - 0.035 * max(0.0, 0.55 - f_norm)))
    return eta_is, eta_v, eta_em


def _solve_single_bin(
    design: SystemDesign,
    bin_idx: int,
    op_bin: OperatingBin,
    n_seg: int = 12,
    htc_mult: float = 1.0,
    sbg_mult: float = 1.0,
    dp_mult: float = 1.0,
) -> BinEvaluation:
    """Solve coupled compressor + condenser + evaporator + IHX/FTVI balance for one bin."""
    T_out_K = op_bin.T_out_C + 273.15
    T_in_K = op_bin.T_in_C + 273.15
    T_dead = T_out_K  # Ambient dead-state temperature for Gouy-Stodola exergy audit

    # Active composition shift via receiver (if enabled)
    hold_f = design.liq_hold_frac[bin_idx] if design.comp_shift_enabled else 0.0
    if design.comp_shift_enabled and len(design.components) > 1 and hold_f > 1e-4:
        # Receiver sits near outdoor/condenser-subcooled temperature (~T_out + 12K in heat, T_in + 8K in cool)
        T_vessel = (T_out_K + 12.0) if op_bin.is_heating else (T_in_K + 6.0)
        z_circ = compute_composition_shift_via_receiver(
            design.components, design.mole_fractions, hold_f, T_vessel
        )
    else:
        z_circ = design.mole_fractions

    try:
        bk_cond = MixtureBackend(design.components, z_circ)
    except Exception as e:
        return _infeasible_bin(op_bin, f"EOS init failed: {e}", z_circ)

    # Identify which physical coil is condenser vs evaporator
    if op_bin.is_heating:
        cond_coil, evap_coil = design.indoor_coil, design.outdoor_coil
        T_air_cond, RH_cond, v_air_cond = T_in_K, op_bin.RH_in, design.v_air_in_m3s[bin_idx]
        T_air_evap, RH_evap, v_air_evap = T_out_K, 0.80, design.v_air_out_m3s[bin_idx]
    else:
        cond_coil, evap_coil = design.outdoor_coil, design.indoor_coil
        T_air_cond, RH_cond, v_air_cond = T_out_K, 0.45, design.v_air_out_m3s[bin_idx]
        T_air_evap, RH_evap, v_air_evap = T_in_K, op_bin.RH_in, design.v_air_in_m3s[bin_idx]

    sh_K = max(MIN_SUPERHEAT_K, design.superheat_K[bin_idx])
    sc_K = max(MIN_SUBCOOL_K, design.subcool_K[bin_idx])
    V_disp_m3 = design.comp_disp_cc * 1e-6

    # Initial temperature guesses for saturation midpoints
    T_evap_mid = T_air_evap - 11.0
    T_cond_mid = T_air_cond + 12.0

    try:
        P_evap = bk_cond.p_for_t_mid(T_evap_mid)
        P_cond = bk_cond.p_for_t_mid(T_cond_mid)
    except Exception as e:
        return _infeasible_bin(op_bin, f"Initial P bracket failed: {e}", z_circ)

    P_evap = max(0.8e5, min(28e5, P_evap))
    P_cond = max(P_evap * 1.45, min(MAX_P_HIGH_BAR * 1e5 * 1.05, 0.92 * bk_cond.Pc, P_cond))

    f_hz = 60.0
    use_ftvi = design.topology in ("FTVI_ECON", "FTVI_IHX")
    use_ihx = design.topology in ("IHX_SUBCOOL", "FTVI_IHX")
    ihx_eff = max(0.0, min(0.75, design.ihx_eff)) if use_ihx else 0.0

    # Exact 2D bounded Newton-Raphson solver with macro finite-difference Jacobian
    # (step = 0.01 in ln P, strictly larger than isobar cache quantization -> exact smooth J)
    cache_eval = {"state": None}
    ln_pe_lo = math.log(0.7e5)
    ln_pe_hi = math.log(28e5)
    ln_pc_hi = math.log(min(MAX_P_HIGH_BAR * 1e5 * 1.12, 0.92 * bk_cond.Pc))

    def eval_residual(x_vec: np.ndarray) -> np.ndarray:
        ln_pe = float(np.clip(x_vec[0], ln_pe_lo, ln_pe_hi))
        ln_pc = float(np.clip(x_vec[1], ln_pe + math.log(1.30), ln_pc_hi))
        Pe = math.exp(ln_pe)
        Pc = math.exp(ln_pc)

        sl_cond = bk_cond.get_isobar(Pc, nq=9)
        T_cond_out_target = sl_cond.Tb - sc_K
        h_cond_out_target, s_cond_out, rho_cond_out = bk_cond.liquid_state_PT(Pc, T_cond_out_target)

        if use_ftvi:
            P_mid = design.ft_pr_bias * math.sqrt(Pe * Pc)
            P_mid = max(Pe * 1.15, min(Pc * 0.85, P_mid))
            T_ft, q_ft_mass, x_liq_ft, y_vap_ft, h_ft_liq, h_ft_vap = bk_cond.flash_phase_split(P_mid, h_cond_out_target)
            q_ft_mass = max(0.02, min(0.45, q_ft_mass))
            if len(design.components) > 1:
                bk_evap = MixtureBackend(design.components, tuple(float(v) for v in x_liq_ft))
            else:
                bk_evap = bk_cond
            h_pre_eev2 = h_ft_liq
            T_pre_eev2 = T_ft
        else:
            P_mid = Pe
            q_ft_mass = 0.0
            bk_evap = bk_cond
            h_pre_eev2 = h_cond_out_target
            T_pre_eev2 = T_cond_out_target
            h_ft_vap = 0.0
            y_vap_ft = np.array(z_circ)

        sl_evap = bk_evap.get_isobar(Pe, nq=9)
        T_evap_out_target = sl_evap.Td + sh_K
        h_evap_out_target, s_evap_out, rho_evap_out = bk_evap.vapor_state_PT(Pe, T_evap_out_target)

        if use_ihx and ihx_eff > 1e-3 and T_pre_eev2 > T_evap_out_target + 2.0:
            C_vap = sl_evap.cp_v
            C_liq = sl_cond.cp_l
            C_min = min(C_vap, C_liq)
            dh_ihx = ihx_eff * C_min * (T_pre_eev2 - T_evap_out_target)
            dh_ihx = min(dh_ihx, 0.85 * C_liq * max(0.0, T_pre_eev2 - T_evap_out_target - 1.5))
            h_suc = h_evap_out_target + dh_ihx
            T_suc = T_evap_out_target + dh_ihx / max(sl_evap.cp_v, 500.0)
            h_evap_in = h_pre_eev2 - dh_ihx
            _, s_suc, rho_suc = bk_evap.vapor_state_PT(Pe, T_suc)
        else:
            dh_ihx = 0.0
            h_suc = h_evap_out_target
            T_suc = T_evap_out_target
            s_suc = s_evap_out
            rho_suc = rho_evap_out
            h_evap_in = h_pre_eev2

        PR_tot = Pc / Pe
        f_loc = 60.0
        if not use_ftvi:
            h_dis_s, _ = bk_cond.state_PS(Pc, s_suc)
            for _ in range(3):
                eta_is, eta_v, eta_em = _compressor_efficiencies(PR_tot, f_loc)
                w_spec_shaft = max(1e3, (h_dis_s - h_suc) / eta_is)
                h_cond_in = h_suc + w_spec_shaft
                dh_evap_id = max(20e3, h_evap_out_target - h_evap_in)
                dh_cond_id = max(20e3, h_cond_in - h_cond_out_target)
                m_req = op_bin.Q_load_W / (dh_cond_id if op_bin.is_heating else dh_evap_id)
                f_loc = max(COMP_F_MIN_HZ, min(COMP_F_MAX_HZ, m_req / max(rho_suc * V_disp_m3 * eta_v, 1e-9)))
            eta_is, eta_v, eta_em = _compressor_efficiencies(PR_tot, f_loc)
            w_spec_shaft = max(1e3, (h_dis_s - h_suc) / eta_is)
            h_cond_in = h_suc + w_spec_shaft
            T_dis, s_dis, _ = bk_cond.state_PH(Pc, h_cond_in)
            m_dot_evap = rho_suc * V_disp_m3 * f_loc * eta_v
            m_dot_cond = m_dot_evap
            W_comp_elec = m_dot_evap * w_spec_shaft / eta_em
            s_gen_comp = m_dot_evap * max(0.0, s_dis - s_suc)
            s_gen_ihx_ft = 0.0
            if dh_ihx > 10.0:
                T_liq_avg = T_pre_eev2 - 0.5 * dh_ihx / max(sl_cond.cp_l, 500.0)
                T_vap_avg = 0.5 * (T_evap_out_target + T_suc)
                s_gen_ihx_ft += m_dot_evap * dh_ihx * max(0.0, 1.0 / T_vap_avg - 1.0 / T_liq_avg)
        else:
            PR1 = P_mid / Pe
            PR2 = Pc / P_mid
            h_mid_s, _ = bk_evap.state_PS(P_mid, s_suc)
            inj_ratio = q_ft_mass / max(1.0 - q_ft_mass, 0.2)
            for _ in range(3):
                eta_is1, eta_v, eta_em = _compressor_efficiencies(PR1, f_loc)
                eta_is2, _, _ = _compressor_efficiencies(PR2, f_loc)
                w1_shaft = max(500.0, (h_mid_s - h_suc) / eta_is1)
                h_mid_1 = h_suc + w1_shaft
                h_mid_mix = (h_mid_1 + inj_ratio * h_ft_vap) / (1.0 + inj_ratio)
                _, s_mid_mix, _ = bk_cond.state_PH(P_mid, h_mid_mix)
                h_dis_s, _ = bk_cond.state_PS(Pc, s_mid_mix)
                w2_shaft = max(500.0, (h_dis_s - h_mid_mix) / eta_is2)
                h_cond_in = h_mid_mix + w2_shaft
                dh_evap_id = max(20e3, h_evap_out_target - h_evap_in)
                dh_cond_id = max(20e3, h_cond_in - h_cond_out_target)
                m_ev_req = (op_bin.Q_load_W / dh_cond_id) / (1.0 + inj_ratio) if op_bin.is_heating else (op_bin.Q_load_W / dh_evap_id)
                f_loc = max(COMP_F_MIN_HZ, min(COMP_F_MAX_HZ, m_ev_req / max(rho_suc * V_disp_m3 * eta_v, 1e-9)))
            eta_is1, eta_v, eta_em = _compressor_efficiencies(PR1, f_loc)
            eta_is2, _, _ = _compressor_efficiencies(PR2, f_loc)
            w1_shaft = max(500.0, (h_mid_s - h_suc) / eta_is1)
            h_mid_1 = h_suc + w1_shaft
            _, s_mid_1, _ = bk_evap.state_PH(P_mid, h_mid_1)
            h_mid_mix = (h_mid_1 + inj_ratio * h_ft_vap) / (1.0 + inj_ratio)
            T_mid_mix, s_mid_mix, _ = bk_cond.state_PH(P_mid, h_mid_mix)
            try:
                bk_vap_ft = MixtureBackend(design.components, tuple(float(v) for v in y_vap_ft)) if len(design.components) > 1 else bk_cond
                _, s_ft_vap, _ = bk_vap_ft.state_PH(P_mid, h_ft_vap)
            except Exception:
                s_ft_vap = s_mid_mix
            h_dis_s, _ = bk_cond.state_PS(Pc, s_mid_mix)
            w2_shaft = max(500.0, (h_dis_s - h_mid_mix) / eta_is2)
            h_cond_in = h_mid_mix + w2_shaft
            T_dis, s_dis, _ = bk_cond.state_PH(Pc, h_cond_in)
            m_dot_evap = rho_suc * V_disp_m3 * f_loc * eta_v
            m_dot_cond = m_dot_evap * (1.0 + inj_ratio)
            W_comp_elec = (m_dot_evap * w1_shaft + m_dot_cond * w2_shaft) / eta_em
            s_gen_comp = m_dot_evap * max(0.0, s_mid_1 - s_suc) + m_dot_cond * max(0.0, s_dis - s_mid_mix)
            m_dot_inj = m_dot_cond - m_dot_evap
            s_gen_ihx_ft = max(0.0, m_dot_cond * s_mid_mix - m_dot_evap * s_mid_1 - m_dot_inj * s_ft_vap)
            if dh_ihx > 10.0:
                T_liq_avg = T_pre_eev2 - 0.5 * dh_ihx / max(sl_cond.cp_l, 500.0)
                T_vap_avg = 0.5 * (T_evap_out_target + T_suc)
                s_gen_ihx_ft += m_dot_evap * dh_ihx * max(0.0, 1.0 / T_vap_avg - 1.0 / T_liq_avg)

        res_cond = solve_microchannel_hx(
            coil=cond_coil, backend=bk_cond, is_condenser=True,
            is_heating_mode=op_bin.is_heating,
            P_in=Pc, h_in=h_cond_in, m_dot_ref=m_dot_cond,
            T_air_in=T_air_cond, RH_air_in=RH_cond, v_dot_air=v_air_cond,
            n_seg=n_seg, htc_mult=htc_mult, sbg_mult=sbg_mult, dp_mult=dp_mult,
        )
        res_evap = solve_microchannel_hx(
            coil=evap_coil, backend=bk_evap, is_condenser=False,
            is_heating_mode=op_bin.is_heating,
            P_in=Pe, h_in=h_evap_in, m_dot_ref=m_dot_evap,
            T_air_in=T_air_evap, RH_air_in=RH_evap, v_dot_air=v_air_evap,
            n_seg=n_seg, htc_mult=htc_mult, sbg_mult=sbg_mult, dp_mult=dp_mult,
        )

        Q_cond_req = m_dot_cond * (h_cond_in - h_cond_out_target)
        Q_evap_req = m_dot_evap * (h_evap_out_target - h_evap_in)
        err_evap = (Q_evap_req - res_evap.Q) / max(Q_evap_req, 500.0)
        err_cond = (Q_cond_req - res_cond.Q) / max(Q_cond_req, 500.0)

        cache_eval["state"] = (
            sl_cond, sl_evap, res_cond, res_evap,
            W_comp_elec, T_dis, f_loc, m_dot_evap, m_dot_cond,
            s_gen_comp, s_gen_ihx_ft, s_cond_out, h_evap_in, h_cond_out_target,
            bk_evap, Pe, Pc, err_evap, err_cond,
        )
        return np.array([err_evap, err_cond], dtype=float)

    try:
        x = np.array([math.log(P_evap), math.log(P_cond)], dtype=float)
        r = eval_residual(x)
        converged = False
        fd_h = 0.012
        for it in range(12):
            norm_r = max(abs(r[0]), abs(r[1]))
            if norm_r < 0.005:
                converged = True
                break
            # Compute explicit 2x2 macro Jacobian
            r1 = eval_residual(np.array([x[0] + fd_h, x[1]]))
            r2 = eval_residual(np.array([x[0], x[1] + fd_h]))
            J = np.column_stack([(r1 - r) / fd_h, (r2 - r) / fd_h])
            try:
                dx = -np.linalg.solve(J + 1e-4 * np.eye(2), r)
            except np.linalg.LinAlgError:
                dx = np.array([-0.25 * r[0], 0.25 * r[1]])
            # Trust-region cap on step size in ln(P)
            max_step = max(abs(dx[0]), abs(dx[1]))
            if max_step > 0.25:
                dx *= 0.25 / max_step
            # Backtracking line search
            accepted = False
            for alpha in (1.0, 0.5, 0.25):
                x_try = np.array([
                    float(np.clip(x[0] + alpha * dx[0], ln_pe_lo, ln_pe_hi)),
                    float(np.clip(x[1] + alpha * dx[1], x[0] + alpha * dx[0] + math.log(1.30), ln_pc_hi)),
                ])
                x_try[1] = max(x_try[0] + math.log(1.30), min(ln_pc_hi, x_try[1]))
                r_try = eval_residual(x_try)
                if max(abs(r_try[0]), abs(r_try[1])) < norm_r * 0.98 or alpha == 0.25:
                    x, r = x_try, r_try
                    accepted = True
                    break
            if not accepted:
                break
        if not converged:
            # Final re-evaluation at best x to ensure cache consistency
            r = eval_residual(x)
            converged = max(abs(r[0]), abs(r[1])) < 0.015
    except Exception as e:
        return _infeasible_bin(op_bin, f"Cycle solver error: {e}", z_circ)

    if cache_eval["state"] is None or not converged:
        return _infeasible_bin(op_bin, "Cycle-HX balance did not converge", z_circ)

    (
        sl_cond, sl_evap, res_cond, res_evap,
        W_comp_elec, T_dis, f_hz, m_dot_evap, m_dot_cond,
        s_gen_comp, s_gen_ihx_ft, s_cond_out, h_evap_in, h_cond_out_target,
        bk_evap, P_evap, P_cond, _, _,
    ) = cache_eval["state"]

    # Account for suction-side pressure drop penalty on compressor work!
    # Compressor actually pulls from P_suc = P_evap - dP_evap, raising actual shaft power
    P_suc_actual = max(0.5e5, P_evap - res_evap.dP_ref)
    dp_comp_factor = math.log(P_cond / P_suc_actual) / max(0.15, math.log(P_cond / P_evap))
    W_comp_elec *= max(1.0, min(1.35, dp_comp_factor))

    # Delivered useful heat & EN 14825 backup resistance penalty
    Q_delivered = res_cond.Q if op_bin.is_heating else res_evap.Q
    Q_deficit = max(0.0, op_bin.Q_load_W - Q_delivered)
    if op_bin.is_heating:
        # EN 14825 exact rule: electric resistance backup covers shortfall at COP=1.0
        W_backup = Q_deficit
    else:
        # Cooling mode: cannot use resistance heater to cool! Require >= 92% of load!
        if Q_delivered < 0.92 * op_bin.Q_load_W:
            return _infeasible_bin(
                op_bin,
                f"Cooling capacity shortfall ({Q_delivered:.0f}W < 92% of {op_bin.Q_load_W:.0f}W)",
                z_circ,
            )
        W_backup = 0.0

    # Total system refrigerant charge (indoor + outdoor + 6m liquid/suction lines + receiver/FT)
    v_liq_line = math.pi * (0.0032**2) * 5.0    # 6.4 mm ID, 5 m liquid line
    v_suc_line = math.pi * (0.0065**2) * 5.0    # 13 mm ID, 5 m suction line
    m_lines_kg = v_liq_line * sl_cond.rho_l + v_suc_line * sl_evap.rho_v
    m_ft_kg = (0.00018 * (0.35 * sl_cond.rho_l + 0.65 * sl_evap.rho_v)) if use_ftvi else 0.0
    m_base_kg = res_cond.charge_kg + res_evap.charge_kg + m_lines_kg + m_ft_kg
    # Active composition-shift receiver adds trapped liquid mass proportional to hold_f!
    if design.comp_shift_enabled and hold_f > 1e-4:
        m_receiver_kg = m_base_kg * (hold_f / max(0.45, 1.0 - hold_f))
    else:
        m_receiver_kg = 0.0
    total_charge_g = 1000.0 * (m_base_kg + m_receiver_kg)

    # Fan power mapping
    W_fan_in = res_cond.W_fan if op_bin.is_heating else res_evap.W_fan
    W_fan_out = res_evap.W_fan if op_bin.is_heating else res_cond.W_fan

    W_tot = W_comp_elec + W_fan_in + W_fan_out + W_backup
    Q_eff = op_bin.Q_load_W if op_bin.is_heating else Q_delivered
    COP_sys = Q_eff / max(W_tot, 10.0)

    # Ideal Carnot reference between outdoor & indoor air inlet temperatures
    if op_bin.is_heating:
        COP_carnot = T_in_K / max(1.0, T_in_K - T_out_K)
    else:
        COP_carnot = T_in_K / max(1.0, T_out_K - T_in_K)

    # Exact expansion valve entropy generation (1-stage EEV or 2-stage EEV1 + EEV2)
    try:
        if not use_ftvi:
            # Account for subcooling inside IHX before expansion valve!
            _, s_pre_eev, _ = bk_cond.state_PH(P_cond, h_evap_in)
            _, s_exp_out, _ = bk_evap.state_PH(P_evap, h_evap_in)
            s_gen_exp = m_dot_evap * max(0.0, s_exp_out - s_pre_eev)
        else:
            P_mid_f = max(P_evap * 1.15, min(P_cond * 0.85, design.ft_pr_bias * math.sqrt(P_evap * P_cond)))
            _, s_ft_mix, _ = bk_cond.state_PH(P_mid_f, h_cond_out_target)
            _, s_ft_l, _ = bk_evap.state_PH(P_mid_f, h_evap_in)
            _, s_exp2_out, _ = bk_evap.state_PH(P_evap, h_evap_in)
            s_gen_exp = m_dot_cond * max(0.0, s_ft_mix - s_cond_out) + m_dot_evap * max(0.0, s_exp2_out - s_ft_l)
    except Exception:
        s_gen_exp = 0.0

    # Gouy-Stodola exergy destruction (W)
    I_comp = T_dead * s_gen_comp + W_comp_elec * 0.09  # include motor/inverter loss
    I_cond = T_dead * res_cond.s_gen
    I_evap = T_dead * res_evap.s_gen
    I_exp = T_dead * s_gen_exp
    I_ihx_ft = T_dead * s_gen_ihx_ft
    I_fans = W_fan_in + W_fan_out

    return BinEvaluation(
        bin_name=op_bin.name,
        feasible=True,
        reason="OK",
        Q_delivered_W=Q_delivered,
        Q_load_W=op_bin.Q_load_W,
        W_comp_W=W_comp_elec,
        W_fan_in_W=W_fan_in,
        W_fan_out_W=W_fan_out,
        W_backup_W=W_backup,
        COP_sys=COP_sys,
        COP_carnot_ref=COP_carnot,
        P_evap_bar=P_evap / 1e5,
        P_cond_bar=P_cond / 1e5,
        T_dis_C=T_dis - 273.15,
        f_comp_Hz=f_hz,
        charge_g=total_charge_g,
        glide_evap_K=sl_evap.glide,
        glide_cond_K=sl_cond.glide,
        sbg_evap=res_evap.sbg_penalty_avg,
        sbg_cond=res_cond.sbg_penalty_avg,
        r_air_evap=res_evap.r_air_share,
        r_air_cond=res_cond.r_air_share,
        dP_evap_kPa=res_evap.dP_ref / 1e3,
        dP_cond_kPa=res_cond.dP_ref / 1e3,
        z_circ=z_circ,
        I_comp_W=I_comp,
        I_cond_W=I_cond,
        I_evap_W=I_evap,
        I_exp_W=I_exp,
        I_ihx_ft_W=I_ihx_ft,
        I_fans_W=I_fans,
    )


def _infeasible_bin(op_bin: OperatingBin, reason: str, z_circ: Tuple[float, ...]) -> BinEvaluation:
    return BinEvaluation(
        bin_name=op_bin.name, feasible=False, reason=reason,
        Q_delivered_W=0.0, Q_load_W=op_bin.Q_load_W,
        W_comp_W=1e6, W_fan_in_W=0.0, W_fan_out_W=0.0, W_backup_W=0.0,
        COP_sys=0.0, COP_carnot_ref=1.0,
        P_evap_bar=0.0, P_cond_bar=0.0, T_dis_C=999.0, f_comp_Hz=0.0,
        charge_g=9999.0, glide_evap_K=0.0, glide_cond_K=0.0,
        sbg_evap=1.0, sbg_cond=1.0, r_air_evap=0.0, r_air_cond=0.0,
        dP_evap_kPa=0.0, dP_cond_kPa=0.0, z_circ=z_circ,
        I_comp_W=0.0, I_cond_W=0.0, I_evap_W=0.0, I_exp_W=0.0, I_ihx_ft_W=0.0, I_fans_W=0.0,
    )


def evaluate_system(
    design: SystemDesign,
    enforce_natural_gwp: bool = True,
    max_glide_K: float = 12.0,
    n_seg: int = 12,
    htc_mult: float = 1.0,
    sbg_mult: float = 1.0,
    dp_mult: float = 1.0,
) -> SystemEvaluation:
    """
    Immutable full 6-bin dual-climate benchmark evaluator with hard physical gates.
    """
    violations: List[str] = []

    # Tier 0: Hard hardware & safety gates
    ic, oc = design.indoor_coil, design.outdoor_coil
    if ic.width > INDOOR_W + 1e-6 or ic.height > INDOOR_H + 1e-6:
        violations.append("Indoor coil exceeds face envelope")
    if oc.width > OUTDOOR_W + 1e-6 or oc.height > OUTDOOR_H + 1e-6:
        violations.append("Outdoor coil exceeds face envelope")
    if ic.total_depth > INDOOR_MAX_DEPTH + 1e-6:
        violations.append(f"Indoor depth {ic.total_depth*1e3:.1f}mm > {INDOOR_MAX_DEPTH*1e3:.1f}mm")
    if oc.total_depth > OUTDOOR_MAX_DEPTH + 1e-6:
        violations.append(f"Outdoor depth {oc.total_depth*1e3:.1f}mm > {OUTDOOR_MAX_DEPTH*1e3:.1f}mm")
    if ic.aluminum_mass > INDOOR_MAX_AL_KG:
        violations.append(f"Indoor Al {ic.aluminum_mass:.2f}kg > {INDOOR_MAX_AL_KG:.2f}kg")
    if oc.aluminum_mass > OUTDOOR_MAX_AL_KG:
        violations.append(f"Outdoor Al {oc.aluminum_mass:.2f}kg > {OUTDOOR_MAX_AL_KG:.2f}kg")
    if ic.port_dh < 0.00065 or oc.port_dh < 0.00065:
        violations.append("Port Dh < 0.65mm (oil plugging / maldistribution risk)")

    try:
        bk0 = MixtureBackend(design.components, design.mole_fractions)
    except Exception as e:
        return SystemEvaluation(
            label=design.label, feasible=False, violations=[f"Invalid fluid/EOS: {e}"],
            ASPF=0.0, SEER_cool=0.0, SCOP_heat=0.0, COP_46C_T3=0.0, COP_m7C_cold=0.0,
            max_charge_g=9999.0, flammable_charge_g=9999.0,
            max_T_dis_C=999.0, max_P_high_bar=999.0, max_glide_K=99.0,
            gwp=9999.0, al_mass_in_kg=ic.aluminum_mass, al_mass_out_kg=oc.aluminum_mass,
            inc_bom_usd=design.incremental_bom_usd(), bins=[],
        )
    if enforce_natural_gwp and bk0.gwp > 5.0:
        violations.append(f"GWP {bk0.gwp:.1f} > 5.0")

    bin_results: List[BinEvaluation] = []
    for idx, op_bin in enumerate(STANDARD_BINS):
        br = _solve_single_bin(
            design, idx, op_bin, n_seg=n_seg,
            htc_mult=htc_mult, sbg_mult=sbg_mult, dp_mult=dp_mult,
        )
        bin_results.append(br)
        if not br.feasible:
            violations.append(f"{op_bin.name}: {br.reason}")

    if any(not b.feasible for b in bin_results):
        return SystemEvaluation(
            label=design.label, feasible=False, violations=violations,
            ASPF=0.0, SEER_cool=0.0, SCOP_heat=0.0, COP_46C_T3=0.0, COP_m7C_cold=0.0,
            max_charge_g=9999.0, flammable_charge_g=9999.0,
            max_T_dis_C=999.0, max_P_high_bar=999.0, max_glide_K=99.0,
            gwp=bk0.gwp, al_mass_in_kg=ic.aluminum_mass, al_mass_out_kg=oc.aluminum_mass,
            inc_bom_usd=design.incremental_bom_usd(), bins=bin_results,
        )

    max_charge_g = max(b.charge_g for b in bin_results)
    flam_charge_g = max_charge_g * bk0.flam_mass_frac
    max_T_dis_C = max(b.T_dis_C for b in bin_results)
    max_P_high_bar = max(b.P_cond_bar for b in bin_results)
    actual_max_glide_K = max(max(b.glide_evap_K, b.glide_cond_K) for b in bin_results)

    # IEC 60335-2-40 class-specific flammability charge cap:
    # A3 natural hydrocarbons: strict 500 g cap (conservative vs 988 g Ed.7 indoor limit)
    # A2L (R-32 baseline): 1800 g; A1 (R-410A/CO2): unrestricted
    has_a3 = any(FLUID_METADATA[c]["safety"] == "A3" and z > 1e-4 for c, z in zip(design.components, design.mole_fractions))
    has_a2l = any(FLUID_METADATA[c]["safety"] == "A2L" and z > 1e-4 for c, z in zip(design.components, design.mole_fractions))
    if has_a3 and flam_charge_g > MAX_FLAMMABLE_CHARGE_G:
        violations.append(f"A3 flammable charge {flam_charge_g:.0f}g > {MAX_FLAMMABLE_CHARGE_G:.0f}g")
    elif (not has_a3) and has_a2l and flam_charge_g > 1800.0:
        violations.append(f"A2L flammable charge {flam_charge_g:.0f}g > 1800g")
    if max_T_dis_C > MAX_T_DISCHARGE_C:
        violations.append(f"Max T_dis {max_T_dis_C:.1f}C > {MAX_T_DISCHARGE_C:.1f}C")
    if max_P_high_bar > MAX_P_HIGH_BAR:
        violations.append(f"Max P_high {max_P_high_bar:.1f}bar > {MAX_P_HIGH_BAR:.1f}bar")
    if actual_max_glide_K > max_glide_K:
        violations.append(f"Max glide {actual_max_glide_K:.1f}K > {max_glide_K:.1f}K limit")

    # Compute exact load-weighted Seasonal Performance Factors
    Q_all, W_all = 0.0, 0.0
    Q_cool, W_cool = 0.0, 0.0
    Q_heat, W_heat = 0.0, 0.0
    for op_bin, b in zip(STANDARD_BINS, bin_results):
        W_b = b.W_comp_W + b.W_fan_in_W + b.W_fan_out_W + b.W_backup_W
        Q_b = op_bin.Q_load_W if op_bin.is_heating else b.Q_delivered_W
        Q_all += op_bin.weight * Q_b
        W_all += op_bin.weight * W_b
        if op_bin.is_heating:
            Q_heat += op_bin.weight * Q_b
            W_heat += op_bin.weight * W_b
        else:
            Q_cool += op_bin.weight * Q_b
            W_cool += op_bin.weight * W_b

    ASPF = Q_all / max(W_all, 1.0)
    SEER_cool = Q_cool / max(W_cool, 1.0)
    SCOP_heat = Q_heat / max(W_heat, 1.0)

    return SystemEvaluation(
        label=design.label,
        feasible=(len(violations) == 0),
        violations=violations,
        ASPF=ASPF,
        SEER_cool=SEER_cool,
        SCOP_heat=SCOP_heat,
        COP_46C_T3=bin_results[0].COP_sys,
        COP_m7C_cold=bin_results[5].COP_sys,
        max_charge_g=max_charge_g,
        flammable_charge_g=flam_charge_g,
        max_T_dis_C=max_T_dis_C,
        max_P_high_bar=max_P_high_bar,
        max_glide_K=actual_max_glide_K,
        gwp=bk0.gwp,
        al_mass_in_kg=ic.aluminum_mass,
        al_mass_out_kg=oc.aluminum_mass,
        inc_bom_usd=design.incremental_bom_usd(),
        bins=bin_results,
    )


def max_glide_K_guard(design: SystemDesign, limit_K: float) -> float:
    return limit_K
