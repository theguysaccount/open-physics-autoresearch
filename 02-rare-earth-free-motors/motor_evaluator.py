"""
motor_evaluator.py — Locked Multi-Point WLTC Traction Drive-Cycle Evaluator
Coupling Exact Nonlinear 2D Maxwell FEA, 2D Centrifugal Elasticity FEA,
Worst-Case Cold (-30C Ferrite) / Hot (+140C NdFeB) Demagnetization,
Harmonic Current Injection, SiC Inverter Losses & kVA Sizing, and Total System Cost.
"""

from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Dict, List, Tuple, Optional
import math
import numpy as np
from scipy.interpolate import RectBivariateSpline

from materials import (
    STEEL_RHO, STEEL_COST_PER_KG, CU_DENSITY, CU_COST_PER_KG,
    MAGNET_DB, V_DC_BUS, V_PH_PEAK_MAX,
    INVERTER_USD_PER_KVA, INVERTER_FIXED_USD, INVERTER_FSW_HZ,
)
from motor_geometry import (
    L_STK, K_FILL_MAX, N_THETA_GAP,
    MotorDesign, MotorMeshBundle, build_motor_mesh, validate_design_geometry,
)
from maxwell_fea import MotorFEASolver, solve_centrifugal_stress, GAMMA_E0_RAD


# ---------------------------------------------------------------------------
# Locked 5-Point EV Traction Drive-Cycle Benchmark (80 kW peak / 40 kW cont, 400 V)
# ---------------------------------------------------------------------------
WLTC_OPS = [
    # (name, rpm, target_shaft_torque_Nm, weight)
    ("OP1_URBAN",    2000.0, 40.0, 0.30),   # 8.38 kW urban cruise
    ("OP2_SUBURBAN", 4000.0, 75.0, 0.30),   # 31.42 kW suburban acceleration/cruise
    ("OP3_HIGHWAY",  7500.0, 45.0, 0.25),   # 35.34 kW highway cruise
    ("OP4_HIGHSPD", 11500.0, 30.0, 0.15),   # 36.13 kW deep field-weakening cruise
]
PEAK_OP_RPM = 2500.0
PEAK_TORQUE_MIN_NM = 185.0    # >= 185 Nm @ 2500 rpm (48.4 kW -> ~77-85 kW @ 4000-4400 rpm!)
J_RMS_CONT_MAX = 12.5         # A/mm^2 continuous thermal cap at OP1..OP4
J_RMS_PEAK_MAX = 24.0         # A/mm^2 10-sec peak launch thermal cap at OP5
OVERSPEED_RPM = 14400.0       # 1.2 * 12,000 rpm structural overspeed test
SF_STRESS_MIN = 1.50          # Min von Mises safety factor (sigma_y=420 MPa -> <= 280 MPa)
RIPPLE_MAX_PCT = 14.0         # Hard cap on WLTC-weighted torque ripple (strong reward < 5%)


def _sic_inverter_loss(I_rms: float, I_pk: float, S_inv_kVA: float) -> float:
    """Realistic 400 V SiC MOSFET traction inverter loss (conduction + 12 kHz switching + aux).
    R_ds,on scales inversely with SiC die area (S_inv_kVA): 3.0 mOhm @ 100 kVA."""
    R_ds = 0.30 / max(S_inv_kVA, 35.0)
    P_cond = 3.0 * (I_rms ** 2) * R_ds
    t_sw_eff = 35.0e-9  # 35 ns effective SiC rise+fall switching time
    P_sw = (6.0 / math.pi) * INVERTER_FSW_HZ * V_DC_BUS * I_pk * t_sw_eff
    P_aux = 25.0
    return P_cond + P_sw + P_aux


def _compute_active_bom(design: MotorDesign, mesh: MotorMeshBundle, l_ew: float) -> dict:
    """Compute exact motor active material masses (with 1.35 lamination stamping scrap) and BOM."""
    n_p = design.n_poles
    m_fe_net = n_p * L_STK * (mesh.area_sta_steel + mesh.area_rot_steel) * STEEL_RHO
    m_fe_gross = 1.35 * m_fe_net  # realistic 35% lamination stamping scrap
    m_cu = n_p * (L_STK + l_ew) * mesh.area_copper_1pole * CU_DENSITY
    mp = MAGNET_DB[design.magnet_grade]
    m_mag = n_p * L_STK * mesh.area_magnets * mp["rho"]

    bom_fe = m_fe_gross * STEEL_COST_PER_KG
    bom_cu = m_cu * CU_COST_PER_KG
    bom_mag = m_mag * mp["cost_per_kg"]
    bom_fixed = 95.0  # shaft, bearings, aluminum water-jacket housing, resolver
    bom_motor = bom_fe + bom_cu + bom_mag + bom_fixed
    return {
        "m_steel_net_kg": m_fe_net,
        "m_cu_kg": m_cu,
        "m_mag_kg": m_mag,
        "m_active_kg": m_fe_net + m_cu + m_mag,
        "bom_steel_usd": bom_fe,
        "bom_cu_usd": bom_cu,
        "bom_mag_usd": bom_mag,
        "bom_motor_usd": bom_motor,
    }


def _build_fast_flux_map(solver: MotorFEASolver, I_pk_max: float):
    """Evaluate a non-uniform clustered 4x4 nonlinear FEA flux grid at step_k=0
    (dense knots in the continuous knee region [0, 0.16, 0.44, 1.0]*I_max)
    and fit exact bicubic splines (kx=3, ky=3) for (lam_d(Id, Iq), lam_q(Id, Iq))."""
    fracs = np.array([0.0, 0.16, 0.44, 1.0], dtype=np.float64)
    id_grid = -I_pk_max * fracs[::-1]   # strictly ascending: [-I_max, ..., 0.0]
    iq_grid = I_pk_max * fracs          # strictly ascending: [0.0, ..., +I_max]
    n_g = len(fracs)
    ld_mat = np.zeros((n_g, n_g), dtype=np.float64)
    lq_mat = np.zeros((n_g, n_g), dtype=np.float64)
    A_prev = None
    ge = GAMMA_E0_RAD

    for i_d, Id in enumerate(id_grid):
        for i_q, Iq in enumerate(iq_grid):
            ia = Id * math.cos(ge) - Iq * math.sin(ge)
            ib = Id * math.cos(ge - 2.0 * math.pi / 3.0) - Iq * math.sin(ge - 2.0 * math.pi / 3.0)
            ic = -(ia + ib)
            res = solver.solve_step(0, (ia, ib, ic), A_init=A_prev, tol=5e-5, max_iter=14)
            A_prev = res["A_free"]
            ld_mat[i_d, i_q] = res["lam_d"] + solver.L_ew * Id
            lq_mat[i_d, i_q] = res["lam_q"] + solver.L_ew * Iq

    interp_ld = RectBivariateSpline(id_grid, iq_grid, ld_mat, kx=3, ky=3)
    interp_lq = RectBivariateSpline(id_grid, iq_grid, lq_mat, kx=3, ky=3)
    return id_grid, iq_grid, interp_ld, interp_lq


def _solve_mtpa_fw_from_map(
    solver: MotorFEASolver, interp_ld, interp_lq,
    rpm: float, T_target: float, I_pk_lim: float,
    V_lim: float = V_PH_PEAK_MAX * 0.975,
) -> Tuple[float, float, float]:
    """Find minimum-loss (MTPA / Field-Weakening) (Id, Iq) inside I <= I_pk_lim and V <= V_lim
    achieving shaft torque T_target (or maximum feasible torque if target unreachable)."""
    omega_m = rpm * 2.0 * math.pi / 60.0
    omega_e = solver.p_pairs * omega_m
    P_mech = 12.0 * (rpm / 1000.0) + 0.18 * ((rpm / 1000.0) ** 2.6)
    T_em_target = T_target + P_mech / max(omega_m, 1.0)

    id_v = np.linspace(-I_pk_lim, 0.0, 51)
    iq_v = np.linspace(0.0, I_pk_lim, 51)
    ID, IQ = np.meshgrid(id_v, iq_v, indexing="ij")
    I_mag = np.hypot(ID, IQ)
    LD = interp_ld(id_v, iq_v, grid=True)
    LQ = interp_lq(id_v, iq_v, grid=True)
    T_em = 1.5 * solver.p_pairs * (LD * IQ - LQ * ID)
    VD = solver.R_ph * ID - omega_e * LQ
    VQ = solver.R_ph * IQ + omega_e * LD
    V_mag = np.hypot(VD, VQ)

    feas = (I_mag <= I_pk_lim * 1.002) & (V_mag <= V_lim)
    hit = feas & (T_em >= T_em_target * 0.995)

    if np.any(hit):
        # Minimize combined copper + iron flux-squared proxy among feasible solutions hitting torque
        lam2_max = max(float(np.max(LD ** 2 + LQ ** 2)), 1e-6)
        score = I_mag ** 2 + 0.12 * ((omega_e / 1500.0) ** 2) * (LD ** 2 + LQ ** 2) * (I_pk_lim ** 2) / lam2_max
        score = np.where(hit, score, 1e30)
        idx = np.unravel_index(np.argmin(score), score.shape)
        id0, iq0 = float(ID[idx]), float(IQ[idx])
    elif np.any(feas):
        t_masked = np.where(feas, T_em, -1e9)
        idx = np.unravel_index(np.argmax(t_masked), t_masked.shape)
        return float(ID[idx]), float(IQ[idx]), float(T_em[idx])
    else:
        return -0.7 * I_pk_lim, 0.3 * I_pk_lim, 0.0

    # Local fine 15x15 refinement around (id0, iq0) for <0.2% torque & voltage accuracy
    d_step = I_pk_lim / 50.0
    id_f = np.linspace(max(-I_pk_lim, id0 - 1.6 * d_step), min(0.0, id0 + 1.6 * d_step), 15)
    iq_f = np.linspace(max(0.0, iq0 - 1.6 * d_step), min(I_pk_lim, iq0 + 1.6 * d_step), 15)
    IDf, IQf = np.meshgrid(id_f, iq_f, indexing="ij")
    Imf = np.hypot(IDf, IQf)
    LDf = interp_ld(id_f, iq_f, grid=True)
    LQf = interp_lq(id_f, iq_f, grid=True)
    Temf = 1.5 * solver.p_pairs * (LDf * IQf - LQf * IDf)
    Vf = np.hypot(solver.R_ph * IDf - omega_e * LQf, solver.R_ph * IQf + omega_e * LDf)
    feas_f = (Imf <= I_pk_lim * 1.002) & (Vf <= V_lim)
    hit_f = feas_f & (Temf >= T_em_target * 0.998)
    if np.any(hit_f):
        sc_f = np.where(hit_f, Imf ** 2 + 3.0 * (Temf - T_em_target) ** 2, 1e30)
        idx = np.unravel_index(np.argmin(sc_f), sc_f.shape)
        return float(IDf[idx]), float(IQf[idx]), float(Temf[idx])
    return id0, iq0, float(T_em[idx])


def evaluate_motor_design(
    design: MotorDesign,
    mesh_scale: float = 1.15,
    n_steps: int = 6,
    run_full_demag: bool = True,
) -> dict:
    """Full locked multiphysics evaluation of a MotorDesign across:
      1. Geometric manufacturability checks (bridges >= 0.80 mm, gap >= 0.65 mm)
      2. 2D plane-stress centrifugal elasticity FEA @ 14,400 rpm (SF >= 1.50)
      3. 4-point WLTC drive-cycle exact multi-step nonlinear 2D Maxwell FEA (HOT 140C)
         including cut-edge 1.75x Bertotti iron loss, AC proximity copper loss,
         3D end-winding resistance & inductance, harmonic current injection,
         SiC inverter loss & kVA sizing
      4. Peak launch overload @ 2,500 rpm (T >= 185 Nm within J <= 24 A/mm^2 & V <= 219.4 V)
      5. Worst-case magnet demagnetization test (-30C cold soak for Ferrite; +140C hot for NdFeB)
    """
    ok, msg, min_br, min_rib = validate_design_geometry(design)
    if not ok:
        return {"feasible": False, "reason": msg, "score": -1e6, "name": design.name}

    try:
        mesh = build_motor_mesh(design, mesh_scale=mesh_scale)
    except Exception as exc:
        return {"feasible": False, "reason": f"mesh_fail: {exc}", "score": -1e6, "name": design.name}

    # 1. Exact 2D Centrifugal Stress FEA @ 14,400 rpm (fast ~15 ms gate!)
    try:
        stress = solve_centrifugal_stress(design, mesh, overspeed_rpm=OVERSPEED_RPM)
    except Exception as exc:
        return {"feasible": False, "reason": f"stress_fail: {exc}", "score": -1e6, "name": design.name}

    vm_MPa = stress["vm_nodal_max_MPa"]
    sf_stress = stress["safety_factor"]

    # 2. Initialize Nonlinear Maxwell Solver @ HOT 140C (worst-case torque/flux temp!)
    solver = MotorFEASolver(design, mesh, mag_temp_mode="HOT_140C")
    slot_area = solver.slot_area_mean
    I_rms_cont_lim = J_RMS_CONT_MAX * (K_FILL_MAX * slot_area * 1e6) / design.N_turns
    I_rms_peak_lim = J_RMS_PEAK_MAX * (K_FILL_MAX * slot_area * 1e6) / design.N_turns
    I_pk_cont_lim = I_rms_cont_lim * math.sqrt(2.0)
    I_pk_peak_lim = I_rms_peak_lim * math.sqrt(2.0)

    # 3. Fast 4x4 nonlinear FEA flux map (~80 ms) to find exact MTPA / Field-Weakening currents
    try:
        _, _, interp_ld, interp_lq = _build_fast_flux_map(solver, I_pk_peak_lim)
    except Exception as exc:
        return {"feasible": False, "reason": f"fluxmap_fail: {exc}", "score": -1e6, "name": design.name}

    # 4. Solve Peak Launch OP5 @ 2,500 rpm (exact multi-step FEA)
    id_pk, iq_pk, _ = _solve_mtpa_fw_from_map(
        solver, interp_ld, interp_lq, PEAK_OP_RPM, 300.0, I_pk_peak_lim,
        V_lim=V_PH_PEAK_MAX * 0.99
    )
    # Also find minimum-current vector for exact PEAK_TORQUE_MIN_NM (so inverter kVA is sized fairly!)
    id_pk_req, iq_pk_req, t_map_pk = _solve_mtpa_fw_from_map(
        solver, interp_ld, interp_lq, PEAK_OP_RPM, PEAK_TORQUE_MIN_NM, I_pk_peak_lim,
        V_lim=V_PH_PEAK_MAX * 0.99
    )
    use_id_pk, use_iq_pk = (id_pk_req, iq_pk_req) if t_map_pk >= PEAK_TORQUE_MIN_NM * 0.97 else (id_pk, iq_pk)
    op5 = solver.solve_dq_sweep(
        use_id_pk, use_iq_pk, rpm=PEAK_OP_RPM, n_steps=n_steps,
        k_h5=design.k_h5, phi_h5=design.phi_h5_rad,
        k_h7=design.k_h7, phi_h7=design.phi_h7_rad,
    )
    # If slight shortfall (<=10%) due to 6-step saturation vs 4x4 map, scale current once up to I_pk_peak_lim
    if 0.85 * PEAK_TORQUE_MIN_NM < op5["T_avg"] < PEAK_TORQUE_MIN_NM * 0.995:
        scale = min(1.12, math.sqrt(PEAK_TORQUE_MIN_NM / max(op5["T_avg"], 1.0)))
        cur_mag = math.hypot(use_id_pk, use_iq_pk)
        if cur_mag * scale <= I_pk_peak_lim * 1.005:
            use_id_pk *= scale
            use_iq_pk *= scale
            op5 = solver.solve_dq_sweep(
                use_id_pk, use_iq_pk, rpm=PEAK_OP_RPM, n_steps=n_steps,
                k_h5=design.k_h5, phi_h5=design.phi_h5_rad,
                k_h7=design.k_h7, phi_h7=design.phi_h7_rad,
            )

    T_peak_shaft = op5["P_shaft_W"] / (PEAK_OP_RPM * 2.0 * math.pi / 60.0)

    # 5. Worst-Case Demagnetization Check at Peak Launch (COLD -30C for Ferrite, HOT +140C for NdFeB!)
    mp = MAGNET_DB[design.magnet_grade]
    B_knee = mp["Bk_demag_worst"]
    demag_margin_T = 99.0
    p02_B_demag = 99.0
    if run_full_demag and design.magnet_grade != "NONE":
        if design.magnet_grade == "FERRITE_FB12B":
            # Worst-case demagnetization for ceramic ferrite is COLD SOAK (-30C) at peak -Id!
            solver.set_magnet_temperature("COLD_M30C")
            # Evaluate worst demag angle under peak current + 15% fault surge (1.15x peak launch!)
            ge = GAMMA_E0_RAD
            id_f = 1.15 * use_id_pk
            iq_f = 1.15 * use_iq_pk
            ia = id_f * math.cos(ge) - iq_f * math.sin(ge)
            ib = id_f * math.cos(ge - 2.0 * math.pi / 3.0) - iq_f * math.sin(ge - 2.0 * math.pi / 3.0)
            ic = -(ia + ib)
            r_dem = solver.solve_step(0, (ia, ib, ic))
            p02_B_demag = r_dem["p02_B_mag"]
            solver.set_magnet_temperature("HOT_140C")
        else:
            # Worst-case demag for NdFeB is HOT +140C at 1.15x peak overload!
            ge = GAMMA_E0_RAD
            id_f = 1.15 * use_id_pk
            iq_f = 1.15 * use_iq_pk
            ia = id_f * math.cos(ge) - iq_f * math.sin(ge)
            ib = id_f * math.cos(ge - 2.0 * math.pi / 3.0) - iq_f * math.sin(ge - 2.0 * math.pi / 3.0)
            ic = -(ia + ib)
            r_dem = solver.solve_step(0, (ia, ib, ic))
            p02_B_demag = r_dem["p02_B_mag"]
        demag_margin_T = p02_B_demag - B_knee

    # 6. Evaluate All 4 WLTC Operating Points with Exact Multi-Step Nonlinear FEA
    op_results = {}
    max_I_pk = op5["I_pk"]
    max_V_pk_wltc = op5["V_ph_pk"]
    max_J_cont = 0.0
    wltc_torque_shortfall = 0.0

    for op_name, rpm, T_req, w in WLTC_OPS:
        id_op, iq_op, _ = _solve_mtpa_fw_from_map(
            solver, interp_ld, interp_lq, rpm, T_req, I_pk_cont_lim,
            V_lim=V_PH_PEAK_MAX * 0.98
        )
        # Realistic EV FOC controller: 5th/7th harmonic NVH cancellation active in low/mid-speed
        # audible torque band (rpm <= 5000) and cleanly ramped off in high-speed field weakening!
        use_k5 = design.k_h5 if rpm <= 5000.0 else 0.0
        use_k7 = design.k_h7 if rpm <= 5000.0 else 0.0
        sw = solver.solve_dq_sweep(
            id_op, iq_op, rpm=rpm, n_steps=n_steps,
            k_h5=use_k5, phi_h5=design.phi_h5_rad,
            k_h7=use_k7, phi_h7=design.phi_h7_rad,
        )
        # Single-step exact torque trim if FEA differs from target by >2% (keeps exact shaft power comparison!)
        T_sh = sw["P_shaft_W"] / (rpm * 2.0 * math.pi / 60.0)
        if 0.65 * T_req < T_sh < 1.35 * T_req and abs(T_sh - T_req) / T_req > 0.025:
            ratio = (T_req / max(T_sh, 1.0)) ** 0.65
            id_try, iq_try = id_op * ratio, iq_op * ratio
            if math.hypot(id_try, iq_try) <= I_pk_cont_lim * 1.02:
                sw2 = solver.solve_dq_sweep(
                    id_try, iq_try, rpm=rpm, n_steps=n_steps,
                    k_h5=use_k5, phi_h5=design.phi_h5_rad,
                    k_h7=use_k7, phi_h7=design.phi_h7_rad,
                )
                if sw2["V_ph_pk"] <= V_PH_PEAK_MAX * 1.01:
                    sw = sw2
                    id_op, iq_op = id_try, iq_try
                    T_sh = sw["P_shaft_W"] / (rpm * 2.0 * math.pi / 60.0)

        sw["Id"] = id_op
        sw["Iq"] = iq_op
        sw["T_shaft"] = T_sh
        sw["T_target"] = T_req
        sw["rpm"] = rpm
        sw["weight"] = w
        op_results[op_name] = sw

        max_I_pk = max(max_I_pk, sw["I_pk"])
        max_V_pk_wltc = max(max_V_pk_wltc, sw["V_ph_pk"])
        max_J_cont = max(max_J_cont, sw["J_rms_Amm2"])
        if T_sh < T_req * 0.95:
            wltc_torque_shortfall += (T_req - T_sh) / T_req

    # 7. Exact SiC Inverter Sizing (kVA rating at 400 V bus determined by required peak phase current!)
    S_inv_kVA = 1.5 * V_PH_PEAK_MAX * max_I_pk / 1000.0
    inv_cost_usd = INVERTER_FIXED_USD + INVERTER_USD_PER_KVA * S_inv_kVA

    # 8. Aggregate WLTC-Weighted Drive System Efficiency (Motor + SiC Inverter), cos(phi), Ripple, Losses
    P_shaft_w = 0.0
    P_elec_mot_w = 0.0
    P_elec_sys_w = 0.0
    P_cu_dc_w = 0.0
    P_cu_ac_w = 0.0
    P_fe_sta_w = 0.0
    P_fe_rot_w = 0.0
    P_inv_w = 0.0
    P_mech_w = 0.0
    cos_phi_w = 0.0
    ripple_w = 0.0

    for op_name, rpm, T_req, w in WLTC_OPS:
        sw = op_results[op_name]
        p_inv = _sic_inverter_loss(sw["I_rms"], sw["I_pk"], S_inv_kVA)
        sw["P_inv_loss_W"] = p_inv
        sw["eta_sys_pct"] = float(100.0 * sw["P_shaft_W"] / max(sw["P_elec_W"] + p_inv, 1.0)) if sw["P_shaft_W"] > 0 else 0.0

        P_shaft_w += w * max(sw["P_shaft_W"], 0.0)
        P_elec_mot_w += w * max(sw["P_elec_W"], 1.0)
        P_elec_sys_w += w * max(sw["P_elec_W"] + p_inv, 1.0)
        P_cu_dc_w += w * sw["P_cu_dc_W"]
        P_cu_ac_w += w * sw["P_cu_ac_W"]
        P_fe_sta_w += w * sw["P_fe_sta_W"]
        P_fe_rot_w += w * sw["P_fe_rot_W"]
        P_inv_w += w * p_inv
        P_mech_w += w * sw["P_mech_W"]
        cos_phi_w += w * sw["cos_phi"]
        ripple_w += w * sw["T_ripple_pct"]

    eta_motor_wltc = 100.0 * P_shaft_w / max(P_elec_mot_w, 1.0)
    eta_sys_wltc = 100.0 * P_shaft_w / max(P_elec_sys_w, 1.0)
    bom = _compute_active_bom(design, mesh, solver.l_ew)
    total_sys_cost_usd = bom["bom_motor_usd"] + inv_cost_usd

    # 9. Hard Anti-Reward-Hacking Feasibility Checks
    violations = []
    if sf_stress < SF_STRESS_MIN:
        violations.append(f"stress_SF={sf_stress:.2f}<{SF_STRESS_MIN}")
    if T_peak_shaft < PEAK_TORQUE_MIN_NM * 0.97:
        violations.append(f"T_peak={T_peak_shaft:.1f}<{PEAK_TORQUE_MIN_NM}")
    if op5["J_rms_Amm2"] > J_RMS_PEAK_MAX * 1.02:
        violations.append(f"J_peak={op5['J_rms_Amm2']:.1f}>{J_RMS_PEAK_MAX}")
    if max_J_cont > J_RMS_CONT_MAX * 1.02:
        violations.append(f"J_cont={max_J_cont:.1f}>{J_RMS_CONT_MAX}")
    if max_V_pk_wltc > V_PH_PEAK_MAX * 1.02:
        violations.append(f"V_pk={max_V_pk_wltc:.1f}>{V_PH_PEAK_MAX:.1f}")
    if wltc_torque_shortfall > 0.05:
        violations.append(f"wltc_shortfall={wltc_torque_shortfall:.2f}")
    if demag_margin_T < 0.0:
        violations.append(f"demag_margin={demag_margin_T:.3f}T<0")
    if ripple_w > RIPPLE_MAX_PCT:
        violations.append(f"ripple={ripple_w:.1f}%>{RIPPLE_MAX_PCT}%")

    feasible = (len(violations) == 0)

    # 10. Composite System Objective for Autoresearch:
    # Maximizes WLTC Drive-System Efficiency (Motor + SiC Inverter) while rewarding low System Cost
    # (Motor BOM + SiC Inverter kVA) and smooth NVH torque ripple (<5%), with smooth penalty ramp if infeasible!
    # Calibration: +1.0% WLTC sys efficiency ~ -$85 system cost ~ -3% ripple above 4%
    ripple_excess = max(0.0, ripple_w - 4.0)
    base_obj = (
        eta_sys_wltc
        - 0.0085 * (total_sys_cost_usd - 1000.0)
        - 0.14 * ripple_excess
        + 1.20 * (cos_phi_w - 0.80)
    )
    if not feasible:
        # Smooth continuous penalty so evolution climbs straight into the feasible basin!
        pen = 0.0
        if sf_stress < SF_STRESS_MIN:
            pen += 12.0 * (SF_STRESS_MIN - sf_stress)
        if T_peak_shaft < PEAK_TORQUE_MIN_NM:
            pen += 0.10 * (PEAK_TORQUE_MIN_NM - T_peak_shaft)
        if op5["J_rms_Amm2"] > J_RMS_PEAK_MAX:
            pen += 0.8 * (op5["J_rms_Amm2"] - J_RMS_PEAK_MAX)
        if max_J_cont > J_RMS_CONT_MAX:
            pen += 1.2 * (max_J_cont - J_RMS_CONT_MAX)
        if max_V_pk_wltc > V_PH_PEAK_MAX:
            pen += 0.15 * (max_V_pk_wltc - V_PH_PEAK_MAX)
        if wltc_torque_shortfall > 0:
            pen += 20.0 * wltc_torque_shortfall
        if demag_margin_T < 0.0:
            pen += 40.0 * abs(demag_margin_T)
        if ripple_w > RIPPLE_MAX_PCT:
            pen += 0.35 * (ripple_w - RIPPLE_MAX_PCT)
        score = base_obj - 25.0 - pen
    else:
        score = base_obj

    return {
        "name": design.name,
        "topology": f"{design.n_poles}p/{design.q_slots}s",
        "magnet_grade": design.magnet_grade,
        "n_barriers": len(design.barriers),
        "feasible": feasible,
        "violations": violations,
        "score": float(score),
        "eta_sys_wltc_pct": float(eta_sys_wltc),
        "eta_motor_wltc_pct": float(eta_motor_wltc),
        "cos_phi_wltc": float(cos_phi_w),
        "cos_phi_peak": float(op5["cos_phi"]),
        "ripple_wltc_pct": float(ripple_w),
        "ripple_op2_pct": float(op_results["OP2_SUBURBAN"]["T_ripple_pct"]),
        "T_peak_Nm": float(T_peak_shaft),
        "torque_density_Nm_L": float(T_peak_shaft / (math.pi * 0.100 ** 2 * L_STK * 1000.0)),
        "S_inv_kVA": float(S_inv_kVA),
        "I_inv_pk_A": float(max_I_pk),
        "V_pk_max_V": float(max_V_pk_wltc),
        "J_cont_max_Amm2": float(max_J_cont),
        "J_peak_Amm2": float(op5["J_rms_Amm2"]),
        "vm_stress_MPa": float(vm_MPa),
        "sf_stress": float(sf_stress),
        "demag_p02_B_T": float(p02_B_demag),
        "demag_margin_T": float(demag_margin_T),
        "min_bridge_mm": float(min_br),
        "min_rib_mm": float(min_rib),
        "R_ph_mOhm": float(solver.R_ph * 1e3),
        "L_ew_uH": float(solver.L_ew * 1e6),
        "bom": bom,
        "inv_cost_usd": float(inv_cost_usd),
        "total_sys_cost_usd": float(total_sys_cost_usd),
        "losses_wltc_W": {
            "Cu_DC": float(P_cu_dc_w),
            "Cu_AC": float(P_cu_ac_w),
            "Fe_Stator": float(P_fe_sta_w),
            "Fe_Rotor": float(P_fe_rot_w),
            "SiC_Inverter": float(P_inv_w),
            "Mech": float(P_mech_w),
            "Total": float(P_cu_dc_w + P_cu_ac_w + P_fe_sta_w + P_fe_rot_w + P_inv_w + P_mech_w),
        },
        "ops": {
            k: {
                "rpm": v["rpm"], "T_shaft": round(v["T_shaft"], 2), "T_target": v["T_target"],
                "Id": round(v["Id"], 1), "Iq": round(v["Iq"], 1),
                "V_ph_pk": round(v["V_ph_pk"], 1), "I_rms": round(v["I_rms"], 1),
                "cos_phi": round(v["cos_phi"], 3), "ripple_pct": round(v["T_ripple_pct"], 2),
                "eta_mot_pct": round(v["eta_motor_pct"], 2), "eta_sys_pct": round(v["eta_sys_pct"], 2),
                "P_cu_dc": round(v["P_cu_dc_W"], 1), "P_cu_ac": round(v["P_cu_ac_W"], 1),
                "P_fe_sta": round(v["P_fe_sta_W"], 1), "P_fe_rot": round(v["P_fe_rot_W"], 1),
                "P_inv": round(v["P_inv_loss_W"], 1),
            }
            for k, v in op_results.items()
        },
        "op5_peak": {
            "T_shaft": round(T_peak_shaft, 2), "Id": round(use_id_pk, 1), "Iq": round(use_iq_pk, 1),
            "V_ph_pk": round(op5["V_ph_pk"], 1), "I_pk": round(op5["I_pk"], 1),
            "J_rms": round(op5["J_rms_Amm2"], 2), "cos_phi": round(op5["cos_phi"], 3),
            "ripple_pct": round(op5["T_ripple_pct"], 2),
        },
    }
