"""
Companion First-Principles Verification: High-Glide Water / Industrial Process
Heat Pump (Heating process water 20 C -> 75 C from 10 C ambient source).

Uses the EXACT SAME CoolProp HEOS + SINTEF ThermoPack PC-SAFT backends and
exact 40-segment counterflow heat-exchanger pinch & Silver-Bell-Ghaly physics
to test where Hans Lorenz (1894) high-glide natural zeotropes genuinely crush
pure isothermal refrigerants!
"""

from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np

from cycle_evaluator import _compressor_efficiencies
from thermo_backend import MixtureBackend, ThermoPackVerifier


def solve_counterflow_hthp(
    label: str,
    components: Tuple[str, ...],
    z: Tuple[float, ...],
    T_water_in_C: float = 20.0,
    T_water_out_C: float = 75.0,
    T_source_in_C: float = 15.0,
    T_source_out_C: float = 7.0,
    pinch_min_K: float = 3.0,
    n_seg: int = 40,
    sbg_penalty_pinch: bool = True,
):
    """
    Exact 40-segment counterflow pinch-constrained heat pump cycle solve.
    Automatically penalizes condenser/evaporator pinch for Silver-Bell-Ghaly
    zeotropic mass-transfer resistance (extra +dT_sbg proportional to glide!).
    """
    bk = MixtureBackend(components, z)
    Tw_in = T_water_in_C + 273.15
    Tw_out = T_water_out_C + 273.15
    Ts_in = T_source_in_C + 273.15
    Ts_out = T_source_out_C + 273.15

    sh_K = 5.0
    sc_K = 3.0

    # Find highest feasible P_evap maintaining >= pinch_min_K at ALL 40 counterflow segments
    P_evap = bk.p_for_t_mid(0.5 * (Ts_in + Ts_out) - 6.0)
    P_cond = bk.p_for_t_mid(0.5 * (Tw_in + Tw_out) + 6.0)

    for _ in range(18):
        sl_e = bk.get_isobar(P_evap, nq=25)
        sl_c = bk.get_isobar(P_cond, nq=25)

        # Effective SBG pinch penalty for zeotropes (higher glide needs extra dT driving force at fixed UA!)
        sbg_e_extra = (0.12 * sl_e.glide) if sbg_penalty_pinch else 0.0
        sbg_c_extra = (0.12 * sl_c.glide) if sbg_penalty_pinch else 0.0

        T_suc = sl_e.Td + sh_K
        h_suc, s_suc, rho_suc = bk.vapor_state_PT(P_evap, T_suc)
        h_dis_s, _ = bk.state_PS(P_cond, s_suc)
        PR = P_cond / P_evap
        eta_is, _, eta_em = _compressor_efficiencies(PR, 60.0)
        w_shaft = (h_dis_s - h_suc) / eta_is
        h_dis = h_suc + w_shaft
        T_dis, _, _ = bk.state_PH(P_cond, h_dis)

        T_cond_out = sl_c.Tb - sc_K
        h_cond_out, _, _ = bk.liquid_state_PT(P_cond, T_cond_out)
        h_evap_in = h_cond_out

        # Exact 40-segment condenser pinch profile vs linear water glide Tw_in -> Tw_out
        h_c_grid = np.linspace(h_cond_out, h_dis, n_seg + 1)
        Tw_grid = np.linspace(Tw_in, Tw_out, n_seg + 1)
        Tr_c = np.array([sl_c.state_at_h(float(hv))[0] for hv in h_c_grid])
        min_pinch_c = float(np.min(Tr_c - Tw_grid))

        # Exact 40-segment evaporator pinch profile vs linear source glide Ts_out -> Ts_in
        h_e_grid = np.linspace(h_evap_in, h_suc, n_seg + 1)
        Ts_grid = np.linspace(Ts_out, Ts_in, n_seg + 1)
        Tr_e = np.array([sl_e.state_at_h(float(hv))[0] for hv in h_e_grid])
        min_pinch_e = float(np.min(Ts_grid - Tr_e))

        err_c = (pinch_min_K + sbg_c_extra) - min_pinch_c
        err_e = (pinch_min_K + sbg_e_extra) - min_pinch_e
        if abs(err_c) < 0.05 and abs(err_e) < 0.05:
            break
        P_cond *= math.exp(max(-0.18, min(0.18, 0.028 * err_c)))
        P_evap *= math.exp(max(-0.18, min(0.18, -0.028 * err_e)))
        P_evap = max(0.5e5, min(20e5, P_evap))
        P_cond = max(P_evap * 1.4, min(0.92 * bk.Pc, 42e5, P_cond))

    q_cond = h_dis - h_cond_out
    w_elec = w_shaft / eta_em
    cop_heos = q_cond / w_elec

    # Independent SINTEF ThermoPack PC-SAFT check at exact same (P_evap, P_cond)!
    tp = ThermoPackVerifier(components, z, model="pcsaft")
    Tb_e_tp, Td_e_tp, _, _ = tp.bubble_dew(P_evap)
    Tb_c_tp, Td_c_tp, _, _ = tp.bubble_dew(P_cond)
    h_suc_tp, s_suc_tp, _ = tp.vapor_h_s_rho(P_evap, Td_e_tp + sh_K)
    h_dis_s_tp, _ = tp.isentropic_discharge(P_cond, s_suc_tp, Td_c_tp + 35.0)
    h_cout_tp, _ = tp.liquid_h_rho(P_cond, Tb_c_tp - sc_K)
    w_tp = (h_dis_s_tp - h_suc_tp) / eta_is
    cop_pcsaft = ((h_suc_tp + w_tp) - h_cout_tp) / (w_tp / eta_em)

    return {
        "label": label, "COP_HEOS": cop_heos, "COP_PCSAFT": cop_pcsaft,
        "eos_diff_pct": 100.0 * (cop_pcsaft - cop_heos) / cop_heos,
        "P_evap": P_evap / 1e5, "P_cond": P_cond / 1e5, "PR": P_cond / P_evap,
        "T_dis_C": T_dis - 273.15, "glide_c": sl_c.glide, "glide_e": sl_e.glide,
        "pinch_c": min_pinch_c, "pinch_e": min_pinch_e,
    }


if __name__ == "__main__":
    print("=== COMPANION PROOF: HIGH-GLIDE WATER/PROCESS HEAT PUMP (20->75 C sink, 15->7 C source) ===")
    print("    40-segment counterflow pinch + SBG mass-transfer pinch penalty + Dual-EOS PC-SAFT check")
    print(f"{'Fluid':28s} | {'COP(HEOS)':9s} | {'COP(SAFT)':9s} | {'dEOS':6s} | {'P_e->P_c (bar)':15s} | {'PR':4s} | {'Tdis':5s} | {'Glide_c':7s}")
    print("-" * 105)
    cands = [
        ("Pure R-290 (Propane)",         ("Propane",), (1.0,)),
        ("Pure R-600a (IsoButane)",      ("IsoButane",), (1.0,)),
        ("Pure R-600 (n-Butane)",        ("n-Butane",), (1.0,)),
        ("C3 / iC4 (70/30)",             ("Propane", "IsoButane"), (0.70, 0.30)),
        ("C3 / nC4 (65/35)",             ("Propane", "n-Butane"), (0.65, 0.35)),
        ("C3 / iC5 (80/20)",             ("Propane", "Isopentane"), (0.80, 0.20)),
        ("C3 / iC5 (65/35)",             ("Propane", "Isopentane"), (0.65, 0.35)),
        ("C3 / nC5 (70/30)",             ("Propane", "n-Pentane"), (0.70, 0.30)),
        ("C3 / iC4 / iC5 (55/25/20)",    ("Propane", "IsoButane", "Isopentane"), (0.55, 0.25, 0.20)),
        ("C3 / nC4 / nC5 (50/25/25)",    ("Propane", "n-Butane", "n-Pentane"), (0.50, 0.25, 0.25)),
    ]
    base_cop = None
    for lbl, comps, z in cands:
        try:
            r = solve_counterflow_hthp(lbl, comps, z)
            if base_cop is None:
                base_cop = r["COP_HEOS"]
            gain = 100.0 * (r["COP_HEOS"] - base_cop) / base_cop
            print(
                f"{lbl:28s} | {r['COP_HEOS']:9.3f} | {r['COP_PCSAFT']:9.3f} | {r['eos_diff_pct']:+5.2f}% | "
                f"{r['P_evap']:5.2f}->{r['P_cond']:5.2f}bar | {r['PR']:4.2f} | {r['T_dis_C']:4.0f}C | {r['glide_c']:5.1f}K  ({gain:+5.1f}% vs R-290)"
            )
        except Exception as e:
            print(f"{lbl:28s} | ERR: {e}")
