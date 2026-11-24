"""
Stage 5 (Regime B): Independent 5-Test Verification Gauntlet for all 3
Industrial HTHP Duties (DAIRY_90C, DISTRICT_110C, CHEMICAL_135C) +
2x2 Publication Figure Generator + Engineering Spec Sheet Exporter.

Verifies every Stage-4 autoresearch winner via:
  Test 1: Independent Dual-EOS check (CoolProp HEOS vs SINTEF ThermoPack PC-SAFT)
  Test 2: BPHE Finite-Volume mesh convergence (n_seg = 20 -> 40 -> 80)
  Test 3: 20% Vapor-Phase Rayleigh Leak Fractionation + post-leak ATEX suction check
  Test 4: 64-draw Correlated Monte Carlo (HTC +/-25%, SBG 0.75..1.50x, dP 0.80..1.30x)
  Test 5: Gouy-Stodola Exergy Destruction Audit (kW) + Publication Figure
"""

from __future__ import annotations

import json
import math
import os
import random
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from brazed_plate_hx import BPHEGeometry
from hthp_evaluator import (
    INDUSTRIAL_DUTIES,
    HTHPDesign,
    HTHPEvaluation,
    _industrial_comp_eff,
    evaluate_hthp,
)
from thermo_backend import FLUID_METADATA, MixtureBackend, ThermoPackVerifier


def hthp_from_dict(d: Dict) -> HTHPDesign:
    def mk_hx(h: Dict) -> BPHEGeometry:
        return BPHEGeometry(
            name=h["name"], area_m2=float(h["area_m2"]),
            plate_length_m=float(h["plate_length_m"]), plate_width_m=float(h["plate_width_m"]),
            chevron_deg=float(h["chevron_deg"]), channel_gap_mm=float(h["channel_gap_mm"]),
            plate_thick_mm=float(h.get("plate_thick_mm", 0.40)),
        )

    return HTHPDesign(
        label=d["label"], duty_name=d["duty_name"],
        topology=d["topology"], ihx_eff=float(d["ihx_eff"]), ft_pr_bias=float(d["ft_pr_bias"]),
        components=tuple(d["components"]), mole_fractions=tuple(d["mole_fractions"]),
        cond_area_frac=float(d["cond_area_frac"]),
        cond_hx=mk_hx(d["cond_hx"]), evap_hx=mk_hx(d["evap_hx"]),
        superheat_K=float(d["superheat_K"]), subcool_K=float(d["subcool_K"]),
    )


def fluid_pretty(d: HTHPDesign) -> str:
    m = {
        "R1233zd(E)": "R-1233zd(E)", "R1336mzz(Z)": "R-1336mzz(Z)", "R245fa": "R-245fa",
        "CarbonDioxide": "CO2", "Propane": "R-290", "IsoButane": "R-600a", "n-Butane": "R-600",
        "Isopentane": "R-601a", "n-Pentane": "R-601", "n-Hexane": "R-602", "Ammonia": "R-717",
    }
    if len(d.components) == 1:
        return m.get(d.components[0], d.components[0])
    return " / ".join(f"{m[c]} {100*z:.1f}%" for c, z in zip(d.components, d.mole_fractions))


def dual_eos_pcsaft_check(d: HTHPDesign, ev: HTHPEvaluation) -> Dict:
    """Independently recompute thermodynamic COP & condenser/evaporator glides in SINTEF PC-SAFT."""
    if any(FLUID_METADATA[c]["tp_name"] is None for c in d.components):
        return {"supported": False, "reason": f"No ThermoPack map for {d.components}"}

    Pe = ev.P_evap_bar * 1e5
    Pc = ev.P_cond_bar * 1e5
    sh_K = max(5.0, d.superheat_K)
    sc_K = max(2.0, d.subcool_K)

    # Condenser bulk composition & evaporator composition (supports FTVI split!)
    z_cond = d.mole_fractions
    z_evap = ev.z_evap if ev.z_evap else d.mole_fractions

    tp_c = ThermoPackVerifier(d.components, z_cond, model="pcsaft")
    tp_e = ThermoPackVerifier(d.components, z_evap, model="pcsaft")
    bk_c = MixtureBackend(d.components, z_cond)
    bk_e = MixtureBackend(d.components, z_evap)

    Tb_e_tp, Td_e_tp, _, _ = tp_e.bubble_dew(Pe)
    Tb_c_tp, Td_c_tp, _, _ = tp_c.bubble_dew(Pc)
    sl_e_cp = bk_e.get_isobar(Pe, nq=17)
    sl_c_cp = bk_c.get_isobar(Pc, nq=17)

    # Independent 1-stage reference thermodynamic cycle COP at identical (Pe, Pc)
    # (isolates pure EOS enthalpy/entropy agreement without fitting parameters)
    h_suc_tp, s_suc_tp, _ = tp_c.vapor_h_s_rho(Pe, Td_e_tp + sh_K)
    h_dis_s_tp, _ = tp_c.isentropic_discharge(Pc, s_suc_tp, Td_c_tp + 35.0)
    h_cout_tp, _ = tp_c.liquid_h_rho(Pc, Tb_c_tp - sc_K)

    h_suc_cp, s_suc_cp, _ = bk_c.vapor_state_PT(Pe, sl_e_cp.Td + sh_K)
    h_dis_s_cp, _ = bk_c.state_PS(Pc, s_suc_cp)
    h_cout_cp, _, _ = bk_c.liquid_state_PT(Pc, sl_c_cp.Tb - sc_K)

    eta_is, _, eta_em = _industrial_comp_eff(Pc / Pe)
    w_tp = (h_dis_s_tp - h_suc_tp) / eta_is
    w_cp = (h_dis_s_cp - h_suc_cp) / eta_is
    cop_tp = ((h_suc_tp + w_tp) - h_cout_tp) / (w_tp / eta_em)
    cop_cp = ((h_suc_cp + w_cp) - h_cout_cp) / (w_cp / eta_em)

    return {
        "supported": True,
        "cop_heos_thermo": round(cop_cp, 3),
        "cop_pcsaft_thermo": round(cop_tp, 3),
        "cop_diff_pct": round(100.0 * (cop_tp - cop_cp) / cop_cp, 2),
        "glide_cond_heos_K": round(sl_c_cp.glide, 2),
        "glide_cond_pcsaft_K": round(Td_c_tp - Tb_c_tp, 2),
        "glide_evap_heos_K": round(sl_e_cp.glide, 2),
        "glide_evap_pcsaft_K": round(Td_e_tp - Tb_e_tp, 2),
    }


def _mc_worker(args):
    design, enforce_nat, hm, sm, dm = args
    ev = evaluate_hthp(
        design, enforce_natural_zero_pfas=enforce_nat, enforce_above_atm=True,
        n_seg=24, htc_mult=hm, sbg_mult=sm, dp_mult=dm,
    )
    return ev.COP_sys, ev.LCOH_usd_MWh, ev.feasible


def main():
    with open("hthp_stage4_winners.json") as f:
        data = json.load(f)

    print("=" * 114)
    print("STAGE 5 (REGIME B): INDEPENDENT 5-TEST INDUSTRIAL HTHP VERIFICATION GAUNTLET")
    print("=" * 114)

    report: Dict[str, Dict] = {}
    evals_80: Dict[str, Dict[str, HTHPEvaluation]] = {}
    designs: Dict[str, Dict[str, HTHPDesign]] = {}

    rng = random.Random(2026)
    draws = [(rng.uniform(0.75, 1.25), rng.uniform(0.75, 1.50), rng.uniform(0.80, 1.30)) for _ in range(64)]

    with ProcessPoolExecutor(max_workers=16) as pool:
        for duty_name in ("DAIRY_90C", "DISTRICT_110C", "CHEMICAL_135C"):
            print(f"\n{'#' * 114}\n# DUTY: {duty_name} (Sink {INDUSTRIAL_DUTIES[duty_name].T_sink_in_C:.0f}->{INDUSTRIAL_DUTIES[duty_name].T_sink_out_C:.0f}C | Source {INDUSTRIAL_DUTIES[duty_name].T_src_in_C:.0f}->{INDUSTRIAL_DUTIES[duty_name].T_src_out_C:.0f}C)\n{'#' * 114}")
            report[duty_name] = {}
            evals_80[duty_name] = {}
            designs[duty_name] = {}

            for track in ("PURE_CEILING", "NAT_ZEO_COP", "NAT_ZEO_LCOH"):
                d = hthp_from_dict(data[duty_name][track]["design"])
                designs[duty_name][track] = d
                enf_nat = (track != "PURE_CEILING")

                # Test 2: Mesh convergence n_seg = 20 -> 40 -> 80
                e20 = evaluate_hthp(d, enforce_natural_zero_pfas=enf_nat, n_seg=20)
                e40 = evaluate_hthp(d, enforce_natural_zero_pfas=enf_nat, n_seg=40)
                e80 = evaluate_hthp(d, enforce_natural_zero_pfas=enf_nat, n_seg=80)
                evals_80[duty_name][track] = e80
                drift = 100.0 * (e80.COP_sys - e20.COP_sys) / max(e20.COP_sys, 1e-6)

                # Test 1: Independent Dual-EOS PC-SAFT check
                tp_chk = dual_eos_pcsaft_check(d, e80)

                # Test 3: 20% Vapor-Phase Rayleigh Leak Fractionation + ATEX suction check
                leak_info = None
                if len(d.components) > 1:
                    bk = MixtureBackend(d.components, d.mole_fractions)
                    z_lk = bk.simulate_vapor_leak(0.20, T_leak=303.15)
                    d_lk = replace(d, label=f"{d.label}_LEAK20", mole_fractions=z_lk)
                    ev_lk = evaluate_hthp(d_lk, enforce_natural_zero_pfas=enf_nat, enforce_above_atm=False, n_seg=40)
                    leak_info = {
                        "z_orig": [round(v, 3) for v in d.mole_fractions],
                        "z_leaked": [round(v, 3) for v in z_lk],
                        "COP_orig": round(e80.COP_sys, 3),
                        "COP_leaked": round(ev_lk.COP_sys, 3),
                        "COP_delta_pct": round(100.0 * (ev_lk.COP_sys - e80.COP_sys) / max(e80.COP_sys, 1e-6), 2),
                        "P_suc_leaked_bar": round(ev_lk.P_evap_out_bar, 2),
                        "still_above_atm": bool(ev_lk.P_evap_out_bar >= 1.05),
                    }

                # Test 4: 64-draw Correlated Monte Carlo
                mc_raw = list(pool.map(_mc_worker, [(d, enf_nat, hm, sm, dm) for hm, sm, dm in draws]))
                cops = np.array([r[0] for r in mc_raw if r[0] > 0.5])
                lcohs = np.array([r[1] for r in mc_raw if r[0] > 0.5])

                report[duty_name][track] = {
                    "fluid": fluid_pretty(d),
                    "topology": d.topology,
                    "mesh_n20_n40_n80": [round(e20.COP_sys, 4), round(e40.COP_sys, 4), round(e80.COP_sys, 4)],
                    "mesh_drift_pct": round(drift, 3),
                    "COP_n80": round(e80.COP_sys, 4),
                    "eta_lorenz_pct": round(100.0 * e80.eta_lorenz, 2),
                    "LCOH_n80": round(e80.LCOH_usd_MWh, 2),
                    "CAPEX_usd": round(e80.CAPEX_usd, 0),
                    "VHC_kJ_m3": round(e80.VHC_kJ_m3, 0),
                    "P_suc_bar": round(e80.P_evap_out_bar, 2),
                    "P_cond_bar": round(e80.P_cond_bar, 2),
                    "T_dis_C": round(e80.T_dis_C, 1),
                    "glide_cond_K": round(e80.glide_cond_K, 1),
                    "glide_evap_K": round(e80.glide_evap_K, 1),
                    "dual_eos_pcsaft": tp_chk,
                    "rayleigh_leak_20pct": leak_info,
                    "mc_n_conv": int(len(cops)),
                    "mc_COP_mean": round(float(np.mean(cops)), 4),
                    "mc_COP_std": round(float(np.std(cops)), 4),
                    "mc_COP_p05_p95": [round(float(np.percentile(cops, 5)), 4), round(float(np.percentile(cops, 95)), 4)],
                    "mc_LCOH_mean": round(float(np.mean(lcohs)), 2),
                    "mc_all_cops": [r[0] for r in mc_raw],
                    "mc_all_lcohs": [r[1] for r in mc_raw],
                    "exergy_kW": {
                        "Comp": round(e80.I_comp_kW, 2), "Cond": round(e80.I_cond_kW, 2),
                        "Evap": round(e80.I_evap_kW, 2), "Exp": round(e80.I_exp_kW, 2),
                        "IHX_FT": round(e80.I_ihx_ft_kW, 2), "Pumps": round(e80.I_pumps_kW, 2),
                        "Total": round(e80.I_comp_kW + e80.I_cond_kW + e80.I_evap_kW + e80.I_exp_kW + e80.I_ihx_ft_kW + e80.I_pumps_kW, 2),
                    },
                }

                tp_str = (
                    f"PC-SAFT dCOP={tp_chk['cop_diff_pct']:+.2f}%, gl_c={tp_chk['glide_cond_heos_K']:.1f}/{tp_chk['glide_cond_pcsaft_K']:.1f}K"
                    if tp_chk["supported"] else "N/A (synthetic)"
                )
                lk_str = (
                    f"Leak20%: COP {leak_info['COP_orig']:.2f}->{leak_info['COP_leaked']:.2f} ({leak_info['COP_delta_pct']:+.1f}%), Psuc={leak_info['P_suc_leaked_bar']:.2f}b"
                    if leak_info else "Pure fluid"
                )
                print(
                    f"  {track:12s} | {d.topology:11s} {fluid_pretty(d):28s}\n"
                    f"    Mesh n20->n80: {e20.COP_sys:.3f}->{e80.COP_sys:.3f} ({drift:+.2f}%) | {100*e80.eta_lorenz:.1f}% Lorenz | LCOH=${e80.LCOH_usd_MWh:.2f}/MWh | VHC={e80.VHC_kJ_m3:.0f} | Psuc={e80.P_evap_out_bar:.2f}b | Tdis={e80.T_dis_C:.0f}C\n"
                    f"    Dual-EOS: {tp_str} | {lk_str}\n"
                    f"    MC(64): COP={np.mean(cops):.3f}+/-{np.std(cops):.3f} | Exergy={report[duty_name][track]['exergy_kW']['Total']:.1f}kW "
                    f"(Comp={e80.I_comp_kW:.1f} Cond={e80.I_cond_kW:.1f} Evap={e80.I_evap_kW:.1f} Exp={e80.I_exp_kW:.1f} IHX/FT={e80.I_ihx_ft_kW:.1f})"
                )

            # Pairwise MC win-rates vs PURE_CEILING
            p_cops = report[duty_name]["PURE_CEILING"]["mc_all_cops"]
            p_lcohs = report[duty_name]["PURE_CEILING"]["mc_all_lcohs"]
            c_cops = report[duty_name]["NAT_ZEO_COP"]["mc_all_cops"]
            l_lcohs = report[duty_name]["NAT_ZEO_LCOH"]["mc_all_lcohs"]
            cop_pairs = [(a, b) for a, b in zip(c_cops, p_cops) if a > 0.5 and b > 0.5]
            lcoh_pairs = [(a, b) for a, b in zip(l_lcohs, p_lcohs) if a < 500 and b < 500]
            w_cop = 100.0 * np.mean([a > b for a, b in cop_pairs])
            w_lcoh = 100.0 * np.mean([a < b for a, b in lcoh_pairs])
            d_cop_pct = 100.0 * (evals_80[duty_name]["NAT_ZEO_COP"].COP_sys - evals_80[duty_name]["PURE_CEILING"].COP_sys) / evals_80[duty_name]["PURE_CEILING"].COP_sys
            d_lcoh_pct = 100.0 * (evals_80[duty_name]["NAT_ZEO_LCOH"].LCOH_usd_MWh - evals_80[duty_name]["PURE_CEILING"].LCOH_usd_MWh) / evals_80[duty_name]["PURE_CEILING"].LCOH_usd_MWh
            report[duty_name]["pairwise_summary"] = {
                "COP_gain_vs_pure_opt_pct": round(d_cop_pct, 2),
                "LCOH_change_vs_pure_opt_pct": round(d_lcoh_pct, 2),
                "MC_COP_win_rate_pct": round(float(w_cop), 1),
                "MC_LCOH_win_rate_pct": round(float(w_lcoh), 1),
            }
            print(f"  ==> PAIRWISE VERDICT ({duty_name}): NAT_ZEO_COP {d_cop_pct:+.1f}% COP (MC win={w_cop:.1f}%) | NAT_ZEO_LCOH {d_lcoh_pct:+.1f}% LCOH (MC win={w_lcoh:.1f}%)")

    # Strip raw MC arrays before saving clean JSON
    clean_rep = {}
    for du, ddict in report.items():
        clean_rep[du] = {}
        for tk, v in ddict.items():
            if isinstance(v, dict) and "mc_all_cops" in v:
                clean_rep[du][tk] = {kk: vv for kk, vv in v.items() if not kk.startswith("mc_all_")}
            else:
                clean_rep[du][tk] = v
    with open("hthp_verification_report.json", "w") as f:
        json.dump(clean_rep, f, indent=2)

    # --- Generate 2x2 Industrial HTHP Publication Proof Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5), dpi=160)
    plt.rcParams.update({"font.size": 9.5})
    duties = ["DAIRY_90C", "DISTRICT_110C", "CHEMICAL_135C"]
    duty_lbls = ["Dairy / Food\n(35→90°C / 45→20°C)", "District / Drying\n(40→110°C / 55→25°C)", "Chemical / Paper\n(60→135°C / 75→40°C)"]
    x = np.arange(len(duties))
    w = 0.25

    # Panel A: System COP (n_seg=80) across all 3 industrial duties
    ax = axes[0, 0]
    y_p = [evals_80[d]["PURE_CEILING"].COP_sys for d in duties]
    y_c = [evals_80[d]["NAT_ZEO_COP"].COP_sys for d in duties]
    y_l = [evals_80[d]["NAT_ZEO_LCOH"].COP_sys for d in duties]
    ax.bar(x - w, y_p, w, label="Optimized Pure Ceiling (incl. PFAS)", color="#7f7f7f")
    ax.bar(x,     y_c, w, label="Autoresearch Natural Zeotrope (Max COP)", color="#1f77b4")
    ax.bar(x + w, y_l, w, label="Autoresearch Natural Zeotrope (Min LCOH)", color="#2ca02c")
    for i in range(len(duties)):
        ax.text(x[i] - w, y_p[i] + 0.08, f"{y_p[i]:.2f}", ha="center", fontsize=8)
        g_pct = 100.0 * (y_c[i] - y_p[i]) / y_p[i]
        ax.text(x[i], y_c[i] + 0.08, f"{y_c[i]:.2f}\n({g_pct:+.0f}%)", ha="center", fontsize=8, fontweight="bold", color="#1f77b4")
        ax.text(x[i] + w, y_l[i] + 0.08, f"{y_l[i]:.2f}", ha="center", fontsize=8, color="#2ca02c")
    ax.set_xticks(x); ax.set_xticklabels(duty_lbls, fontsize=8.5)
    ax.set_ylabel("System COP (incl. Water Pumps, n_seg=80)")
    ax.set_title("A. 100 kW Industrial HTHP System COP vs. Optimized Pure Ceiling")
    ax.set_ylim(0, max(y_c) * 1.24)
    ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8, loc="upper right")

    # Panel B: Levelized Cost of Heat (LCOH $/MWh_th) vs Industrial Natural Gas Boiler!
    ax = axes[0, 1]
    l_p = [evals_80[d]["PURE_CEILING"].LCOH_usd_MWh for d in duties]
    l_c = [evals_80[d]["NAT_ZEO_COP"].LCOH_usd_MWh for d in duties]
    l_l = [evals_80[d]["NAT_ZEO_LCOH"].LCOH_usd_MWh for d in duties]
    ax.bar(x - w, l_p, w, label="Optimized Pure Ceiling", color="#7f7f7f")
    ax.bar(x,     l_c, w, label="Natural Zeotrope (Max COP)", color="#1f77b4")
    ax.bar(x + w, l_l, w, label="Natural Zeotrope (Min LCOH)", color="#2ca02c")
    ax.axhline(50.0, color="#d62728", linestyle="--", linewidth=1.6, label=r"Natural Gas Boiler (~\$50/MWh @ \$50/t CO$_2$)")
    for i in range(len(duties)):
        ax.text(x[i] - w, l_p[i] + 0.8, f"${l_p[i]:.1f}", ha="center", fontsize=8)
        ax.text(x[i],     l_c[i] + 0.8, f"${l_c[i]:.1f}", ha="center", fontsize=8)
        ax.text(x[i] + w, l_l[i] + 0.8, f"${l_l[i]:.1f}", ha="center", fontsize=8, fontweight="bold", color="#2ca02c")
    ax.set_xticks(x); ax.set_xticklabels(duty_lbls, fontsize=8.5)
    ax.set_ylabel("Levelized Cost of Heat ($/MWh_th)")
    ax.set_title("B. True Industrial LCOH (Compressor VHC CAPEX + BPHE + $85/MWh Elec)")
    ax.set_ylim(0, max(max(l_p), 54.0) * 1.15)
    ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8, loc="upper left")

    # Panel C: Gouy-Stodola Exergy Destruction Breakdown (Pure Ceiling vs Natural Zeotrope)
    ax = axes[1, 0]
    ex_x = np.arange(6)
    ex_lbls = ["Dairy\nPure Opt", "Dairy\nZeo COP", "District\nPure Opt", "District\nZeo COP", "Chem 135\nPure Opt", "Chem 135\nZeo COP"]
    ex_pairs = [(d, t) for d in duties for t in ("PURE_CEILING", "NAT_ZEO_COP")]
    cats = ["Comp", "Cond", "Evap", "Exp", "IHX_FT", "Pumps"]
    cols = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#937860", "#8172b3"]
    bot = np.zeros(6)
    for cat, col in zip(cats, cols):
        vals = np.array([report[du][tk]["exergy_kW"][cat] for du, tk in ex_pairs])
        ax.bar(ex_x, vals, bottom=bot, label=cat, color=col, width=0.58)
        bot += vals
    for i, tv in enumerate(bot):
        ax.text(i, tv + 0.5, f"{tv:.1f}kW", ha="center", fontsize=8, fontweight="bold")
    ax.set_xticks(ex_x); ax.set_xticklabels(ex_lbls, fontsize=8.5)
    ax.set_ylabel("Exergy Destruction (kW per 100 kW_th Output)")
    ax.set_title("C. Gouy-Stodola Exergy Destruction: Pure Ceiling vs. Natural Zeotrope")
    ax.set_ylim(0, max(bot) * 1.18)
    ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8, ncol=3, loc="upper left")

    # Panel D: Independent Dual-EOS Verification (CoolProp HEOS vs SINTEF ThermoPack PC-SAFT)
    ax = axes[1, 1]
    tp_items = []
    for du, short_du in zip(duties, ["Dairy 90", "Dist 110", "Chem 135"]):
        for tk, short_tk in [("NAT_ZEO_COP", "MaxCOP"), ("NAT_ZEO_LCOH", "MinLCOH")]:
            tp = report[du][tk]["dual_eos_pcsaft"]
            if tp.get("supported"):
                tp_items.append((f"{short_du}\n{short_tk}", tp["cop_heos_thermo"], tp["cop_pcsaft_thermo"], tp["cop_diff_pct"]))
    ix = np.arange(len(tp_items))
    w2 = 0.36
    ax.bar(ix - w2 / 2, [it[1] for it in tp_items], w2, label="CoolProp HEOS (NIST)", color="#1f77b4")
    ax.bar(ix + w2 / 2, [it[2] for it in tp_items], w2, label="SINTEF ThermoPack (PC-SAFT)", color="#ff7f0e")
    for i, it in enumerate(tp_items):
        ax.text(i, max(it[1], it[2]) + 0.12, f"{it[3]:+.2f}%", ha="center", fontsize=8, fontweight="bold")
    ax.set_xticks(ix); ax.set_xticklabels([it[0] for it in tp_items], fontsize=8)
    ax.set_ylabel("Thermodynamic Cycle COP (at Identical P_e, P_c)")
    ax.set_title("D. Dual-EOS Cross-Check Across All 6 Industrial Zeotropic Winners")
    ax.set_ylim(0, max(max(it[1], it[2]) for it in tp_items) * 1.20)
    ax.grid(axis="y", alpha=0.3); ax.legend(fontsize=8, loc="upper right")

    plt.tight_layout()
    plt.savefig("hthp_autoresearch_proof.png", dpi=160)
    print("\nSaved hthp_verification_report.json and hthp_autoresearch_proof.png successfully!")


if __name__ == "__main__":
    main()
