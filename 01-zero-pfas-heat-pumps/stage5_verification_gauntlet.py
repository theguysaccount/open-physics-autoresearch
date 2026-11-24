"""
Stage 5: Independent 5-Test Verification Gauntlet & Publication Figure Generator.

Tests the exact Stage-4 autoresearch winners against all industry baselines:
  Test 1: Dual-EOS cross-check (CoolProp HEOS vs SINTEF ThermoPack PC-SAFT)
  Test 2: Finite-volume mesh convergence (n_seg = 10 -> 24 -> 60)
  Test 3: 20% Rayleigh vapor-leak fractionation stress test
  Test 4: +/-25% Two-phase HTC, +50% SBG, & +/-25% dP Monte Carlo robustness (64 draws)
  Test 5: Exact Gouy-Stodola Watt-by-Watt exergy destruction audit + figures
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

from baselines_and_profile import build_baselines
from cycle_evaluator import (
    STANDARD_BINS,
    SystemDesign,
    SystemEvaluation,
    _compressor_efficiencies,
    evaluate_system,
)
from microchannel_hx import CoilGeometry
from thermo_backend import FLUID_METADATA, MixtureBackend, ThermoPackVerifier


def design_from_dict(d: Dict) -> SystemDesign:
    def mk_coil(c: Dict) -> CoilGeometry:
        return CoilGeometry(
            name=c["name"], width=c["width"], height=c["height"],
            n_rows=int(c["n_rows"]), slab_depth=c["slab_depth"],
            tube_pitch=c["tube_pitch"], tube_height=c["tube_height"],
            port_dh=c["port_dh"], fin_pitch=c["fin_pitch"],
            pass_fractions=tuple(c["pass_fractions"]),
            counterflow_bridge=bool(c["counterflow_bridge"]),
            fin_thickness=c.get("fin_thickness", 0.00010),
            wall_thickness=c.get("wall_thickness", 0.00030),
            louver_pitch=c.get("louver_pitch", 0.0012),
            louver_angle_deg=c.get("louver_angle_deg", 27.0),
        )

    return SystemDesign(
        label=d["label"],
        topology=d["topology"],
        ihx_eff=float(d["ihx_eff"]),
        ft_pr_bias=float(d["ft_pr_bias"]),
        comp_shift_enabled=bool(d["comp_shift_enabled"]),
        components=tuple(d["components"]),
        mole_fractions=tuple(d["mole_fractions"]),
        indoor_coil=mk_coil(d["indoor_coil"]),
        outdoor_coil=mk_coil(d["outdoor_coil"]),
        comp_disp_cc=float(d["comp_disp_cc"]),
        v_air_in_m3s=tuple(d["v_air_in_m3s"]),
        v_air_out_m3s=tuple(d["v_air_out_m3s"]),
        superheat_K=tuple(d["superheat_K"]),
        subcool_K=tuple(d["subcool_K"]),
        liq_hold_frac=tuple(d["liq_hold_frac"]),
    )


def fluid_label(d: SystemDesign) -> str:
    short = {
        "R410A": "R-410A", "R32": "R-32", "CarbonDioxide": "CO2", "Ethane": "R-170",
        "Propylene": "R-1270", "Propane": "R-290", "IsoButane": "R-600a", "n-Butane": "R-600",
        "Isopentane": "R-601a", "n-Pentane": "R-601", "DimethylEther": "DME",
    }
    if len(d.components) == 1:
        return short[d.components[0]]
    return " / ".join(f"{short[c]} {100*z:.1f}%" for c, z in zip(d.components, d.mole_fractions))


def test1_dual_eos_crosscheck(d: SystemDesign, ev: SystemEvaluation) -> Dict:
    """Recompute thermodynamic cycle COP & phase glide independently in SINTEF ThermoPack PC-SAFT."""
    if any(FLUID_METADATA[c]["tp_name"] is None for c in d.components):
        return {"supported": False, "reason": f"ThermoPack lacks component mapping for {d.components}"}

    diffs_cop = []
    diffs_glide = []
    bin_rows = []
    for idx, (op_bin, b) in enumerate(zip(STANDARD_BINS, ev.bins)):
        if not b.feasible:
            continue
        z_use = b.z_circ
        tp = ThermoPackVerifier(d.components, z_use, model="pcsaft")
        bk = MixtureBackend(d.components, z_use)

        Pe = b.P_evap_bar * 1e5
        Pc = b.P_cond_bar * 1e5
        sh_K = max(3.5, d.superheat_K[idx])
        sc_K = max(2.0, d.subcool_K[idx])

        # Independent thermodynamic properties at identical (Pe, Pc, sh_K, sc_K)
        Tb_e_tp, Td_e_tp, hb_e_tp, hd_e_tp = tp.bubble_dew(Pe)
        Tb_c_tp, Td_c_tp, hb_c_tp, hd_c_tp = tp.bubble_dew(Pc)
        sl_e_cp = bk.get_isobar(Pe, nq=9)
        sl_c_cp = bk.get_isobar(Pc, nq=9)

        # Ideal 1-stage thermodynamic cycle COP comparison (pure EOS check, zero HX proxy)
        h_suc_tp, s_suc_tp, _ = tp.vapor_h_s_rho(Pe, Td_e_tp + sh_K)
        h_dis_s_tp, _ = tp.isentropic_discharge(Pc, s_suc_tp, Td_c_tp + 35.0)
        h_cond_out_tp, _ = tp.liquid_h_rho(Pc, Tb_c_tp - sc_K)

        h_suc_cp, s_suc_cp, _ = bk.vapor_state_PT(Pe, sl_e_cp.Td + sh_K)
        h_dis_s_cp, _ = bk.state_PS(Pc, s_suc_cp)
        h_cond_out_cp, _, _ = bk.liquid_state_PT(Pc, sl_c_cp.Tb - sc_K)

        eta_is, _, _ = _compressor_efficiencies(Pc / Pe, b.f_comp_Hz)
        w_tp = (h_dis_s_tp - h_suc_tp) / eta_is
        w_cp = (h_dis_s_cp - h_suc_cp) / eta_is
        q_tp = ((h_suc_tp + w_tp) - h_cond_out_tp) if op_bin.is_heating else (h_suc_tp - h_cond_out_tp)
        q_cp = ((h_suc_cp + w_cp) - h_cond_out_cp) if op_bin.is_heating else (h_suc_cp - h_cond_out_cp)
        cop_tp = q_tp / w_tp
        cop_cp = q_cp / w_cp

        rel_cop_pct = 100.0 * (cop_tp - cop_cp) / cop_cp
        glide_diff_K = (Td_e_tp - Tb_e_tp) - sl_e_cp.glide
        diffs_cop.append(rel_cop_pct)
        diffs_glide.append(glide_diff_K)
        bin_rows.append({
            "bin": op_bin.name,
            "cop_heos": round(cop_cp, 3),
            "cop_pcsaft": round(cop_tp, 3),
            "delta_pct": round(rel_cop_pct, 2),
            "glide_heos_K": round(sl_e_cp.glide, 2),
            "glide_pcsaft_K": round(Td_e_tp - Tb_e_tp, 2),
            "glide_diff_K": round(glide_diff_K, 2),
        })

    return {
        "supported": True,
        "mean_abs_cop_diff_pct": float(np.mean(np.abs(diffs_cop))),
        "max_abs_cop_diff_pct": float(np.max(np.abs(diffs_cop))),
        "mean_abs_glide_diff_K": float(np.mean(np.abs(diffs_glide))),
        "bins": bin_rows,
    }


def _mc_worker(args):
    design, enforce_gwp, htc_m, sbg_m, dp_m = args
    ev = evaluate_system(
        design, enforce_natural_gwp=enforce_gwp, max_glide_K=12.0,
        n_seg=12, htc_mult=htc_m, sbg_mult=sbg_m, dp_mult=dp_m,
    )
    return ev.ASPF, ev.COP_46C_T3, ev.COP_m7C_cold, ev.feasible


def main():
    with open("stage4_winners.json") as f:
        data = json.load(f)

    baselines = build_baselines()
    r32_opt = design_from_dict(data["R32_OPT"]["design"])
    pure_opt = design_from_dict(data["PURE_NAT_OPT"]["design"])
    zeo_opt = design_from_dict(data["ZEOTROPE_FULL_OPT"]["design"])

    # Also grab best ThermoPack-compatible zeotropic candidate (no Propylene) for Dual-EOS PC-SAFT proof
    zeo_tp_cand = None
    for item in data.get("ZEOTROPE_DIVERSE_TOP", []):
        cand = design_from_dict(item["design"])
        if all(FLUID_METADATA[c]["tp_name"] is not None for c in cand.components):
            zeo_tp_cand = cand
            break

    all_models: List[Tuple[str, SystemDesign, bool]] = [
        ("B1_R410A_TXV", baselines[0], False),
        ("B2_R32_TXV", baselines[1], False),
        ("B3_R290_TXV", baselines[2], True),
        ("B4_R290_FTVI_IHX", baselines[3], True),
        ("T1_R32_OPT", r32_opt, False),
        ("T2_PURE_NAT_OPT", pure_opt, True),
        ("T3_ZEOTROPE_OPT", zeo_opt, True),
    ]
    if zeo_tp_cand is not None and zeo_tp_cand.components != zeo_opt.components:
        all_models.append(("T3b_ZEO_PCSAFT_CHK", zeo_tp_cand, True))

    print("=" * 110)
    print("STAGE 5: INDEPENDENT 5-TEST VERIFICATION GAUNTLET")
    print("=" * 110)

    # --- TEST 2: Mesh Convergence (n_seg = 10 -> 24 -> 60) ---
    print("\n[TEST 2] Finite-Volume Mesh Convergence Audit (n_seg = 10 -> 24 -> 60):")
    evals_60: Dict[str, SystemEvaluation] = {}
    mesh_report: Dict[str, Dict] = {}
    for key, d, gwp_flag in all_models:
        e10 = evaluate_system(d, enforce_natural_gwp=gwp_flag, max_glide_K=10.0, n_seg=10)
        e24 = evaluate_system(d, enforce_natural_gwp=gwp_flag, max_glide_K=10.0, n_seg=24)
        e60 = evaluate_system(d, enforce_natural_gwp=gwp_flag, max_glide_K=10.0, n_seg=60)
        evals_60[key] = e60
        drift = 100.0 * (e60.ASPF - e10.ASPF) / max(e10.ASPF, 1e-6)
        mesh_report[key] = {
            "n10": e10.ASPF, "n24": e24.ASPF, "n60": e60.ASPF, "drift_pct": drift,
            "SEER": e60.SEER_cool, "SCOP": e60.SCOP_heat,
            "COP_46C": e60.COP_46C_T3, "COP_m7C": e60.COP_m7C_cold,
            "charge_g": e60.flammable_charge_g, "glide_K": e60.max_glide_K,
            "T_dis_C": e60.max_T_dis_C, "P_high_bar": e60.max_P_high_bar,
            "Al_kg": e60.al_mass_in_kg + e60.al_mass_out_kg, "inc_bom_usd": e60.inc_bom_usd,
            "feasible": e60.feasible, "fluid": fluid_label(d), "topology": d.topology,
        }
        print(
            f"  {key:18s} ({fluid_label(d):24s} | {d.topology:11s}): "
            f"n10={e10.ASPF:.4f} -> n24={e24.ASPF:.4f} -> n60={e60.ASPF:.4f} (drift={drift:+.2f}%) | "
            f"SEER={e60.SEER_cool:.3f} SCOP={e60.SCOP_heat:.3f} 46C={e60.COP_46C_T3:.2f} -7C={e60.COP_m7C_cold:.2f} "
            f"chrg={e60.flammable_charge_g:.0f}g gl={e60.max_glide_K:.1f}K +BOM=${e60.inc_bom_usd:.0f}"
        )

    # --- TEST 1: Dual-EOS Cross-Check (CoolProp HEOS vs SINTEF ThermoPack PC-SAFT) ---
    print("\n[TEST 1] Independent Dual-EOS Verification (CoolProp HEOS vs SINTEF ThermoPack PC-SAFT):")
    dual_eos_results = {}
    for key, d, _ in all_models:
        if key in ("B3_R290_TXV", "T2_PURE_NAT_OPT", "T3_ZEOTROPE_OPT", "T3b_ZEO_PCSAFT_CHK"):
            res_tp = test1_dual_eos_crosscheck(d, evals_60[key])
            dual_eos_results[key] = res_tp
            if res_tp["supported"]:
                print(
                    f"  {key:18s} ({fluid_label(d):24s}): "
                    f"Mean |dCOP| = {res_tp['mean_abs_cop_diff_pct']:.2f}% (max {res_tp['max_abs_cop_diff_pct']:.2f}%), "
                    f"Mean |dGlide| = {res_tp['mean_abs_glide_diff_K']:.2f} K"
                )
            else:
                print(f"  {key:18s}: skipped ({res_tp['reason']})")

    # --- TEST 3: 20% Vapor-Phase Rayleigh Leak Fractionation Test ---
    print("\n[TEST 3] 20% Vapor-Phase Rayleigh Leak Fractionation Stress Test:")
    leak_results = {}
    for key, d, gwp_flag in [("T3_ZEOTROPE_OPT", zeo_opt, True)] + (
        [("T3b_ZEO_PCSAFT_CHK", zeo_tp_cand, True)] if zeo_tp_cand is not None and zeo_tp_cand.components != zeo_opt.components else []
    ):
        bk = MixtureBackend(d.components, d.mole_fractions)
        z_leaked = bk.simulate_vapor_leak(leak_fraction=0.20, T_leak=298.15)
        d_leaked = replace(d, label=f"{d.label}_LEAK20", mole_fractions=z_leaked)
        ev_orig = evals_60[key]
        ev_leak = evaluate_system(d_leaked, enforce_natural_gwp=gwp_flag, max_glide_K=15.0, n_seg=24)
        drop_pct = 100.0 * (ev_leak.ASPF - ev_orig.ASPF) / max(ev_orig.ASPF, 1e-6)
        leak_results[key] = {
            "z_orig": d.mole_fractions, "z_leaked": z_leaked,
            "ASPF_orig": ev_orig.ASPF, "ASPF_leaked": ev_leak.ASPF, "delta_pct": drop_pct,
        }
        print(
            f"  {key:18s}: z_orig={[round(v,3) for v in d.mole_fractions]} -> z_leaked={[round(v,3) for v in z_leaked]} | "
            f"ASPF: {ev_orig.ASPF:.4f} -> {ev_leak.ASPF:.4f} ({drop_pct:+.2f}%)"
        )

    # --- TEST 4: +/-25% Two-Phase HTC, +50% SBG, & +/-25% dP Monte Carlo Audit ---
    print("\n[TEST 4] Correlated Monte Carlo Robustness Audit (64 draws: HTC +/-25%, SBG 0.75..1.50x, dP 0.80..1.30x):")
    rng = random.Random(2026)
    draws = [(rng.uniform(0.75, 1.25), rng.uniform(0.75, 1.50), rng.uniform(0.80, 1.30)) for _ in range(64)]
    mc_summary = {}
    with ProcessPoolExecutor(max_workers=16) as pool:
        for key, d, gwp_flag in [("B2_R32_TXV", baselines[1], False), ("B3_R290_TXV", baselines[2], True),
                                 ("B4_R290_FTVI_IHX", baselines[3], True),
                                 ("T1_R32_OPT", r32_opt, False), ("T2_PURE_NAT_OPT", pure_opt, True),
                                 ("T3_ZEOTROPE_OPT", zeo_opt, True)]:
            res = list(pool.map(_mc_worker, [(d, gwp_flag, hm, sm, dm) for hm, sm, dm in draws]))
            aspfs = np.array([r[0] for r in res if r[3]])
            n_ok = len(aspfs)
            if n_ok == 0:
                aspfs = np.array([0.0])
            mc_summary[key] = {
                "n_conv": n_ok,
                "mean": float(np.mean(aspfs)), "p05": float(np.percentile(aspfs, 5)),
                "p95": float(np.percentile(aspfs, 95)), "std": float(np.std(aspfs)),
                "all_vals": [r[0] for r in res],
            }
            print(
                f"  {key:18s} ({n_ok:2d}/64 conv): ASPF = {np.mean(aspfs):.4f} +/- {np.std(aspfs):.4f} "
                f"(90% CI: [{np.percentile(aspfs,5):.4f}, {np.percentile(aspfs,95):.4f}])"
            )

    # Correlated pairwise win-rates (only over mutually converged draws for strict fairness!)
    def _pair_win(k_a: str, k_b: str) -> Tuple[float, float, int]:
        va, vb = mc_summary[k_a]["all_vals"], mc_summary[k_b]["all_vals"]
        pairs = [(a, b) for a, b in zip(va, vb) if a > 0.5 and b > 0.5]
        if not pairs:
            return 0.0, 0.0, 0
        wins = 100.0 * np.mean([a > b for a, b in pairs])
        dm = float(np.mean([a - b for a, b in pairs]))
        return wins, dm, len(pairs)

    for opp in ("B2_R32_TXV", "B3_R290_TXV", "B4_R290_FTVI_IHX", "T1_R32_OPT", "T2_PURE_NAT_OPT"):
        w, dm, n_p = _pair_win("T3_ZEOTROPE_OPT", opp)
        print(f"  Correlated Win-Rate T3_ZEOTROPE_OPT vs {opp:18s}: {w:5.1f}% (mean delta={dm:+.4f} across {n_p} paired draws)")

    # --- TEST 5: Gouy-Stodola Exergy Destruction Breakdown (Load-Weighted Annual Average W) ---
    print("\n[TEST 5] Load-Weighted Annual Gouy-Stodola Exergy Destruction Audit (W):")
    exergy_summary = {}
    for key in ("B1_R410A_TXV", "B2_R32_TXV", "B3_R290_TXV", "B4_R290_FTVI_IHX", "T1_R32_OPT", "T2_PURE_NAT_OPT", "T3_ZEOTROPE_OPT"):
        ev = evals_60[key]
        comps = {k: 0.0 for k in ("Comp", "Cond", "Evap", "Exp", "IHX_FT", "Fans", "Backup")}
        for op_bin, b in zip(STANDARD_BINS, ev.bins):
            w = op_bin.weight
            comps["Comp"] += w * b.I_comp_W
            comps["Cond"] += w * b.I_cond_W
            comps["Evap"] += w * b.I_evap_W
            comps["Exp"] += w * b.I_exp_W
            comps["IHX_FT"] += w * b.I_ihx_ft_W
            comps["Fans"] += w * b.I_fans_W
            comps["Backup"] += w * b.W_backup_W
        tot = sum(comps.values())
        comps["Total"] = tot
        exergy_summary[key] = comps
        print(
            f"  {key:18s}: Total={tot:6.1f}W | Comp={comps['Comp']:5.1f} Cond={comps['Cond']:5.1f} "
            f"Evap={comps['Evap']:5.1f} Exp={comps['Exp']:5.1f} IHX/FT={comps['IHX_FT']:4.1f} Fans={comps['Fans']:5.1f}"
        )

    # Save full verification JSON
    out_json = {
        "mesh_convergence_n60": mesh_report,
        "dual_eos_pcsaft": dual_eos_results,
        "rayleigh_leak_20pct": leak_results,
        "monte_carlo_robustness": {k: {mk: mv for mk, mv in v.items() if mk != "all_vals"} for k, v in mc_summary.items()},
        "exergy_destruction_W": exergy_summary,
    }
    with open("verification_report.json", "w") as f:
        json.dump(out_json, f, indent=2)

    # Generate 2x2 publication-grade summary plot (both regimes + exergy + autoresearch curve)
    from test_high_glide_hthp import solve_counterflow_hthp

    hthp_cands = [
        ("Pure R-600a", ("IsoButane",), (1.0,)),
        ("Pure R-600", ("n-Butane",), (1.0,)),
        ("Pure R-290", ("Propane",), (1.0,)),
        ("C3/nC4\n55/45", ("Propane", "n-Butane"), (0.55, 0.45)),
        ("C3/iC4/iC5\n45/25/30", ("Propane", "IsoButane", "Isopentane"), (0.45, 0.25, 0.30)),
        ("C3/nC4/nC5\n45/25/30", ("Propane", "n-Butane", "n-Pentane"), (0.45, 0.25, 0.30)),
        ("C3/iC5\n60/40", ("Propane", "Isopentane"), (0.60, 0.40)),
    ]
    hthp_res = [
        solve_counterflow_hthp(lbl, c, z, T_water_in_C=35.0, T_water_out_C=85.0, T_source_in_C=45.0, T_source_out_C=15.0)
        for lbl, c, z in hthp_cands
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.5), dpi=160)
    plt.rcParams.update({"font.size": 9.5})

    # Panel A: ASPF + 46C Heatwave + -7C Cold Heating Comparison
    ax = axes[0, 0]
    keys_plot = ["B1_R410A_TXV", "B2_R32_TXV", "B3_R290_TXV", "B4_R290_FTVI_IHX", "T1_R32_OPT", "T2_PURE_NAT_OPT", "T3_ZEOTROPE_OPT"]
    labels_short = ["R-410A\n(Std TXV)", "R-32\n(Std TXV)", "R-290\n(Std TXV)", "R-290\n(FTVI+IHX)", "R-32\n(AutoOpt)", "Pure Nat.\n(AutoOpt)", "4-Layer\nNatural"]
    x = np.arange(len(keys_plot))
    w_bar = 0.26
    ax.bar(x - w_bar, [evals_60[k].ASPF for k in keys_plot], w_bar, label="Blended ASPF", color="#1f77b4")
    ax.bar(x, [evals_60[k].COP_46C_T3 for k in keys_plot], w_bar, label="46°C Heatwave COP", color="#d62728")
    ax.bar(x + w_bar, [evals_60[k].COP_m7C_cold for k in keys_plot], w_bar, label="-7°C Cold Heat COP", color="#2ca02c")
    ax.set_xticks(x)
    ax.set_xticklabels(labels_short, fontsize=8.5)
    ax.set_ylabel("System COP (incl. Fans & Backup)")
    ax.set_title("A. Air-to-Air Reversible Split AC/HP (6-Bin, n_seg=60)")
    ax.set_ylim(0, 6.0)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    for idx_k, k in enumerate(keys_plot):
        ax.text(idx_k - w_bar, evals_60[k].ASPF + 0.08, f"{evals_60[k].ASPF:.2f}", ha="center", fontsize=8, fontweight="bold")

    # Panel B: Stacked Annual Exergy Destruction (Gouy-Stodola Theorem)
    ax = axes[0, 1]
    ex_keys = ["B1_R410A_TXV", "B2_R32_TXV", "B3_R290_TXV", "B4_R290_FTVI_IHX", "T2_PURE_NAT_OPT", "T3_ZEOTROPE_OPT"]
    ex_labels = ["R-410A\nStd", "R-32\nStd", "R-290\nStd", "R-290\nFTVI", "Pure Nat\nAutoOpt", "4-Layer\nNatural"]
    cats = ["Comp", "Cond", "Evap", "Exp", "Fans", "IHX_FT"]
    colors = ["#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b3", "#937860"]
    bottom = np.zeros(len(ex_keys))
    for cat, col in zip(cats, colors):
        vals = np.array([exergy_summary[k][cat] for k in ex_keys])
        ax.bar(np.arange(len(ex_keys)), vals, bottom=bottom, label=cat, color=col, width=0.56)
        bottom += vals
    for i_x, tot_v in enumerate(bottom):
        ax.text(i_x, tot_v + 15, f"{tot_v:.0f}W", ha="center", fontsize=8, fontweight="bold")
    ax.set_xticks(np.arange(len(ex_keys)))
    ax.set_xticklabels(ex_labels, fontsize=8.5)
    ax.set_ylabel("Annual Weighted Exergy Destruction (W)")
    ax.set_title("B. Gouy-Stodola Exergy Destruction Audit (Exact Post-IHX)")
    ax.set_ylim(0, max(bottom) * 1.15)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, ncol=3, loc="upper right")

    # Panel C: Autoresearch Evolution Trajectory from autoresearch_log.tsv
    ax = axes[1, 0]
    track_pts: Dict[str, Tuple[List[int], List[float]]] = {}
    if os.path.exists("autoresearch_log.tsv"):
        with open("autoresearch_log.tsv") as f:
            _ = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                tr, ev_n, feas, aspf = parts[0], int(parts[2]), int(parts[3]), float(parts[4])
                if feas == 1:
                    track_pts.setdefault(tr, ([], []))
                    track_pts[tr][0].append(ev_n)
                    track_pts[tr][1].append(aspf)
    style_map = {
        "R32_OPT": ("#ff7f0e", "Track 1: R-32 AutoOpt (GWP 675, 1120g)"),
        "PURE_NAT_OPT": ("#2ca02c", "Track 2: Pure Natural AutoOpt (GWP ≤3, ≤500g)"),
        "ZEOTROPE_FULL_OPT": ("#1f77b4", "Track 3: 4-Layer Co-Design (GWP 2.8, 479g)"),
    }
    for tr, (xs, ys) in track_pts.items():
        col, lbl = style_map.get(tr, ("gray", tr))
        ax.step(xs, ys, where="post", label=lbl, color=col, linewidth=2.2)
        ax.scatter(xs, ys, color=col, s=22)
    ax.axhline(evals_60["B2_R32_TXV"].ASPF, color="gray", linestyle="--", linewidth=1.2, label=f"Std R-32 Baseline ({evals_60['B2_R32_TXV'].ASPF:.2f})")
    ax.set_xlabel("Cumulative Exact CoolProp HEOS Evaluations")
    ax.set_ylabel("Feasible Blended ASPF")
    ax.set_title("C. Stage-4 Autoresearch Optimization Trajectory (2,382 Evals)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")

    # Panel D: Dual-EOS High-Glide Industrial Waste-Heat Recovery Lorenz Proof
    ax = axes[1, 1]
    hx_x = np.arange(len(hthp_res))
    w2 = 0.36
    ax.bar(hx_x - w2 / 2, [r["COP_HEOS"] for r in hthp_res], w2, label="CoolProp HEOS (NIST)", color="#1f77b4")
    ax.bar(hx_x + w2 / 2, [r["COP_PCSAFT"] for r in hthp_res], w2, label="SINTEF ThermoPack (PC-SAFT)", color="#ff7f0e")
    ax.set_xticks(hx_x)
    ax.set_xticklabels([r["label"] for r in hthp_res], fontsize=8)
    ax.set_ylabel("Heating COP (incl. SBG Pinch Penalty)")
    ax.set_title("D. Dual-Glide Industrial Heat Pump (35→85°C / 45→15°C): Lorenz +35%")
    ax.set_ylim(0, 4.7)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    for idx_h, r in enumerate(hthp_res):
        ax.text(idx_h, max(r["COP_HEOS"], r["COP_PCSAFT"]) + 0.08, f"{r['COP_HEOS']:.2f}\n({r['glide_c']:.0f}K)", ha="center", fontsize=7.5, fontweight="bold")

    plt.tight_layout()
    plt.savefig("lorenz_autoresearch_proof.png", dpi=160)
    print("\nSaved verification_report.json and 2x2 lorenz_autoresearch_proof.png successfully!")


if __name__ == "__main__":
    main()
