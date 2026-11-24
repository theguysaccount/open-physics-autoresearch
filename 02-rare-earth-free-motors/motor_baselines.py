"""
motor_baselines.py — Stage-3 Industry Baselines & Profiling Micro-Experiments.
Benchmarks B1 (NdFeB V-IPM), B2 (Pure SynRM 4p), B3 (Ferrite PMa-SynRM 6p),
and B4 (Ferrite PMa-SynRM 4p) on the locked 5-point WLTC + SiC Inverter + Stress + Demag harness,
and quantifies loss breakdowns, wall-clock evaluation speed, and hypothesis micro-experiments.
"""

from __future__ import annotations
import copy
import json
import time
from motor_geometry import (
    preset_ndfeb_vipm, preset_pure_synrm_4p,
    preset_ferrite_pmasynrm_6p, preset_ferrite_pmasynrm_4p,
)
from motor_evaluator import evaluate_motor_design


def print_row(r: dict, dt_s: float):
    feas = "OK" if r["feasible"] else "FAIL:" + ",".join(r["violations"][:3])
    L = r["losses_wltc_W"]
    print(
        f"{r['name']:27s} | {r['topology']:6s} | eta_sys={r['eta_sys_wltc_pct']:5.2f}% "
        f"(mot={r['eta_motor_wltc_pct']:5.2f}%) | cos={r['cos_phi_wltc']:.3f} | Rip={r['ripple_wltc_pct']:4.1f}% "
        f"| Tpk={r['T_peak_Nm']:5.1f}Nm | Inv={r['S_inv_kVA']:5.1f}kVA | SysCost=${r['total_sys_cost_usd']:6.0f} "
        f"| SF={r['sf_stress']:.2f} | Demag={r['demag_margin_T']:+.2f}T | {dt_s:.2f}s | {feas}"
    )
    print(
        f"    Losses(W): CuDC={L['Cu_DC']:5.0f} CuAC={L['Cu_AC']:4.0f} FeSta={L['Fe_Stator']:5.0f} "
        f"FeRot={L['Fe_Rotor']:4.0f} SiC={L['SiC_Inverter']:4.0f} Mech={L['Mech']:3.0f} | "
        f"BOM: Fe=${r['bom']['bom_steel_usd']:.0f} Cu=${r['bom']['bom_cu_usd']:.0f} "
        f"Mag=${r['bom']['bom_mag_usd']:.0f} ({r['bom']['m_mag_kg']:.2f}kg) Inv=${r['inv_cost_usd']:.0f}"
    )


def main():
    print("=== STAGE-3 BASELINES & PROFILING (Exact 5-Point Nonlinear FEA + 2D Stress + Cold/Hot Demag) ===\n")
    baselines = [
        preset_ndfeb_vipm(),
        preset_pure_synrm_4p(),
        preset_ferrite_pmasynrm_6p(),
        preset_ferrite_pmasynrm_4p(),
    ]

    results = {}
    for d in baselines:
        t0 = time.perf_counter()
        r = evaluate_motor_design(d, mesh_scale=1.15, n_steps=6)
        dt = time.perf_counter() - t0
        results[d.name] = r
        print_row(r, dt)

    print("\n=== STAGE-3 MICRO-EXPERIMENTS (Structural Rib Repair + Asymmetric Barriers + 5th/7th Harmonic Injection) ===")
    # Exp 1: Strengthen B3 (6p Ferrite) center ribs to pass 14.4krpm SF >= 1.50
    b3_safe = copy.deepcopy(preset_ferrite_pmasynrm_6p())
    b3_safe.name = "M1_B3_6p_SafeRibs"
    b3_safe.barriers[0].rib_mm = 1.9
    b3_safe.barriers[1].rib_mm = 1.5
    b3_safe.barriers[2].rib_mm = 1.0
    b3_safe.barriers[0].bridge_mm = 1.25
    t0 = time.perf_counter()
    r_m1 = evaluate_motor_design(b3_safe, mesh_scale=1.15, n_steps=6)
    print_row(r_m1, time.perf_counter() - t0)

    # Exp 2: Asymmetric Barrier Wing Shift on M1 (asymmetric rotor!)
    b3_asym = copy.deepcopy(b3_safe)
    b3_asym.name = "M2_B3_6p_AsymBarriers"
    b3_asym.barriers[0].beta_L_frac -= 0.04
    b3_asym.barriers[0].beta_R_frac += 0.03
    b3_asym.barriers[1].beta_L_frac += 0.05
    b3_asym.barriers[1].beta_R_frac -= 0.04
    b3_asym.barriers[2].beta_L_frac -= 0.05
    b3_asym.barriers[2].beta_R_frac += 0.05
    t0 = time.perf_counter()
    r_m2 = evaluate_motor_design(b3_asym, mesh_scale=1.15, n_steps=6)
    print_row(r_m2, time.perf_counter() - t0)

    # Exp 3: Add 5th/7th Harmonic Current Injection on M2
    b3_harm = copy.deepcopy(b3_asym)
    b3_harm.name = "M3_B3_6p_Asym+Harmonic"
    b3_harm.k_h5 = 0.05
    b3_harm.phi_h5_rad = 1.2
    b3_harm.k_h7 = 0.03
    b3_harm.phi_h7_rad = -0.8
    t0 = time.perf_counter()
    r_m3 = evaluate_motor_design(b3_harm, mesh_scale=1.15, n_steps=6)
    print_row(r_m3, time.perf_counter() - t0)

    # Exp 4: 4p/48s (q=4!) Ferrite PMa-SynRM with strong ribs
    b4_48 = copy.deepcopy(preset_ferrite_pmasynrm_4p())
    b4_48.name = "M4_Ferrite_4p48s_Safe"
    b4_48.q_spp = 4
    b4_48.N_turns = 3.0
    b4_48.barriers[0].rib_mm = 2.2
    b4_48.barriers[1].rib_mm = 1.8
    b4_48.barriers[2].rib_mm = 1.4
    b4_48.barriers[3].rib_mm = 1.0
    t0 = time.perf_counter()
    r_m4 = evaluate_motor_design(b4_48, mesh_scale=1.15, n_steps=6)
    print_row(r_m4, time.perf_counter() - t0)

    results["M1_B3_6p_SafeRibs"] = r_m1
    results["M2_B3_6p_AsymBarriers"] = r_m2
    results["M3_B3_6p_Asym+Harmonic"] = r_m3
    results["M4_Ferrite_4p48s_Safe"] = r_m4

    with open("motor_baselines.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nSaved motor_baselines.json")


if __name__ == "__main__":
    main()
