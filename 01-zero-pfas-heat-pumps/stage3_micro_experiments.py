"""
Stage 3: Controlled Hypothesis-Verification Micro-Experiments & Profiling.

Quantitatively tests each of the 132-year physical hypotheses (H1..H6) on the
EXACT SAME locked hardware envelope, compressor physics, Silver-Bell-Ghaly
degradation, Friedel/MSH dP, fan power, and IEC 60335-2-40 500g A3 charge cap.
"""

from __future__ import annotations

import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from typing import Dict, List, Tuple

import numpy as np

from cycle_evaluator import (
    INDOOR_H,
    INDOOR_W,
    OUTDOOR_H,
    OUTDOOR_W,
    SystemDesign,
    SystemEvaluation,
    evaluate_system,
)
from microchannel_hx import CoilGeometry


def make_2row_coils(bridge: bool = False) -> Tuple[CoilGeometry, CoilGeometry]:
    """Balanced 2-row microchannel coils strictly inside the locked Al & depth caps."""
    in_coil = CoilGeometry(
        name="indoor_2row",
        width=INDOOR_W, height=INDOOR_H,
        n_rows=2, slab_depth=0.0155,
        tube_pitch=0.0105, tube_height=0.0016,
        port_dh=0.00080, fin_pitch=0.00145,
        pass_fractions=(0.55, 0.45),
        counterflow_bridge=bridge,
    )
    out_coil = CoilGeometry(
        name="outdoor_2row",
        width=OUTDOOR_W, height=OUTDOOR_H,
        n_rows=2, slab_depth=0.0170,
        tube_pitch=0.0108, tube_height=0.0016,
        port_dh=0.00085, fin_pitch=0.00150,
        pass_fractions=(0.60, 0.40),
        counterflow_bridge=bridge,
    )
    return in_coil, out_coil


# Reasonable variable-speed airflows tuned for 2-row coils
BASE_V_IN = (0.30, 0.30, 0.21, 0.19, 0.25, 0.31)
BASE_V_OUT = (0.85, 0.80, 0.52, 0.48, 0.65, 0.82)
BASE_SH = (4.5, 4.5, 4.5, 4.5, 4.5, 4.5)
BASE_SC = (5.0, 5.0, 4.0, 4.0, 4.5, 5.0)
ZERO_HOLD = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _eval_worker(args):
    design, enforce_gwp = args
    return design.label, evaluate_system(design, enforce_natural_gwp=enforce_gwp)


def run_batch(tasks: List[Tuple[SystemDesign, bool]], max_workers: int = 16) -> List[Tuple[str, SystemEvaluation]]:
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(_eval_worker, tasks))


def build_design(
    label: str,
    comps: Tuple[str, ...],
    z: Tuple[float, ...],
    topology: str = "FTVI_IHX",
    bridge: bool = False,
    comp_shift: bool = False,
    hold: Tuple[float, ...] = ZERO_HOLD,
    disp_cc: float = 34.0,
    ihx_eff: float = 0.55,
    n_rows: int = 2,
) -> SystemDesign:
    if n_rows == 1:
        ic = CoilGeometry("in_1r", INDOOR_W, INDOOR_H, 1, 0.020, 0.0095, 0.0018, 0.0010, 0.0014, (0.55, 0.45), bridge)
        oc = CoilGeometry("out_1r", OUTDOOR_W, OUTDOOR_H, 1, 0.022, 0.0095, 0.0018, 0.0011, 0.0015, (0.60, 0.40), bridge)
    else:
        ic, oc = make_2row_coils(bridge=bridge)
    return SystemDesign(
        label=label, topology=topology, ihx_eff=ihx_eff, ft_pr_bias=1.0,
        comp_shift_enabled=comp_shift, components=comps, mole_fractions=z,
        indoor_coil=ic, outdoor_coil=oc, comp_disp_cc=disp_cc,
        v_air_in_m3s=BASE_V_IN, v_air_out_m3s=BASE_V_OUT,
        superheat_K=BASE_SH, subcool_K=BASE_SC, liq_hold_frac=hold,
    )


def main():
    t_start = time.perf_counter()

    # Check 2-row aluminum mass first
    ic2, oc2 = make_2row_coils()
    print(f"2-Row Coil Al Mass: indoor={ic2.aluminum_mass:.2f}kg (cap 3.20), outdoor={oc2.aluminum_mass:.2f}kg (cap 7.50)")
    print(f"2-Row Internal Vol: indoor={ic2.internal_volume*1e6:.0f}mL, outdoor={oc2.internal_volume*1e6:.0f}mL\n")

    tasks: List[Tuple[SystemDesign, bool]] = []

    # --- H1: Pure fluids (1-row vs 2-row, no bridge vs bridge) ---
    tasks.append((build_design("H1a: R-32 1-row TXV (std)", ("R32",), (1.0,), "TXV_BASIC", False, disp_cc=20.0, n_rows=1), False))
    tasks.append((build_design("H1b: R-32 2-row FTVI+IHX", ("R32",), (1.0,), "FTVI_IHX", False, disp_cc=20.0, n_rows=2), False))
    tasks.append((build_design("H1c: Pure R-290 1-row TXV", ("Propane",), (1.0,), "TXV_BASIC", False, disp_cc=32.0, n_rows=1), True))
    tasks.append((build_design("H1d: Pure R-290 2-row FTVI+IHX (no bridge)", ("Propane",), (1.0,), "FTVI_IHX", False, disp_cc=32.0, n_rows=2), True))
    tasks.append((build_design("H1e: Pure R-290 2-row FTVI+IHX (+bridge)", ("Propane",), (1.0,), "FTVI_IHX", True, disp_cc=32.0, n_rows=2), True))
    tasks.append((build_design("H1f: Pure R-1270 2-row FTVI+IHX (+bridge)", ("Propylene",), (1.0,), "FTVI_IHX", True, disp_cc=26.0, n_rows=2), True))

    # --- H2 & H3: Binary natural zeotropes WITHOUT vs WITH counterflow bridge ---
    binaries = [
        ("C3/iC4", ("Propane", "IsoButane"), [0.15, 0.25, 0.35, 0.45], 38.0),
        ("C3/nC4", ("Propane", "n-Butane"), [0.10, 0.20, 0.30, 0.40], 39.0),
        ("C3H6/iC4", ("Propylene", "IsoButane"), [0.15, 0.25, 0.35, 0.45], 32.0),
        ("CO2/C3", ("CarbonDioxide", "Propane"), [0.04, 0.08, 0.12, 0.18], 26.0),
        ("C2/C3", ("Ethane", "Propane"), [0.06, 0.12, 0.18, 0.25], 27.0),
        ("C3/iC5", ("Propane", "Isopentane"), [0.04, 0.08, 0.12, 0.16], 37.0),
    ]
    for tag, comps, fracs, disp in binaries:
        for f2 in fracs:
            # For CO2/C3 and C2/C3, f2 is the light component (index 0); for others, heavy component (index 1)
            if comps[0] in ("CarbonDioxide", "Ethane"):
                z = (f2, 1.0 - f2)
                lbl_z = f"{int(round(f2*100))}/{int(round((1-f2)*100))}"
            else:
                z = (1.0 - f2, f2)
                lbl_z = f"{int(round((1-f2)*100))}/{int(round(f2*100))}"
            tasks.append((build_design(f"H2-NoBridge {tag} {lbl_z}", comps, z, "FTVI_IHX", False, disp_cc=disp), True))
            tasks.append((build_design(f"H3-Bridge   {tag} {lbl_z}", comps, z, "FTVI_IHX", True, disp_cc=disp), True))

    # --- H4: Topology ablation on top binary (C3H6/iC4 70/30 & C3/nC4 80/20) ---
    for topo in ("TXV_BASIC", "IHX_SUBCOOL", "FTVI_ECON", "FTVI_IHX"):
        tasks.append((build_design(f"H4-Topo {topo:11s} C3/nC4 80/20", ("Propane", "n-Butane"), (0.80, 0.20), topo, True, disp_cc=38.0), True))
        tasks.append((build_design(f"H4-Topo {topo:11s} C3H6/iC4 70/30", ("Propylene", "IsoButane"), (0.70, 0.30), topo, True, disp_cc=32.0), True))

    # --- H5: Active Composition-Shifting Receiver (Carr 1949 / Kruse 1981) ---
    # Hold heavy-rich liquid in cold heating bins (bins 4, 5) to enrich loop in volatile component!
    hold_winter = (0.0, 0.0, 0.0, 0.10, 0.22, 0.34)
    for tag, comps, z, disp in [
        ("CO2/C3 10/90", ("CarbonDioxide", "Propane"), (0.10, 0.90), 26.0),
        ("C3/nC4 75/25", ("Propane", "n-Butane"), (0.75, 0.25), 38.0),
        ("C3H6/iC4 65/35", ("Propylene", "IsoButane"), (0.65, 0.35), 32.0),
        ("CO2/C3/iC4 6/69/25", ("CarbonDioxide", "Propane", "IsoButane"), (0.06, 0.69, 0.25), 32.0),
    ]:
        tasks.append((build_design(f"H5-NoShift {tag}", comps, z, "FTVI_IHX", True, False, ZERO_HOLD, disp), True))
        tasks.append((build_design(f"H5-Shift   {tag}", comps, z, "FTVI_IHX", True, True, hold_winter, disp), True))

    # --- H6: Ternary natural blends (non-linear T(h) glide shaping) ---
    ternaries = [
        ("CO2/C3/iC4 4/71/25", ("CarbonDioxide", "Propane", "IsoButane"), (0.04, 0.71, 0.25), 33.0),
        ("CO2/C3/iC4 8/62/30", ("CarbonDioxide", "Propane", "IsoButane"), (0.08, 0.62, 0.30), 31.0),
        ("CO2/C3/nC4 5/75/20", ("CarbonDioxide", "Propane", "n-Butane"), (0.05, 0.75, 0.20), 33.0),
        ("C2/C3/iC4 8/67/25", ("Ethane", "Propane", "IsoButane"), (0.08, 0.67, 0.25), 32.0),
        ("C2/C3/nC4 10/70/20", ("Ethane", "Propane", "n-Butane"), (0.10, 0.70, 0.20), 32.0),
        ("C3H6/C3/iC4 40/35/25", ("Propylene", "Propane", "IsoButane"), (0.40, 0.35, 0.25), 33.0),
        ("C3/iC4/nC4 72/18/10", ("Propane", "IsoButane", "n-Butane"), (0.72, 0.18, 0.10), 38.0),
        ("C3/iC4/iC5 78/16/6", ("Propane", "IsoButane", "Isopentane"), (0.78, 0.16, 0.06), 38.0),
    ]
    for tag, comps, z, disp in ternaries:
        tasks.append((build_design(f"H6-Ternary {tag}", comps, z, "FTVI_IHX", True, False, ZERO_HOLD, disp), True))

    print(f"Running {len(tasks)} controlled micro-experiments across 16 Apple Silicon cores...")
    results = run_batch(tasks, max_workers=16)
    dt = time.perf_counter() - t_start
    print(f"Completed {len(tasks)} full 6-bin HEOS evaluations in {dt:.2f}s ({len(tasks)/dt:.1f} evals/s)!\n")

    print(f"{'Experiment':40s} | {'Feas':4s} | {'ASPF':6s} | {'SEER':6s} | {'SCOP':6s} | {'46C':5s} | {'-7C':5s} | {'Glide':5s} | {'Chrg':5s} | {'Tdis':5s} | Notes")
    print("-" * 135)
    for lbl, ev in results:
        note = "" if ev.feasible else "; ".join(ev.violations[:2])
        print(
            f"{lbl:40s} | {'YES' if ev.feasible else 'NO ':4s} | "
            f"{ev.ASPF:6.3f} | {ev.SEER_cool:6.3f} | {ev.SCOP_heat:6.3f} | "
            f"{ev.COP_46C_T3:5.2f} | {ev.COP_m7C_cold:5.2f} | "
            f"{ev.max_glide_K:4.1f}K | {ev.flammable_charge_g:4.0f}g | {ev.max_T_dis_C:4.0f}C | {note}"
        )

    # Print ranked Top-10 Feasible
    feas = [(lbl, ev) for lbl, ev in results if ev.feasible]
    feas.sort(key=lambda x: x[1].ASPF, reverse=True)
    print("\n=== STAGE 3 TOP 10 FEASIBLE MICRO-EXPERIMENTS ===")
    for i, (lbl, ev) in enumerate(feas[:10], 1):
        print(
            f"  #{i:2d} {lbl:38s}: ASPF={ev.ASPF:.3f} (SEER={ev.SEER_cool:.3f}, SCOP={ev.SCOP_heat:.3f}, "
            f"46C={ev.COP_46C_T3:.2f}, -7C={ev.COP_m7C_cold:.2f}, glide={ev.max_glide_K:.1f}K, chrg={ev.flammable_charge_g:.0f}g)"
        )


if __name__ == "__main__":
    main()
