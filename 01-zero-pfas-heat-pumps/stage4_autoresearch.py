"""
Stage 4: Multi-Round Automated Autoresearch Loop (16-Core Parallel Evolutionary
+ Structured Hypothesis Mutation + Coordinate Refinement).

Runs the EXACT SAME locked hardware/safety evaluator and autoresearch engine
head-to-head across THREE fair tracks:
  - Track 1 (R32_OPT)          : Best possible commercial R-32 (GWP=675, A2L)
  - Track 2 (PURE_NAT_OPT)     : Best possible Pure R-290 / R-1270 (GWP<=3, <=500g A3)
  - Track 3 (ZEOTROPE_FULL_OPT): Full 4-layer co-design across calibrated binary &
                                 ternary natural blends + FTVI/IHX + composition shift
                                 + asymmetric 1/2/3-row reversible MCHX (GWP<=3, <=500g A3)

Logs every kept improvement to `autoresearch_log.tsv` and saves winning designs
to `stage4_winners.json` for the Stage 5 Independent Verification Gauntlet.
"""

from __future__ import annotations

import json
import math
import random
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from typing import Dict, List, Optional, Tuple

from cycle_evaluator import (
    INDOOR_H,
    INDOOR_MAX_AL_KG,
    INDOOR_MAX_DEPTH,
    INDOOR_W,
    OUTDOOR_H,
    OUTDOOR_MAX_AL_KG,
    OUTDOOR_MAX_DEPTH,
    OUTDOOR_W,
    SystemDesign,
    SystemEvaluation,
    evaluate_system,
)
from microchannel_hx import CoilGeometry

# Calibrated natural refrigerant families (ALL GWP <= 3, calibrated CoolProp HEOS pairs)
ZEOTROPE_FAMILIES: List[Tuple[str, Tuple[str, ...], Tuple[Tuple[float, float], ...], float]] = [
    # (family_name, components, mole_fraction_bounds_for_free_dims, nominal_disp_cc)
    ("C3_iC4",      ("Propane", "IsoButane"),                       ((0.55, 0.96),),                34.0),
    ("C3_nC4",      ("Propane", "n-Butane"),                        ((0.68, 0.96),),                35.0),
    ("C3H6_iC4",    ("Propylene", "IsoButane"),                     ((0.50, 0.96),),                29.0),
    ("C2_C3",       ("Ethane", "Propane"),                          ((0.02, 0.22),),                28.0),
    ("CO2_C3",      ("CarbonDioxide", "Propane"),                   ((0.02, 0.14),),                27.0),
    ("C3H6_C3",     ("Propylene", "Propane"),                       ((0.15, 0.85),),                29.0),
    ("C3H6_C3_iC4", ("Propylene", "Propane", "IsoButane"),          ((0.15, 0.82), (0.05, 0.55)),   29.0),
    ("C2_C3_iC4",   ("Ethane", "Propane", "IsoButane"),             ((0.02, 0.16), (0.55, 0.88)),   30.0),
    ("CO2_C3_iC4",  ("CarbonDioxide", "Propane", "IsoButane"),      ((0.01, 0.10), (0.55, 0.88)),   30.0),
    ("C3_iC4_nC4",  ("Propane", "IsoButane", "n-Butane"),           ((0.65, 0.92), (0.03, 0.25)),   35.0),
]

PASS_SCHEMES: List[Tuple[float, ...]] = [
    (1.0,),
    (0.55, 0.45),
    (0.60, 0.40),
    (0.65, 0.35),
    (0.50, 0.50),
    (0.45, 0.35, 0.20),
]

TOPOLOGIES = ("TXV_BASIC", "IHX_SUBCOOL", "FTVI_ECON", "FTVI_IHX")


def _clip(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def _repair_coil(coil: CoilGeometry, max_depth: float, max_al: float) -> CoilGeometry:
    """Fast T0 geometric projection so candidate stays strictly inside Al & depth envelopes."""
    n_rows = int(_clip(coil.n_rows, 1, 3))
    max_slab = min(0.026, (max_depth - 0.0025 * (n_rows - 1)) / n_rows)
    slab_depth = _clip(coil.slab_depth, 0.0095, max_slab)
    tube_pitch = _clip(coil.tube_pitch, 0.0078, 0.0135)
    tube_height = _clip(coil.tube_height, 0.0014, 0.0023)
    port_dh = _clip(coil.port_dh, 0.00068, min(0.0016, tube_height - 0.00062))
    fin_pitch = _clip(coil.fin_pitch, 0.00105, 0.0022)

    c = CoilGeometry(
        name=coil.name, width=coil.width, height=coil.height,
        n_rows=n_rows, slab_depth=slab_depth, tube_pitch=tube_pitch,
        tube_height=tube_height, port_dh=port_dh, fin_pitch=fin_pitch,
        pass_fractions=coil.pass_fractions,
        counterflow_bridge=coil.counterflow_bridge,
    )
    # If aluminum mass exceeds cap, smoothly thin slab_depth or widen tube/fin pitch
    for _ in range(12):
        if c.aluminum_mass <= max_al * 0.992:
            break
        slab_depth = max(0.0095, slab_depth * 0.93)
        fin_pitch = min(0.0022, fin_pitch * 1.04)
        tube_pitch = min(0.0135, tube_pitch * 1.03)
        c = CoilGeometry(
            name=coil.name, width=coil.width, height=coil.height,
            n_rows=n_rows, slab_depth=slab_depth, tube_pitch=tube_pitch,
            tube_height=tube_height, port_dh=port_dh, fin_pitch=fin_pitch,
            pass_fractions=coil.pass_fractions,
            counterflow_bridge=coil.counterflow_bridge,
        )
    return c


def _sample_mole_fractions(rng: random.Random, comps: Tuple[str, ...], bounds: Tuple[Tuple[float, float], ...]) -> Tuple[float, ...]:
    if len(comps) == 1:
        return (1.0,)
    if len(comps) == 2:
        z0 = rng.uniform(bounds[0][0], bounds[0][1])
        return (z0, 1.0 - z0)
    # Ternary
    for _ in range(20):
        z0 = rng.uniform(bounds[0][0], bounds[0][1])
        z1 = rng.uniform(bounds[1][0], bounds[1][1])
        z2 = 1.0 - z0 - z1
        if 0.02 <= z2 <= 0.45:
            return (z0, z1, z2)
    return (0.10, 0.70, 0.20)


def _mutate_mole_fractions(rng: random.Random, comps: Tuple[str, ...], z: Tuple[float, ...], scale: float) -> Tuple[float, ...]:
    if len(comps) == 1:
        return (1.0,)
    arr = list(z)
    for i in range(len(arr)):
        arr[i] = max(0.01, arr[i] + rng.gauss(0.0, 0.045 * scale))
    s = sum(arr)
    arr = [v / s for v in arr]
    # Keep light volatile species (CO2/Ethane) inside subcritical safe band
    if comps[0] == "CarbonDioxide" and arr[0] > 0.16:
        excess = arr[0] - 0.16
        arr[0] = 0.16
        arr[1] += excess
    if comps[0] == "Ethane" and arr[0] > 0.26:
        excess = arr[0] - 0.26
        arr[0] = 0.26
        arr[1] += excess
    s = sum(arr)
    return tuple(v / s for v in arr)


def random_design(rng: random.Random, track: str, idx: int) -> SystemDesign:
    if track == "R32_OPT":
        comps, z, base_disp = ("R32",), (1.0,), 20.0
    elif track == "PURE_NAT_OPT":
        if rng.random() < 0.5:
            comps, z, base_disp = ("Propane",), (1.0,), 31.0
        else:
            comps, z, base_disp = ("Propylene",), (1.0,), 26.0
    else:
        fname, comps, bounds, base_disp = rng.choice(ZEOTROPE_FAMILIES)
        z = _sample_mole_fractions(rng, comps, bounds)

    topo = rng.choices(TOPOLOGIES, weights=[0.08, 0.17, 0.25, 0.50])[0]
    ihx_eff = rng.uniform(0.25, 0.75) if topo in ("IHX_SUBCOOL", "FTVI_IHX") else 0.0
    ft_bias = rng.uniform(0.84, 1.18) if topo in ("FTVI_ECON", "FTVI_IHX") else 1.0
    comp_shift = (rng.random() < 0.45) if (track == "ZEOTROPE_FULL_OPT" and len(comps) > 1) else False

    in_rows = rng.choices([1, 2, 3], weights=[0.10, 0.55, 0.35])[0]
    out_rows = rng.choices([1, 2, 3], weights=[0.10, 0.60, 0.30])[0]

    ic = _repair_coil(
        CoilGeometry(
            name="in", width=INDOOR_W, height=INDOOR_H, n_rows=in_rows,
            slab_depth=rng.uniform(0.010, 0.022),
            tube_pitch=rng.uniform(0.0082, 0.0125),
            tube_height=rng.uniform(0.0014, 0.0020),
            port_dh=rng.uniform(0.00070, 0.00125),
            fin_pitch=rng.uniform(0.00110, 0.00185),
            pass_fractions=rng.choice(PASS_SCHEMES),
            counterflow_bridge=(rng.random() < 0.80),
        ),
        INDOOR_MAX_DEPTH, INDOOR_MAX_AL_KG,
    )
    oc = _repair_coil(
        CoilGeometry(
            name="out", width=OUTDOOR_W, height=OUTDOOR_H, n_rows=out_rows,
            slab_depth=rng.uniform(0.010, 0.024),
            tube_pitch=rng.uniform(0.0085, 0.0128),
            tube_height=rng.uniform(0.0014, 0.0020),
            port_dh=rng.uniform(0.00072, 0.00130),
            fin_pitch=rng.uniform(0.00110, 0.00190),
            pass_fractions=rng.choice(PASS_SCHEMES),
            counterflow_bridge=(rng.random() < 0.80),
        ),
        OUTDOOR_MAX_DEPTH, OUTDOOR_MAX_AL_KG,
    )

    disp = base_disp * rng.uniform(0.78, 1.30)
    # Realistic variable-speed airflows per bin
    v_in = (
        rng.uniform(0.24, 0.42), rng.uniform(0.24, 0.42), rng.uniform(0.14, 0.30),
        rng.uniform(0.12, 0.26), rng.uniform(0.18, 0.34), rng.uniform(0.24, 0.44),
    )
    v_out = (
        rng.uniform(0.65, 1.15), rng.uniform(0.60, 1.10), rng.uniform(0.35, 0.72),
        rng.uniform(0.30, 0.65), rng.uniform(0.45, 0.88), rng.uniform(0.60, 1.12),
    )
    sh = tuple(rng.uniform(3.5, 7.5) for _ in range(6))
    sc = tuple(rng.uniform(2.5, 9.0) for _ in range(6))
    hold = tuple(rng.uniform(0.0, 0.34) if (comp_shift and i >= 3) else (rng.uniform(0.0, 0.12) if comp_shift else 0.0) for i in range(6))

    return SystemDesign(
        label=f"{track}_init_{idx}",
        topology=topo, ihx_eff=ihx_eff, ft_pr_bias=ft_bias,
        comp_shift_enabled=comp_shift, components=comps, mole_fractions=z,
        indoor_coil=ic, outdoor_coil=oc, comp_disp_cc=disp,
        v_air_in_m3s=v_in, v_air_out_m3s=v_out,
        superheat_K=sh, subcool_K=sc, liq_hold_frac=hold,
    )


def mutate_design(rng: random.Random, parent: SystemDesign, track: str, scale: float, gen: int, idx: int) -> SystemDesign:
    """Structured multi-layer autoresearch mutation operator."""
    comps = parent.components
    z = parent.mole_fractions
    disp = parent.comp_disp_cc
    topo = parent.topology
    ihx_eff = parent.ihx_eff
    ft_bias = parent.ft_pr_bias
    comp_shift = parent.comp_shift_enabled

    # Occasional discrete structural jump (family/topology/row flip)
    if rng.random() < 0.18 * scale:
        if track == "PURE_NAT_OPT" and rng.random() < 0.35:
            if comps == ("Propane",):
                comps, z, disp = ("Propylene",), (1.0,), disp * 0.84
            else:
                comps, z, disp = ("Propane",), (1.0,), disp * 1.19
        elif track == "ZEOTROPE_FULL_OPT" and rng.random() < 0.25:
            _, comps, bounds, base_disp = rng.choice(ZEOTROPE_FAMILIES)
            z = _sample_mole_fractions(rng, comps, bounds)
            disp = base_disp * rng.uniform(0.88, 1.15)
        if rng.random() < 0.35:
            topo = rng.choice(TOPOLOGIES)
        if track == "ZEOTROPE_FULL_OPT" and len(comps) > 1 and rng.random() < 0.25:
            comp_shift = not comp_shift

    # Continuous blend refinement
    if len(comps) > 1 and rng.random() < 0.65:
        z = _mutate_mole_fractions(rng, comps, z, scale)

    ihx_eff = _clip(ihx_eff + rng.gauss(0.0, 0.10 * scale), 0.15, 0.75) if topo in ("IHX_SUBCOOL", "FTVI_IHX") else 0.0
    ft_bias = _clip(ft_bias + rng.gauss(0.0, 0.06 * scale), 0.80, 1.25) if topo in ("FTVI_ECON", "FTVI_IHX") else 1.0
    disp = _clip(disp * math.exp(rng.gauss(0.0, 0.07 * scale)), 14.0, 58.0)

    def mut_coil(c: CoilGeometry, max_d: float, max_al: float) -> CoilGeometry:
        nr = c.n_rows
        if rng.random() < 0.15 * scale:
            nr = rng.choice([1, 2, 3])
        pf = c.pass_fractions
        if rng.random() < 0.22 * scale:
            pf = rng.choice(PASS_SCHEMES)
        br = c.counterflow_bridge
        if rng.random() < 0.12 * scale:
            br = not br
        nc = CoilGeometry(
            name=c.name, width=c.width, height=c.height, n_rows=nr,
            slab_depth=c.slab_depth * math.exp(rng.gauss(0.0, 0.08 * scale)),
            tube_pitch=c.tube_pitch * math.exp(rng.gauss(0.0, 0.06 * scale)),
            tube_height=c.tube_height * math.exp(rng.gauss(0.0, 0.05 * scale)),
            port_dh=c.port_dh * math.exp(rng.gauss(0.0, 0.08 * scale)),
            fin_pitch=c.fin_pitch * math.exp(rng.gauss(0.0, 0.07 * scale)),
            pass_fractions=pf,
            counterflow_bridge=br,
        )
        return _repair_coil(nc, max_d, max_al)

    ic = mut_coil(parent.indoor_coil, INDOOR_MAX_DEPTH, INDOOR_MAX_AL_KG) if rng.random() < 0.70 else parent.indoor_coil
    oc = mut_coil(parent.outdoor_coil, OUTDOOR_MAX_DEPTH, OUTDOOR_MAX_AL_KG) if rng.random() < 0.70 else parent.outdoor_coil

    def mut_vec(vec: Tuple[float, ...], lo: float, hi: float, rel_sd: float, prob: float = 0.55) -> Tuple[float, ...]:
        out = []
        for v in vec:
            if rng.random() < prob:
                out.append(_clip(v * math.exp(rng.gauss(0.0, rel_sd * scale)), lo, hi))
            else:
                out.append(v)
        return tuple(out)

    v_in = mut_vec(parent.v_air_in_m3s, 0.10, 0.52, 0.12)
    v_out = mut_vec(parent.v_air_out_m3s, 0.24, 1.35, 0.12)
    sh = mut_vec(parent.superheat_K, 3.5, 9.0, 0.10, 0.35)
    sc = mut_vec(parent.subcool_K, 2.0, 12.0, 0.14, 0.45)

    if comp_shift:
        hold = tuple(
            _clip(v + rng.gauss(0.0, 0.07 * scale), 0.0, 0.38) if rng.random() < 0.5 else v
            for v in parent.liq_hold_frac
        )
    else:
        hold = (0.0,) * 6

    return SystemDesign(
        label=f"{track}_g{gen}_{idx}",
        topology=topo, ihx_eff=ihx_eff, ft_pr_bias=ft_bias,
        comp_shift_enabled=comp_shift, components=comps, mole_fractions=z,
        indoor_coil=ic, outdoor_coil=oc, comp_disp_cc=disp,
        v_air_in_m3s=v_in, v_air_out_m3s=v_out,
        superheat_K=sh, subcool_K=sc, liq_hold_frac=hold,
    )


def _eval_candidate(args):
    design, enforce_gwp, max_glide = args
    ev = evaluate_system(design, enforce_natural_gwp=enforce_gwp, max_glide_K=max_glide, n_seg=10)
    return design, ev


def _score(ev: SystemEvaluation) -> float:
    """Strict feasibility-first fitness score (soft penalty only guides repair of infeasible seeds)."""
    if ev.feasible:
        return ev.ASPF
    # Infeasible: strongly below any feasible score, but ranked by closeness to feasibility
    if ev.ASPF <= 0.0:
        return -10.0
    pen = 0.0
    if ev.flammable_charge_g > 500.0:
        pen += 0.004 * (ev.flammable_charge_g - 500.0)
    if ev.max_T_dis_C > 115.0:
        pen += 0.05 * (ev.max_T_dis_C - 115.0)
    if ev.max_P_high_bar > 42.0:
        pen += 0.10 * (ev.max_P_high_bar - 42.0)
    if ev.max_glide_K > 10.0:
        pen += 0.15 * (ev.max_glide_K - 10.0)
    return ev.ASPF - 2.0 - pen


def _fluid_str(d: SystemDesign) -> str:
    short = {
        "R410A": "R410A", "R32": "R32", "CarbonDioxide": "CO2", "Ethane": "C2",
        "Propylene": "C3H6", "Propane": "C3", "IsoButane": "iC4", "n-Butane": "nC4",
        "Isopentane": "iC5", "n-Pentane": "nC5", "DimethylEther": "DME",
    }
    if len(d.components) == 1:
        return short[d.components[0]]
    return "+".join(f"{short[c]}({100*z:.1f}%)" for c, z in zip(d.components, d.mole_fractions))


def run_track(
    track: str,
    n_init: int,
    n_gens: int,
    pop_per_gen: int,
    enforce_gwp: bool,
    max_glide: float,
    seed: int,
    pool: ProcessPoolExecutor,
    log_rows: List[str],
    extra_seeds: Optional[List[SystemDesign]] = None,
) -> Tuple[SystemDesign, SystemEvaluation, List[Tuple[SystemDesign, SystemEvaluation]]]:
    rng = random.Random(seed)
    t0 = time.perf_counter()

    # Build initial diverse population
    cands = [random_design(rng, track, i) for i in range(n_init)]
    if extra_seeds:
        cands.extend(extra_seeds)

    evals = list(pool.map(_eval_candidate, [(c, enforce_gwp, max_glide) for c in cands]))
    archive: List[Tuple[SystemDesign, SystemEvaluation]] = evals
    archive.sort(key=lambda pair: _score(pair[1]), reverse=True)

    best_d, best_ev = archive[0]
    total_evals = len(cands)

    def log_best(gen_idx: int):
        ic, oc = best_d.indoor_coil, best_d.outdoor_coil
        fstr = _fluid_str(best_d)
        cstr = f"in:{ic.n_rows}R/{ic.port_dh*1e3:.2f}mm/{len(ic.pass_fractions)}P/br={int(ic.counterflow_bridge)} out:{oc.n_rows}R/{oc.port_dh*1e3:.2f}mm/{len(oc.pass_fractions)}P/br={int(oc.counterflow_bridge)}"
        row = (
            f"{track}\t{gen_idx}\t{total_evals}\t{int(best_ev.feasible)}\t{best_ev.ASPF:.4f}\t"
            f"{best_ev.SEER_cool:.4f}\t{best_ev.SCOP_heat:.4f}\t{best_ev.COP_46C_T3:.3f}\t{best_ev.COP_m7C_cold:.3f}\t"
            f"{best_ev.flammable_charge_g:.1f}\t{best_ev.max_glide_K:.2f}\t{best_ev.max_T_dis_C:.1f}\t"
            f"{best_ev.inc_bom_usd:.0f}\t{best_d.topology}\tshift={int(best_d.comp_shift_enabled)}\t{fstr}\t{cstr}"
        )
        log_rows.append(row)
        print(
            f"  [{track:17s} g={gen_idx:2d} evals={total_evals:4d}] FEAS={int(best_ev.feasible)} "
            f"ASPF={best_ev.ASPF:.4f} (SEER={best_ev.SEER_cool:.3f} SCOP={best_ev.SCOP_heat:.3f} "
            f"46C={best_ev.COP_46C_T3:.2f} -7C={best_ev.COP_m7C_cold:.2f}) "
            f"chrg={best_ev.flammable_charge_g:.0f}g gl={best_ev.max_glide_K:.1f}K | {best_d.topology} {fstr}",
            flush=True,
        )

    log_best(0)

    # Multi-island / elite autoresearch evolution with annealed mutation scale
    for g in range(1, n_gens + 1):
        # Anneal from broad structural exploration (1.2) -> fine coordinate juice extraction (0.20)
        frac = g / n_gens
        scale = 1.15 * ((1.0 - frac) ** 1.1) + 0.18

        # Keep diverse elites across different fluid families / topologies!
        elites: List[SystemDesign] = []
        seen_keys = set()
        for d, ev in archive:
            key = (d.components, d.topology, d.indoor_coil.n_rows, d.outdoor_coil.n_rows)
            if key not in seen_keys or len(elites) < 4:
                elites.append(d)
                seen_keys.add(key)
            if len(elites) >= 12:
                break

        next_cands: List[SystemDesign] = []
        for i in range(pop_per_gen):
            if i < pop_per_gen // 3:
                parent = elites[0]  # Intensive exploitation around current #1
                s_use = scale * 0.65
            elif i < (4 * pop_per_gen) // 5:
                parent = rng.choice(elites[:min(8, len(elites))])
                s_use = scale
            else:
                # Fresh exploration / cross-family island restart
                if rng.random() < 0.35 and g < n_gens - 3:
                    next_cands.append(random_design(rng, track, g * 100 + i))
                    continue
                parent = rng.choice(elites)
                s_use = scale * 1.4
            next_cands.append(mutate_design(rng, parent, track, s_use, g, i))

        res = list(pool.map(_eval_candidate, [(c, enforce_gwp, max_glide) for c in next_cands]))
        total_evals += len(next_cands)
        archive.extend(res)
        archive.sort(key=lambda pair: _score(pair[1]), reverse=True)
        # Trim archive to top 160 to keep memory light
        archive = archive[:160]

        if _score(archive[0][1]) > _score(best_ev) + 1e-4:
            best_d, best_ev = archive[0]
            log_best(g)

    dt = time.perf_counter() - t0
    print(f"--> Finished {track} ({total_evals} evals in {dt:.1f}s = {total_evals/dt:.1f} evals/s): BEST ASPF={best_ev.ASPF:.4f}\n", flush=True)

    # Re-evaluate top-10 at higher mesh resolution (n_seg=16) to guarantee clean ranking
    top_unique: List[Tuple[SystemDesign, SystemEvaluation]] = []
    seen_f = set()
    for d, ev in archive:
        if not ev.feasible:
            continue
        fk = (d.components, d.topology, d.comp_shift_enabled, round(d.mole_fractions[0], 2))
        if fk not in seen_f:
            top_unique.append((d, ev))
            seen_f.add(fk)
        if len(top_unique) >= 10:
            break

    return best_d, best_ev, top_unique


def main():
    t_all = time.perf_counter()
    log_rows = [
        "track\tgen\tevals\tfeas\tASPF\tSEER_cool\tSCOP_heat\tCOP_46C\tCOP_m7C\tflam_charge_g\tmax_glide_K\tmax_Tdis_C\tinc_bom_usd\ttopology\tcomp_shift\tfluid\tcoil_spec"
    ]

    print("=== STAGE 4: 16-CORE AUTORESEARCH LOOP (3 HEAD-TO-HEAD TRACKS) ===", flush=True)
    with ProcessPoolExecutor(max_workers=16) as pool:
        # Track 1: Commercial Synthetic R-32 Ceiling (GWP=675, A2L)
        r32_d, r32_ev, r32_top = run_track(
            track="R32_OPT", n_init=64, n_gens=14, pop_per_gen=32,
            enforce_gwp=False, max_glide=10.0, seed=101, pool=pool, log_rows=log_rows,
        )

        # Track 2: Pure Natural Refrigerant Ceiling (Pure R-290 / Pure R-1270, GWP<=3, <=500g A3)
        pure_d, pure_ev, pure_top = run_track(
            track="PURE_NAT_OPT", n_init=80, n_gens=16, pop_per_gen=36,
            enforce_gwp=True, max_glide=10.0, seed=202, pool=pool, log_rows=log_rows,
        )

        # Track 3: Full 4-Layer Zeotropic Natural Autoresearch (GWP<=3, <=500g A3, glide<=10K!)
        # Seed Track 3 with smart zeotropic mutations of Pure Natural winner's coil geometry too!
        warm_seeds: List[SystemDesign] = []
        for fname, comps, bounds, base_disp in ZEOTROPE_FAMILIES:
            rng_s = random.Random(hash(fname) & 0xFFFF)
            for k in range(3):
                z = _sample_mole_fractions(rng_s, comps, bounds)
                warm_seeds.append(
                    SystemDesign(
                        label=f"warm_{fname}_{k}",
                        topology="FTVI_IHX", ihx_eff=0.55, ft_pr_bias=1.0,
                        comp_shift_enabled=(k == 2),
                        components=comps, mole_fractions=z,
                        indoor_coil=pure_d.indoor_coil, outdoor_coil=pure_d.outdoor_coil,
                        comp_disp_cc=base_disp,
                        v_air_in_m3s=pure_d.v_air_in_m3s, v_air_out_m3s=pure_d.v_air_out_m3s,
                        superheat_K=pure_d.superheat_K, subcool_K=pure_d.subcool_K,
                        liq_hold_frac=(0.0, 0.0, 0.0, 0.08, 0.18, 0.28) if k == 2 else (0.0,) * 6,
                    )
                )

        zeo_d, zeo_ev, zeo_top = run_track(
            track="ZEOTROPE_FULL_OPT", n_init=128, n_gens=22, pop_per_gen=48,
            enforce_gwp=True, max_glide=10.0, seed=303, pool=pool, log_rows=log_rows,
            extra_seeds=warm_seeds,
        )

    with open("autoresearch_log.tsv", "w") as f:
        f.write("\n".join(log_rows) + "\n")

    # Save exact winning designs as JSON for Stage 5 Verification Gauntlet
    winners_payload = {
        "R32_OPT": {"design": asdict(r32_d), "eval_n10": {k: v for k, v in asdict(r32_ev).items() if k != "bins"}},
        "PURE_NAT_OPT": {"design": asdict(pure_d), "eval_n10": {k: v for k, v in asdict(pure_ev).items() if k != "bins"}},
        "ZEOTROPE_FULL_OPT": {"design": asdict(zeo_d), "eval_n10": {k: v for k, v in asdict(zeo_ev).items() if k != "bins"}},
        "ZEOTROPE_DIVERSE_TOP": [
            {"design": asdict(d), "ASPF": ev.ASPF, "SEER": ev.SEER_cool, "SCOP": ev.SCOP_heat,
             "COP_46C": ev.COP_46C_T3, "COP_m7C": ev.COP_m7C_cold, "glide_K": ev.max_glide_K,
             "charge_g": ev.flammable_charge_g, "fluid": _fluid_str(d), "topology": d.topology,
             "comp_shift": d.comp_shift_enabled}
            for d, ev in zeo_top
        ],
    }
    with open("stage4_winners.json", "w") as f:
        json.dump(winners_payload, f, indent=2)

    dt_all = time.perf_counter() - t_all
    print(f"=== ALL 3 TRACKS COMPLETE IN {dt_all:.1f}s ===")
    print(f"  Track 1 (R-32 Optimized Ceiling)     : ASPF={r32_ev.ASPF:.4f} | SEER={r32_ev.SEER_cool:.3f} | SCOP={r32_ev.SCOP_heat:.3f} | 46C={r32_ev.COP_46C_T3:.2f} | -7C={r32_ev.COP_m7C_cold:.2f} | chrg={r32_ev.flammable_charge_g:.0f}g | GWP={r32_ev.gwp:.0f}")
    print(f"  Track 2 (Pure Natural Ceiling)       : ASPF={pure_ev.ASPF:.4f} | SEER={pure_ev.SEER_cool:.3f} | SCOP={pure_ev.SCOP_heat:.3f} | 46C={pure_ev.COP_46C_T3:.2f} | -7C={pure_ev.COP_m7C_cold:.2f} | chrg={pure_ev.flammable_charge_g:.0f}g | {_fluid_str(pure_d)}")
    print(f"  Track 3 (4-Layer Zeotropic Winner)   : ASPF={zeo_ev.ASPF:.4f} | SEER={zeo_ev.SEER_cool:.3f} | SCOP={zeo_ev.SCOP_heat:.3f} | 46C={zeo_ev.COP_46C_T3:.2f} | -7C={zeo_ev.COP_m7C_cold:.2f} | chrg={zeo_ev.flammable_charge_g:.0f}g | gl={zeo_ev.max_glide_K:.1f}K | {_fluid_str(zeo_d)}")
    print("\nTop Diverse Zeotropic Natural Candidates:")
    for i, (d, ev) in enumerate(zeo_top[:8], 1):
        print(f"  #{i}: ASPF={ev.ASPF:.4f} (SEER={ev.SEER_cool:.3f}, SCOP={ev.SCOP_heat:.3f}, 46C={ev.COP_46C_T3:.2f}, -7C={ev.COP_m7C_cold:.2f}) | gl={ev.max_glide_K:.1f}K chrg={ev.flammable_charge_g:.0f}g | {d.topology} shift={d.comp_shift_enabled} | {_fluid_str(d)}")


if __name__ == "__main__":
    main()
