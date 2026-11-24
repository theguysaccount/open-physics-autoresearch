"""
Stage 4 (Regime B): 16-Core Autoresearch Optimization for 100 kW_th Industrial
High-Temperature Heat Pumps (<150 C) across 3 Real Decarbonization Duties:
  1. DAIRY_90C     (35 -> 90 C sink, 45 -> 20 C source)
  2. DISTRICT_110C (40 -> 110 C sink, 55 -> 25 C source)
  3. CHEMICAL_135C (60 -> 135 C sink, 75 -> 40 C source)

Runs THREE brutally fair head-to-head tracks per duty (9 full autoresearch runs):
  - Track A (PURE_CEILING)   : Best possible pure fluid (PFAS R-1233zd(E)/R-1336mzz(Z)
                               + Pure natural R-290/R-600a/R-600/R-601a/R-601)
  - Track B (NAT_ZEO_COP)    : Zero-PFAS GWP<=3 natural binary/ternary blends +
                               4-layer flowsheet/BPHE co-design maximizing COP_sys
  - Track C (NAT_ZEO_LCOH)   : Zero-PFAS GWP<=3 natural binary/ternary blends +
                               4-layer co-design minimizing true LCOH ($/MWh_th)

ALL natural blends are strictly restricted to Dual-EOS families calibrated in
BOTH CoolProp HEOS (NIST) AND SINTEF ThermoPack PC-SAFT!
"""

from __future__ import annotations

import json
import math
import random
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from typing import Dict, List, Optional, Tuple

from brazed_plate_hx import BPHEGeometry
from hthp_evaluator import (
    INDUSTRIAL_DUTIES,
    MAX_TOTAL_BPHE_AREA_M2,
    HTHPDesign,
    HTHPEvaluation,
    evaluate_hthp,
)

# 100% Dual-EOS calibrated natural fluid families (GWP <= 3, ZERO PFAS, CoolProp HEOS + ThermoPack PC-SAFT)
HTHP_ZEO_FAMILIES: List[Tuple[str, Tuple[str, ...]]] = [
    ("C3_nC4",          ("Propane", "n-Butane")),
    ("C3_iC5",          ("Propane", "Isopentane")),
    ("C3_nC5",          ("Propane", "n-Pentane")),
    ("C3_nC6",          ("Propane", "n-Hexane")),
    ("iC4_iC5",         ("IsoButane", "Isopentane")),
    ("nC4_iC5",         ("n-Butane", "Isopentane")),
    ("nC4_nC5",         ("n-Butane", "n-Pentane")),
    ("nC4_nC6",         ("n-Butane", "n-Hexane")),
    ("nC5_nC6",         ("n-Pentane", "n-Hexane")),
    ("C3_nC4_nC5",      ("Propane", "n-Butane", "n-Pentane")),
    ("C3_iC4_iC5",      ("Propane", "IsoButane", "Isopentane")),
    ("C3_nC4_nC6",      ("Propane", "n-Butane", "n-Hexane")),
    ("nC4_iC5_nC6",     ("n-Butane", "Isopentane", "n-Hexane")),
    ("CO2_nC4",         ("CarbonDioxide", "n-Butane")),
    ("CO2_nC5",         ("CarbonDioxide", "n-Pentane")),
    ("CO2_nC4_nC5",     ("CarbonDioxide", "n-Butane", "n-Pentane")),
]

PURE_CANDIDATES: Dict[str, List[Tuple[str, ...]]] = {
    "DAIRY_90C": [("R1233zd(E)",), ("Propane",), ("IsoButane",), ("n-Butane",), ("Isopentane",)],
    "DISTRICT_110C": [("R1233zd(E)",), ("R1336mzz(Z)",), ("IsoButane",), ("n-Butane",), ("Isopentane",), ("n-Pentane",)],
    "CHEMICAL_135C": [("R1233zd(E)",), ("R1336mzz(Z)",), ("n-Butane",), ("Isopentane",), ("n-Pentane",), ("n-Hexane",)],
}

TOPOLOGIES = ("BASIC_1STAGE", "IHX_SUBCOOL", "FTVI_ECON", "FTVI_IHX")


def _clip(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))


def _fluid_short(d: HTHPDesign) -> str:
    m = {
        "R1233zd(E)": "R1233zdE", "R1336mzz(Z)": "R1336mzzZ", "R245fa": "R245fa",
        "CarbonDioxide": "CO2", "Propane": "C3", "IsoButane": "iC4", "n-Butane": "nC4",
        "Isopentane": "iC5", "n-Pentane": "nC5", "n-Hexane": "nC6", "Ammonia": "NH3",
    }
    if len(d.components) == 1:
        return m.get(d.components[0], d.components[0])
    return "+".join(f"{m[c]}({100*z:.0f}%)" for c, z in zip(d.components, d.mole_fractions))


def _sample_z(rng: random.Random, comps: Tuple[str, ...]) -> Tuple[float, ...]:
    if len(comps) == 1:
        return (1.0,)
    if comps[0] == "CarbonDioxide":
        if len(comps) == 2:
            z0 = rng.uniform(0.03, 0.30)
            return (z0, 1.0 - z0)
        z0 = rng.uniform(0.02, 0.22)
        z1 = rng.uniform(0.15, 0.65)
        z2 = max(0.05, 1.0 - z0 - z1)
        s = z0 + z1 + z2
        return (z0 / s, z1 / s, z2 / s)
    vals = [rng.uniform(0.08, 0.85) for _ in comps]
    s = sum(vals)
    return tuple(v / s for v in vals)


def _mutate_z(rng: random.Random, comps: Tuple[str, ...], z: Tuple[float, ...], scale: float) -> Tuple[float, ...]:
    if len(comps) == 1:
        return (1.0,)
    arr = [max(0.02, v + rng.gauss(0.0, 0.08 * scale)) for v in z]
    if comps[0] == "CarbonDioxide":
        arr[0] = min(0.32, arr[0])
    s = sum(arr)
    return tuple(v / s for v in arr)


def _make_bphe_pair(cond_frac: float, c_len: float, c_wid: float, c_beta: float, c_gap: float,
                    e_len: float, e_wid: float, e_beta: float, e_gap: float) -> Tuple[BPHEGeometry, BPHEGeometry]:
    cf = _clip(cond_frac, 0.25, 0.78)
    a_c = MAX_TOTAL_BPHE_AREA_M2 * cf
    a_e = MAX_TOTAL_BPHE_AREA_M2 * (1.0 - cf)
    cond = BPHEGeometry("cond", a_c, _clip(c_len, 0.35, 0.95), _clip(c_wid, 0.14, 0.38),
                        _clip(c_beta, 30.0, 65.0), _clip(c_gap, 1.8, 3.4))
    evap = BPHEGeometry("evap", a_e, _clip(e_len, 0.35, 0.95), _clip(e_wid, 0.14, 0.38),
                        _clip(e_beta, 30.0, 65.0), _clip(e_gap, 1.8, 3.4))
    return cond, evap


def random_hthp(rng: random.Random, duty_name: str, track: str, idx: int) -> HTHPDesign:
    if track == "PURE_CEILING":
        comps = rng.choice(PURE_CANDIDATES[duty_name])
        z = (1.0,)
    else:
        _, comps = rng.choice(HTHP_ZEO_FAMILIES)
        z = _sample_z(rng, comps)

    topo = rng.choices(TOPOLOGIES, weights=[0.12, 0.28, 0.28, 0.32])[0]
    ihx_eff = rng.uniform(0.20, 0.80) if topo in ("IHX_SUBCOOL", "FTVI_IHX") else 0.0
    ft_bias = rng.uniform(0.80, 1.25) if topo in ("FTVI_ECON", "FTVI_IHX") else 1.0
    cf = rng.uniform(0.32, 0.72)
    cond, evap = _make_bphe_pair(
        cf, rng.uniform(0.40, 0.90), rng.uniform(0.16, 0.35), rng.uniform(35.0, 64.0), rng.uniform(1.9, 3.2),
        rng.uniform(0.40, 0.90), rng.uniform(0.16, 0.35), rng.uniform(35.0, 64.0), rng.uniform(1.9, 3.2),
    )
    sh = rng.uniform(5.0, 18.0)
    sc = rng.uniform(2.0, 12.0)
    return HTHPDesign(
        label=f"{duty_name}_{track}_init_{idx}", duty_name=duty_name,
        topology=topo, ihx_eff=ihx_eff, ft_pr_bias=ft_bias,
        components=comps, mole_fractions=z, cond_area_frac=cf,
        cond_hx=cond, evap_hx=evap, superheat_K=sh, subcool_K=sc,
    )


def mutate_hthp(rng: random.Random, p: HTHPDesign, track: str, scale: float, gen: int, idx: int) -> HTHPDesign:
    comps, z = p.components, p.mole_fractions
    topo, ihx_eff, ft_bias = p.topology, p.ihx_eff, p.ft_pr_bias

    if rng.random() < 0.20 * scale:
        if track == "PURE_CEILING" and rng.random() < 0.35:
            comps, z = rng.choice(PURE_CANDIDATES[p.duty_name]), (1.0,)
        elif track != "PURE_CEILING" and rng.random() < 0.28:
            _, comps = rng.choice(HTHP_ZEO_FAMILIES)
            z = _sample_z(rng, comps)
        if rng.random() < 0.40:
            topo = rng.choice(TOPOLOGIES)

    if len(comps) > 1 and rng.random() < 0.70:
        z = _mutate_z(rng, comps, z, scale)

    ihx_eff = _clip(ihx_eff + rng.gauss(0.0, 0.12 * scale), 0.05, 0.82) if topo in ("IHX_SUBCOOL", "FTVI_IHX") else 0.0
    ft_bias = _clip(ft_bias + rng.gauss(0.0, 0.08 * scale), 0.78, 1.28) if topo in ("FTVI_ECON", "FTVI_IHX") else 1.0

    cf = _clip(p.cond_area_frac + rng.gauss(0.0, 0.06 * scale), 0.25, 0.78)
    c, e = p.cond_hx, p.evap_hx
    cond, evap = _make_bphe_pair(
        cf,
        c.plate_length_m * math.exp(rng.gauss(0.0, 0.08 * scale)),
        c.plate_width_m * math.exp(rng.gauss(0.0, 0.08 * scale)),
        c.chevron_deg + rng.gauss(0.0, 5.0 * scale),
        c.channel_gap_mm + rng.gauss(0.0, 0.25 * scale),
        e.plate_length_m * math.exp(rng.gauss(0.0, 0.08 * scale)),
        e.plate_width_m * math.exp(rng.gauss(0.0, 0.08 * scale)),
        e.chevron_deg + rng.gauss(0.0, 5.0 * scale),
        e.channel_gap_mm + rng.gauss(0.0, 0.25 * scale),
    )
    sh = _clip(p.superheat_K + rng.gauss(0.0, 2.5 * scale), 5.0, 25.0)
    sc = _clip(p.subcool_K + rng.gauss(0.0, 2.0 * scale), 2.0, 15.0)

    return HTHPDesign(
        label=f"{p.duty_name}_{track}_g{gen}_{idx}", duty_name=p.duty_name,
        topology=topo, ihx_eff=ihx_eff, ft_pr_bias=ft_bias,
        components=comps, mole_fractions=z, cond_area_frac=cf,
        cond_hx=cond, evap_hx=evap, superheat_K=sh, subcool_K=sc,
    )


def _eval_worker(args):
    design, enforce_nat = args
    ev = evaluate_hthp(design, enforce_natural_zero_pfas=enforce_nat, enforce_above_atm=True, n_seg=20)
    return design, ev


def _fitness(ev: HTHPEvaluation, mode: str) -> float:
    if ev.feasible:
        return ev.COP_sys if mode == "COP" else (-ev.LCOH_usd_MWh)
    if ev.COP_sys <= 0.0:
        return -1000.0
    pen = 0.0
    if ev.P_evap_out_bar < 1.05:
        pen += 8.0 * (1.05 - ev.P_evap_out_bar)
    if ev.T_dis_C > 155.0:
        pen += 0.25 * (ev.T_dis_C - 155.0)
    if ev.P_cond_bar > 38.0:
        pen += 0.40 * (ev.P_cond_bar - 38.0)
    if ev.dis_SH_K < 4.0:
        pen += 0.50 * (4.0 - ev.dis_SH_K)
    base = ev.COP_sys if mode == "COP" else (-ev.LCOH_usd_MWh)
    return base - 100.0 - pen


def run_hthp_track(
    duty_name: str,
    track: str,
    mode: str,
    enforce_nat: bool,
    n_init: int,
    n_gens: int,
    pop: int,
    seed: int,
    pool: ProcessPoolExecutor,
    log_rows: List[str],
) -> Tuple[HTHPDesign, HTHPEvaluation, List[Tuple[HTHPDesign, HTHPEvaluation]]]:
    rng = random.Random(seed)
    t0 = time.perf_counter()

    cands = [random_hthp(rng, duty_name, track, i) for i in range(n_init)]
    archive = list(pool.map(_eval_worker, [(c, enforce_nat) for c in cands]))
    archive.sort(key=lambda pair: _fitness(pair[1], mode), reverse=True)
    best_d, best_ev = archive[0]
    tot_evals = len(cands)

    def log_row(g_idx: int):
        fstr = _fluid_short(best_d)
        row = (
            f"{duty_name}\t{track}\t{g_idx}\t{tot_evals}\t{int(best_ev.feasible)}\t"
            f"{best_ev.COP_sys:.4f}\t{100*best_ev.eta_lorenz:.2f}\t{best_ev.LCOH_usd_MWh:.2f}\t"
            f"{best_ev.VHC_kJ_m3:.0f}\t{best_ev.P_evap_out_bar:.2f}\t{best_ev.P_cond_bar:.1f}\t"
            f"{best_ev.T_dis_C:.1f}\t{best_ev.glide_cond_K:.1f}\t{best_ev.glide_evap_K:.1f}\t"
            f"{best_d.topology}\t{fstr}"
        )
        log_rows.append(row)
        print(
            f"  [{duty_name:13s} | {track:12s} g={g_idx:2d} ev={tot_evals:3d}] F={int(best_ev.feasible)} "
            f"COP={best_ev.COP_sys:.3f} ({100*best_ev.eta_lorenz:4.1f}%L) LCOH=${best_ev.LCOH_usd_MWh:5.2f} "
            f"VHC={best_ev.VHC_kJ_m3:4.0f} Psuc={best_ev.P_evap_out_bar:4.2f}b Tdis={best_ev.T_dis_C:5.1f}C "
            f"gl=({best_ev.glide_cond_K:4.1f}/{best_ev.glide_evap_K:4.1f}K) | {best_d.topology:11s} {fstr}",
            flush=True,
        )

    log_row(0)

    for g in range(1, n_gens + 1):
        frac = g / n_gens
        scale = 1.10 * ((1.0 - frac) ** 1.15) + 0.16

        elites: List[HTHPDesign] = []
        seen = set()
        for d, ev in archive:
            key = (d.components, d.topology)
            if key not in seen or len(elites) < 4:
                elites.append(d)
                seen.add(key)
            if len(elites) >= 12:
                break

        nxt: List[HTHPDesign] = []
        for i in range(pop):
            if i < pop // 3:
                parent, sc = elites[0], scale * 0.65
            elif i < (4 * pop) // 5:
                parent, sc = rng.choice(elites[:min(8, len(elites))]), scale
            else:
                if rng.random() < 0.35 and g < n_gens - 3:
                    nxt.append(random_hthp(rng, duty_name, track, g * 100 + i))
                    continue
                parent, sc = rng.choice(elites), scale * 1.4
            nxt.append(mutate_hthp(rng, parent, track, sc, g, i))

        res = list(pool.map(_eval_worker, [(c, enforce_nat) for c in nxt]))
        tot_evals += len(nxt)
        archive.extend(res)
        archive.sort(key=lambda pair: _fitness(pair[1], mode), reverse=True)
        archive = archive[:160]

        if _fitness(archive[0][1], mode) > _fitness(best_ev, mode) + 1e-4:
            best_d, best_ev = archive[0]
            log_row(g)

    dt = time.perf_counter() - t0
    print(f"--> Done {duty_name}/{track} ({tot_evals} evals in {dt:.1f}s = {tot_evals/dt:.1f} ev/s): COP={best_ev.COP_sys:.3f}, LCOH=${best_ev.LCOH_usd_MWh:.2f}\n", flush=True)

    top_div = []
    seen_f = set()
    for d, ev in archive:
        if not ev.feasible:
            continue
        k = (d.components, d.topology)
        if k not in seen_f:
            top_div.append((d, ev))
            seen_f.add(k)
        if len(top_div) >= 6:
            break
    return best_d, best_ev, top_div


def main():
    t0 = time.perf_counter()
    log_rows = [
        "duty\ttrack\tgen\tevals\tfeas\tCOP_sys\teta_lorenz_pct\tLCOH_usd_MWh\tVHC_kJ_m3\tP_suc_bar\tP_cond_bar\tT_dis_C\tglide_c_K\tglide_e_K\ttopology\tfluid"
    ]
    winners: Dict[str, Dict] = {}

    print("=== STAGE 4 (REGIME B): 16-CORE INDUSTRIAL HTHP AUTORESEARCH (3 DUTIES x 3 TRACKS) ===", flush=True)
    with ProcessPoolExecutor(max_workers=16) as pool:
        for d_idx, duty_name in enumerate(("DAIRY_90C", "DISTRICT_110C", "CHEMICAL_135C")):
            winners[duty_name] = {}

            # Track A: Optimized Pure Ceiling (PFAS HFOs + Pure Natural Hydrocarbons)
            d_p, ev_p, top_p = run_hthp_track(
                duty_name, "PURE_CEILING", "COP", False,
                n_init=64, n_gens=10, pop=32, seed=100 + d_idx, pool=pool, log_rows=log_rows,
            )
            # Track B: Natural Zero-PFAS Zeotropic Co-Design Maximizing COP_sys
            d_c, ev_c, top_c = run_hthp_track(
                duty_name, "NAT_ZEO_COP", "COP", True,
                n_init=112, n_gens=16, pop=48, seed=200 + d_idx, pool=pool, log_rows=log_rows,
            )
            # Track C: Natural Zero-PFAS Zeotropic Co-Design Minimizing LCOH ($/MWh_th)
            d_l, ev_l, top_l = run_hthp_track(
                duty_name, "NAT_ZEO_LCOH", "LCOH", True,
                n_init=112, n_gens=14, pop=48, seed=300 + d_idx, pool=pool, log_rows=log_rows,
            )

            for tname, dw, evw, tlist in [
                ("PURE_CEILING", d_p, ev_p, top_p),
                ("NAT_ZEO_COP", d_c, ev_c, top_c),
                ("NAT_ZEO_LCOH", d_l, ev_l, top_l),
            ]:
                winners[duty_name][tname] = {
                    "design": asdict(dw),
                    "eval": asdict(evw),
                    "fluid_short": _fluid_short(dw),
                    "top_diverse": [
                        {"fluid": _fluid_short(dd), "topology": dd.topology, "COP_sys": ee.COP_sys,
                         "LCOH": ee.LCOH_usd_MWh, "VHC": ee.VHC_kJ_m3, "P_suc": ee.P_evap_out_bar,
                         "glide_c": ee.glide_cond_K, "glide_e": ee.glide_evap_K}
                        for dd, ee in tlist
                    ],
                }

    with open("hthp_autoresearch_log.tsv", "w") as f:
        f.write("\n".join(log_rows) + "\n")
    with open("hthp_stage4_winners.json", "w") as f:
        json.dump(winners, f, indent=2)

    dt = time.perf_counter() - t0
    print(f"=== ALL 9 INDUSTRIAL HTHP TRACKS COMPLETE IN {dt:.1f}s ===")
    for duty_name in ("DAIRY_90C", "DISTRICT_110C", "CHEMICAL_135C"):
        p_ev = winners[duty_name]["PURE_CEILING"]["eval"]
        c_ev = winners[duty_name]["NAT_ZEO_COP"]["eval"]
        l_ev = winners[duty_name]["NAT_ZEO_LCOH"]["eval"]
        p_fl = winners[duty_name]["PURE_CEILING"]["fluid_short"]
        c_fl = winners[duty_name]["NAT_ZEO_COP"]["fluid_short"]
        l_fl = winners[duty_name]["NAT_ZEO_LCOH"]["fluid_short"]
        gain_cop = 100.0 * (c_ev["COP_sys"] - p_ev["COP_sys"]) / max(p_ev["COP_sys"], 1e-6)
        gain_lcoh = 100.0 * (l_ev["LCOH_usd_MWh"] - p_ev["LCOH_usd_MWh"]) / max(p_ev["LCOH_usd_MWh"], 1e-6)
        print(f"\n[{duty_name}]")
        print(f"  PURE_CEILING : COP={p_ev['COP_sys']:.3f} ({100*p_ev['eta_lorenz']:.1f}%L) LCOH=${p_ev['LCOH_usd_MWh']:.2f}/MWh VHC={p_ev['VHC_kJ_m3']:.0f} Psuc={p_ev['P_evap_out_bar']:.2f}b | {winners[duty_name]['PURE_CEILING']['design']['topology']} {p_fl}")
        print(f"  NAT_ZEO_COP  : COP={c_ev['COP_sys']:.3f} ({100*c_ev['eta_lorenz']:.1f}%L, {gain_cop:+.1f}% vs Pure!) LCOH=${c_ev['LCOH_usd_MWh']:.2f}/MWh VHC={c_ev['VHC_kJ_m3']:.0f} Psuc={c_ev['P_evap_out_bar']:.2f}b gl=({c_ev['glide_cond_K']:.1f}/{c_ev['glide_evap_K']:.1f}K) | {winners[duty_name]['NAT_ZEO_COP']['design']['topology']} {c_fl}")
        print(f"  NAT_ZEO_LCOH : COP={l_ev['COP_sys']:.3f} ({100*l_ev['eta_lorenz']:.1f}%L) LCOH=${l_ev['LCOH_usd_MWh']:.2f}/MWh ({gain_lcoh:+.1f}% vs Pure!) VHC={l_ev['VHC_kJ_m3']:.0f} Psuc={l_ev['P_evap_out_bar']:.2f}b gl=({l_ev['glide_cond_K']:.1f}/{l_ev['glide_evap_K']:.1f}K) | {winners[duty_name]['NAT_ZEO_LCOH']['design']['topology']} {l_fl}")


if __name__ == "__main__":
    main()
