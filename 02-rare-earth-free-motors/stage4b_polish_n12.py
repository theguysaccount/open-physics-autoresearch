"""
stage4b_polish_n12.py — Anti-Aliased Exact High-Resolution Polish (mesh_scale=1.0, n_steps=12)
Starting directly from Stage-4 winners across all 3 tracks (NDFEB_OPT, MAGNETFREE_SYNRM_OPT,
FERRITE_CO_DESIGN_OPT) to eliminate 18th-harmonic Nyquist aliasing on 6p/54s, optimize
low/mid-speed 5th/7th harmonic injection & asymmetric angles, and tighten bus-voltage & demag margins.
"""

from __future__ import annotations
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import copy
import json
import random
import time
from concurrent.futures import ProcessPoolExecutor

from stage4_motor_autoresearch import (
    _worker_eval, design_to_dict, dict_to_design, _mutate_design, _repair_design
)
from motor_evaluator import evaluate_motor_design


def polish_track(pool, winners, track_name: str, n_gens: int, pop_per_gen: int, log_f, seed_offset: int):
    rng = random.Random(202699 + seed_offset)
    t0 = time.perf_counter()
    base_d = dict_to_design(winners[track_name]["design"])
    base_eval_count = winners[track_name]["total_evals"]

    # Evaluate exact starting point at mesh=1.0, n_steps=12
    print(f"\n=== STAGE-4B ANTI-ALIASED POLISH (mesh=1.0, n_steps=12): {track_name} ({n_gens}x{pop_per_gen}={n_gens*pop_per_gen} evals) ===", flush=True)
    r0 = evaluate_motor_design(base_d, mesh_scale=1.0, n_steps=12, run_full_demag=True)
    best_feas = (r0["score"], base_d, r0) if r0["feasible"] else None
    best_any = (r0["score"], base_d, r0)
    archive = [(r0["score"], r0["feasible"], base_d, r0)]
    evals = base_eval_count

    modes = ["ripple_asym_harmonic", "full_polish", "magnet_demag_stress", "stator_bus_opt"]
    mw = [0.42, 0.30, 0.14, 0.14]

    for gen in range(n_gens):
        archive.sort(key=lambda x: x[0], reverse=True)
        elites = [x[2] for x in archive[:4]]
        scale = 0.55 * (1.0 - 0.55 * (gen / max(n_gens - 1, 1)))

        batch = []
        for c_i in range(pop_per_gen):
            parent = elites[0] if c_i < pop_per_gen // 2 else rng.choice(elites)
            mode = rng.choices(modes, weights=mw, k=1)[0]
            child = _mutate_design(parent, rng, scale, mode, track_name, f"{track_name}_P{gen}_C{c_i}")
            # On first polish gen, explicitly try higher N_turns (+0.25/+0.50) & harmonic probes on Ferrite/NdFeB!
            if gen == 0 and c_i == 0:
                child = copy.deepcopy(base_d)
                child.name = f"{track_name}_P0_TurnUp"
                child.N_turns += 0.25
                child = _repair_design(child)
            elif gen == 0 and c_i == 1:
                child = copy.deepcopy(base_d)
                child.name = f"{track_name}_P0_TurnUp2"
                child.N_turns += 0.50
                child = _repair_design(child)
            elif gen == 0 and c_i == 2 and track_name == "FERRITE_CO_DESIGN_OPT":
                child = copy.deepcopy(base_d)
                child.name = f"{track_name}_P0_DemagPad"
                child.barriers[0].thick_mm += 0.5
                child.barriers[1].thick_mm += 0.5
                child = _repair_design(child)
            batch.append(child)

        tasks = [(design_to_dict(d), 1.0, 12) for d in batch]
        results = list(pool.map(_worker_eval, tasks))

        new_best = False
        for d_dict, res in results:
            evals += 1
            d_obj = dict_to_design(d_dict)
            sc = res.get("score", -1e6)
            feas = bool(res.get("feasible", False))
            if sc > -5e5:
                archive.append((sc, feas, d_obj, res))
                if sc > best_any[0]:
                    best_any = (sc, d_obj, res)
                if feas and (best_feas is None or sc > best_feas[0]):
                    best_feas = (sc, d_obj, res)
                    new_best = True
            log_f.write(
                f"{track_name}\t{evals}\t{100+gen}\t{d_obj.name}\t{res.get('topology','?')}\t"
                f"{len(d_obj.barriers)}\t{int(feas)}\t{sc:.3f}\t"
                f"{res.get('eta_sys_wltc_pct',0):.3f}\t{res.get('cos_phi_wltc',0):.4f}\t"
                f"{res.get('ripple_wltc_pct',99):.2f}\t{res.get('T_peak_Nm',0):.1f}\t"
                f"{res.get('S_inv_kVA',0):.1f}\t{res.get('total_sys_cost_usd',0):.1f}\t"
                f"{res.get('sf_stress',0):.2f}\t{res.get('demag_margin_T',-9):.3f}\n"
            )
        log_f.flush()
        archive = sorted(archive, key=lambda x: x[0], reverse=True)[:40]

        cur = best_feas if best_feas is not None else best_any
        sc_b, d_b, r_b = cur
        tag = "FEAS" if r_b["feasible"] else "INFEAS"
        star = " *** NEW BEST!" if new_best else ""
        print(
            f"  [{track_name} P{gen:02d} | {evals:3d} tot | {time.perf_counter()-t0:5.1f}s] {tag:6s} "
            f"{r_b['topology']} {len(d_b.barriers)}B Nt={d_b.N_turns:.2f} | score={sc_b:6.2f} | "
            f"eta_sys={r_b['eta_sys_wltc_pct']:5.2f}% | cos={r_b['cos_phi_wltc']:.3f} | "
            f"Rip={r_b['ripple_wltc_pct']:4.1f}% | Tpk={r_b['T_peak_Nm']:5.1f}Nm | "
            f"Inv={r_b['S_inv_kVA']:5.1f}kVA | Cost=${r_b['total_sys_cost_usd']:6.0f} | "
            f"SF={r_b['sf_stress']:.2f} | Demag={r_b['demag_margin_T']:+.3f}T{star}",
            flush=True
        )

    final_t = best_feas if best_feas is not None else best_any
    _, win_d, win_r = final_t
    winners[track_name]["total_evals"] = evals
    winners[track_name]["design"] = design_to_dict(win_d)
    winners[track_name]["result_n12"] = win_r
    # Also record fast mesh=1.25, n=6 baseline of the polished winner for clean convergence table
    winners[track_name]["result_stage4"] = evaluate_motor_design(win_d, mesh_scale=1.25, n_steps=6, run_full_demag=True)


def main():
    t0 = time.perf_counter()
    winners = json.load(open("motor_stage4_winners.json"))
    n_workers = min(16, os.cpu_count() or 8)
    with open("motor_autoresearch_log.tsv", "a") as log_f:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            polish_track(pool, winners, "NDFEB_OPT", n_gens=5, pop_per_gen=16, log_f=log_f, seed_offset=11)
            polish_track(pool, winners, "MAGNETFREE_SYNRM_OPT", n_gens=4, pop_per_gen=16, log_f=log_f, seed_offset=22)
            polish_track(pool, winners, "FERRITE_CO_DESIGN_OPT", n_gens=8, pop_per_gen=16, log_f=log_f, seed_offset=33)

    with open("motor_stage4_winners.json", "w") as f:
        json.dump(winners, f, indent=2)
    print(f"\n=== STAGE-4B ANTI-ALIASED POLISH COMPLETE in {(time.perf_counter()-t0)/60.0:.2f} min ===", flush=True)


if __name__ == "__main__":
    main()
