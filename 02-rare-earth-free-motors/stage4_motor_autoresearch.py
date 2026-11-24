"""
stage4_motor_autoresearch.py — 16-Core Multi-Layer Evolutionary Autoresearch
for Rare-Earth-Free Electric Motor + SiC Inverter-kVA (cos phi) & Harmonic Control Co-Design.

Runs 3 Head-to-Head Tracks under identical locked physics, 2D centrifugal stress (SF >= 1.50 @ 14.4krpm),
0.5 mm cut-edge lamination degradation, AC proximity copper loss, 3D end-winding R/L,
and worst-case demagnetization (-30C cold soak for Ferrite; +140C hot for NdFeB):
  - Track 1 (NDFEB_OPT):           Rare-earth N42SH NdFeB optimized reference ceiling
  - Track 2 (MAGNETFREE_SYNRM_OPT): Zero-magnet pure reluctance SynRM + asymmetric barriers + harmonic control
  - Track 3 (FERRITE_CO_DESIGN_OPT): Rare-earth-free FB12B ceramic ferrite + asymmetric flux concentration
                                     + pole/slot topology + inverter-kVA (cos phi) + 5th/7th harmonic co-design
"""

from __future__ import annotations
import os
# Prevent BLAS/OpenMP oversubscription inside 16-worker ProcessPoolExecutor
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import copy
import json
import math
import random
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from typing import List, Tuple, Dict

import numpy as np

from motor_geometry import (
    MotorDesign, BarrierSpec, VALID_TOPOLOGIES,
    validate_design_geometry, preset_ndfeb_vipm,
    preset_pure_synrm_4p, preset_ferrite_pmasynrm_6p, preset_ferrite_pmasynrm_4p,
)
from motor_evaluator import evaluate_motor_design


def _worker_eval(args):
    """Isolated process worker (each worker has its own clean Gmsh instance)."""
    design_dict, mesh_scale, n_steps = args
    barriers = [BarrierSpec(**b) for b in design_dict["barriers"]]
    d = MotorDesign(**{k: v for k, v in design_dict.items() if k != "barriers"}, barriers=barriers)
    res = evaluate_motor_design(d, mesh_scale=mesh_scale, n_steps=n_steps, run_full_demag=True)
    return design_dict, res


def design_to_dict(d: MotorDesign) -> dict:
    out = asdict(d)
    return out


def dict_to_design(dct: dict) -> MotorDesign:
    barriers = [BarrierSpec(**b) for b in dct["barriers"]]
    return MotorDesign(**{k: v for k, v in dct.items() if k != "barriers"}, barriers=barriers)


def _repair_design(d: MotorDesign) -> MotorDesign:
    """Clamp & enforce strict ordering, non-overlapping barriers, and >=0.80 mm structural minimums."""
    d.R_ro_mm = float(np.clip(d.R_ro_mm, 53.0, 69.5))
    d.gap_mm = float(np.clip(d.gap_mm, 0.65, 1.15))
    d.k_tooth = float(np.clip(d.k_tooth, 0.37, 0.62))
    # Pole-aware realistic yoke bounds (4p needs thicker yoke, 8p thinner)
    y_min = 18.0 if d.p_pairs == 2 else (13.0 if d.p_pairs == 3 else 10.5)
    y_max = min(28.5, 100.0 - d.R_ro_mm - d.gap_mm - 10.5)
    if y_max < y_min:
        d.R_ro_mm = min(d.R_ro_mm, 100.0 - d.gap_mm - 10.8 - y_min)
        y_max = 100.0 - d.R_ro_mm - d.gap_mm - 10.5
    d.h_yoke_mm = float(np.clip(d.h_yoke_mm, y_min, y_max))

    r_si = d.R_ro_mm + d.gap_mm
    slot_pitch = 2.0 * math.pi * r_si / d.q_slots
    tooth_w = slot_pitch * d.k_tooth
    b_so_max = max(1.8, (slot_pitch - tooth_w) * 0.88)
    d.b_so_mm = float(np.clip(d.b_so_mm, 1.7, min(5.0, b_so_max)))
    d.h_tt_mm = float(np.clip(d.h_tt_mm, 0.8, 2.0))
    # Snap N_turns to realistic quarter-turns (integer turns with 1, 2, or 4 parallel paths)
    d.N_turns = round(float(np.clip(d.N_turns, 1.75, 6.0)) * 4.0) / 4.0

    d.notch_depth_mm = float(np.clip(d.notch_depth_mm, 0.0, 1.4))
    if d.notch_depth_mm < 0.15:
        d.notch_depth_mm = 0.0
    d.notch_angle_frac = float(np.clip(d.notch_angle_frac, 0.55, 0.92))
    d.notch_width_deg = float(np.clip(d.notch_width_deg, 0.8, 3.2))

    d.k_h5 = float(np.clip(d.k_h5, 0.0, 0.10))
    if d.k_h5 < 0.004:
        d.k_h5 = 0.0
    d.phi_h5_rad = float((d.phi_h5_rad + math.pi) % (2.0 * math.pi) - math.pi)
    d.k_h7 = float(np.clip(d.k_h7, 0.0, 0.07))
    if d.k_h7 < 0.004:
        d.k_h7 = 0.0
    d.phi_h7_rad = float((d.phi_h7_rad + math.pi) % (2.0 * math.pi) - math.pi)

    # Sort barriers inner -> outer and enforce strict radial steel channels >= 1.65 mm & angular separation
    n_b = len(d.barriers)
    d.barriers.sort(key=lambda b: b.r_frac)
    r_avail_mm = (d.R_ro_mm - 2.6) - (26.0 + 3.5)

    for idx, b in enumerate(d.barriers):
        b.thick_mm = float(np.clip(b.thick_mm, 2.0, 8.8 if d.magnet_grade == "FERRITE_FB12B" else 6.8))
        b.half_width_frac = float(np.clip(b.half_width_frac, 0.08, 0.68))
        b.bridge_mm = float(np.clip(b.bridge_mm, 0.80, 2.80))
        if b.rib_mm < 0.40:
            b.rib_mm = 0.0
        else:
            b.rib_mm = float(np.clip(b.rib_mm, 0.80, 3.20))
        if d.magnet_grade == "NONE":
            b.mag_center_fill = 0.0
            b.mag_wing_fill = 0.0
        else:
            b.mag_center_fill = 0.0 if b.mag_center_fill < 0.12 else float(np.clip(b.mag_center_fill, 0.15, 0.92))
            b.mag_wing_fill = 0.0 if b.mag_wing_fill < 0.12 else float(np.clip(b.mag_wing_fill, 0.15, 0.86))

    # Enforce monotonic r_frac spacing & non-overlap
    if n_b == 1:
        d.barriers[0].r_frac = float(np.clip(d.barriers[0].r_frac, 0.20, 0.70))
        d.barriers[0].beta_L_frac = float(np.clip(d.barriers[0].beta_L_frac, 0.40, 0.91))
        d.barriers[0].beta_R_frac = float(np.clip(d.barriers[0].beta_R_frac, 0.40, 0.91))
    else:
        # Shrink thicknesses slightly if sum exceeds 64% of available rotor radial span
        tot_t = sum(b.thick_mm for b in d.barriers)
        max_t = max(8.0, r_avail_mm - (n_b - 1) * 2.0 - 4.0)
        if tot_t > max_t:
            sc = max_t / tot_t
            for b in d.barriers:
                b.thick_mm = max(2.0, b.thick_mm * sc)

        d.barriers[0].r_frac = float(np.clip(d.barriers[0].r_frac, 0.10, 0.36))
        for i in range(1, n_b):
            min_dr = (0.5 * (d.barriers[i - 1].thick_mm + d.barriers[i].thick_mm) + 1.80) / max(r_avail_mm, 15.0)
            lo = d.barriers[i - 1].r_frac + min_dr
            hi = min(0.91, 1.0 - (n_b - 1 - i) * 0.17)
            d.barriers[i].r_frac = float(np.clip(max(d.barriers[i].r_frac, lo), lo, max(lo, hi)))

        # Enforce strictly decreasing beta_L and beta_R from inner -> outer with >= 0.13 separation
        d.barriers[0].beta_L_frac = float(np.clip(d.barriers[0].beta_L_frac, 0.62, 0.92))
        d.barriers[0].beta_R_frac = float(np.clip(d.barriers[0].beta_R_frac, 0.62, 0.92))
        for i in range(1, n_b):
            max_L = d.barriers[i - 1].beta_L_frac - 0.13
            max_R = d.barriers[i - 1].beta_R_frac - 0.13
            min_b = max(0.16, 0.16 + (n_b - 1 - i) * 0.13)
            d.barriers[i].beta_L_frac = float(np.clip(d.barriers[i].beta_L_frac, min_b, max(min_b, max_L)))
            d.barriers[i].beta_R_frac = float(np.clip(d.barriers[i].beta_R_frac, min_b, max(min_b, max_R)))

    return d


def _seed_track(track_name: str) -> List[MotorDesign]:
    """Build diverse, physics-informed initial seeds spanning topologies & barrier architectures."""
    seeds: List[MotorDesign] = []

    if track_name == "NDFEB_OPT":
        # 1. Standard 8p/48s V-IPM (tuned N_turns)
        s0 = preset_ndfeb_vipm()
        s0.name = "S0_NdFeB_8p48s"
        seeds.append(s0)
        # 2. 6p/36s V-IPM (lower iron loss at 11.5krpm!)
        s1 = copy.deepcopy(s0)
        s1.name = "S1_NdFeB_6p36s"
        s1.p_pairs, s1.q_spp = 3, 2
        s1.R_ro_mm, s1.h_yoke_mm, s1.N_turns = 62.5, 16.5, 4.25
        seeds.append(s1)
        # 3. 6p/54s Double-V NdFeB (smooth 9-slot/pole ripple!)
        s2 = copy.deepcopy(s1)
        s2.name = "S2_NdFeB_6p54s_2V"
        s2.q_spp, s2.N_turns = 3, 3.0
        s2.barriers[1].mag_wing_fill = 0.55
        s2.barriers[1].rib_mm = 1.1
        seeds.append(s2)
        # 4. 4p/48s V-IPM (ultra-low high-speed iron loss + high saliency!)
        s3 = copy.deepcopy(s0)
        s3.name = "S3_NdFeB_4p48s"
        s3.p_pairs, s3.q_spp = 2, 4
        s3.R_ro_mm, s3.h_yoke_mm, s3.N_turns = 60.0, 21.0, 3.25
        seeds.append(s3)
        # 5. 8p/48s Asymmetric V-IPM with notches
        s4 = copy.deepcopy(s0)
        s4.name = "S4_NdFeB_8p48s_Asym"
        s4.notch_depth_mm = 0.6
        s4.barriers[0].beta_L_frac = 0.80
        s4.barriers[0].beta_R_frac = 0.72
        s4.barriers[1].beta_L_frac = 0.44
        s4.barriers[1].beta_R_frac = 0.35
        seeds.append(s4)

    elif track_name == "MAGNETFREE_SYNRM_OPT":
        # 1. 4p/36s 3-barrier strong-rib SynRM
        s0 = preset_pure_synrm_4p()
        s0.name = "S0_SynRM_4p36s_3B"
        s0.N_turns = 3.25
        s0.barriers[0].rib_mm, s0.barriers[1].rib_mm, s0.barriers[2].rib_mm = 1.9, 1.5, 1.0
        seeds.append(s0)
        # 2. 4p/48s 4-barrier high-saliency SynRM
        s1 = MotorDesign(
            name="S1_SynRM_4p48s_4B", magnet_grade="NONE", p_pairs=2, q_spp=4,
            R_ro_mm=61.0, gap_mm=0.65, k_tooth=0.47, h_yoke_mm=20.0, b_so_mm=2.2, h_tt_mm=1.0, N_turns=2.5,
            barriers=[
                BarrierSpec(0.14, 4.4, 0.44, 0.89, 0.89, 1.05, 1.9, 0.0, 0.0),
                BarrierSpec(0.38, 4.0, 0.38, 0.71, 0.71, 1.00, 1.5, 0.0, 0.0),
                BarrierSpec(0.62, 3.4, 0.32, 0.51, 0.51, 0.95, 1.1, 0.0, 0.0),
                BarrierSpec(0.84, 2.6, 0.24, 0.29, 0.29, 0.90, 0.0, 0.0, 0.0),
            ],
        )
        seeds.append(s1)
        # 3. 4p/36s 4-barrier SynRM
        s2 = copy.deepcopy(s1)
        s2.name = "S2_SynRM_4p36s_4B"
        s2.q_spp, s2.N_turns = 3, 3.25
        seeds.append(s2)
        # 4. 6p/54s 3-barrier SynRM (shorter end-windings, larger R_ro!)
        s3 = MotorDesign(
            name="S3_SynRM_6p54s_3B", magnet_grade="NONE", p_pairs=3, q_spp=3,
            R_ro_mm=64.0, gap_mm=0.65, k_tooth=0.47, h_yoke_mm=15.5, b_so_mm=2.1, h_tt_mm=1.0, N_turns=2.75,
            barriers=[
                BarrierSpec(0.18, 5.2, 0.45, 0.88, 0.88, 1.10, 1.8, 0.0, 0.0),
                BarrierSpec(0.50, 4.2, 0.36, 0.62, 0.62, 1.00, 1.3, 0.0, 0.0),
                BarrierSpec(0.80, 3.0, 0.26, 0.34, 0.34, 0.95, 0.0, 0.0, 0.0),
            ],
        )
        seeds.append(s3)
        # 5. 6p/36s 3-barrier SynRM
        s4 = copy.deepcopy(s3)
        s4.name = "S4_SynRM_6p36s_3B"
        s4.q_spp, s4.N_turns = 2, 4.0
        seeds.append(s4)

    else:  # FERRITE_CO_DESIGN_OPT
        # Diverse seeds explicitly exploring:
        # (a) Demag-shielded air-outer-barrier + thick inner/mid flux-concentrating ferrite!
        # (b) 4p/36s, 4p/48s, 6p/36s, 6p/54s, 8p/48s topologies!
        # (c) 2-, 3-, and 4-barrier configurations with SF >= 1.50 structural ribs!

        # Seed 0: 4p/48s 3-Barrier Demag-Shielded (2 thick inner ferrite + 1 air outer shield!)
        seeds.append(MotorDesign(
            name="S0_Fer_4p48s_3B_Shield", magnet_grade="FERRITE_FB12B", p_pairs=2, q_spp=4,
            R_ro_mm=61.0, gap_mm=0.65, k_tooth=0.47, h_yoke_mm=20.0, b_so_mm=2.2, h_tt_mm=1.0, N_turns=3.25,
            barriers=[
                BarrierSpec(0.18, 6.4, 0.42, 0.88, 0.88, 1.25, 2.2, 0.90, 0.82),
                BarrierSpec(0.52, 5.6, 0.36, 0.62, 0.62, 1.15, 1.8, 0.90, 0.78),
                BarrierSpec(0.82, 3.0, 0.26, 0.34, 0.34, 1.00, 1.0, 0.00, 0.00),  # Air demag shield!
            ],
        ))
        # Seed 1: 4p/36s 3-Barrier Demag-Shielded (thick ferrite)
        s1 = copy.deepcopy(seeds[0])
        s1.name = "S1_Fer_4p36s_3B_Shield"
        s1.q_spp, s1.N_turns = 3, 4.25
        seeds.append(s1)
        # Seed 2: 6p/54s 3-Barrier Demag-Shielded (compact yoke + bigger R_ro = 64.5 mm!)
        seeds.append(MotorDesign(
            name="S2_Fer_6p54s_3B_Shield", magnet_grade="FERRITE_FB12B", p_pairs=3, q_spp=3,
            R_ro_mm=64.0, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=15.5, b_so_mm=2.1, h_tt_mm=1.0, N_turns=3.25,
            barriers=[
                BarrierSpec(0.18, 6.8, 0.44, 0.87, 0.87, 1.30, 2.2, 0.90, 0.84),
                BarrierSpec(0.52, 5.8, 0.36, 0.60, 0.60, 1.20, 1.7, 0.90, 0.78),
                BarrierSpec(0.82, 3.0, 0.26, 0.33, 0.33, 1.05, 1.0, 0.00, 0.00),  # Air demag shield!
            ],
        ))
        # Seed 3: 6p/36s 3-Barrier Demag-Shielded
        s3 = copy.deepcopy(seeds[2])
        s3.name = "S3_Fer_6p36s_3B_Shield"
        s3.q_spp, s3.N_turns = 2, 4.75
        seeds.append(s3)
        # Seed 4: 4p/48s 4-Barrier (3 thick inner ferrite + 1 air outer shield!)
        seeds.append(MotorDesign(
            name="S4_Fer_4p48s_4B_Shield", magnet_grade="FERRITE_FB12B", p_pairs=2, q_spp=4,
            R_ro_mm=61.5, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=19.5, b_so_mm=2.2, h_tt_mm=1.0, N_turns=3.25,
            barriers=[
                BarrierSpec(0.13, 5.6, 0.42, 0.89, 0.89, 1.25, 2.3, 0.90, 0.82),
                BarrierSpec(0.39, 5.2, 0.38, 0.71, 0.71, 1.15, 1.9, 0.90, 0.78),
                BarrierSpec(0.64, 4.6, 0.32, 0.51, 0.51, 1.10, 1.4, 0.88, 0.70),
                BarrierSpec(0.86, 2.6, 0.24, 0.28, 0.28, 1.00, 0.9, 0.00, 0.00),
            ],
        ))
        # Seed 5: 6p/54s 2-Barrier Ultra-Thick Deep-V Ferrite Flux-Concentrator (zero outer demag!)
        seeds.append(MotorDesign(
            name="S5_Fer_6p54s_2B_ThickV", magnet_grade="FERRITE_FB12B", p_pairs=3, q_spp=3,
            R_ro_mm=64.5, gap_mm=0.65, k_tooth=0.48, h_yoke_mm=15.5, b_so_mm=2.1, h_tt_mm=1.0, N_turns=3.25,
            barriers=[
                BarrierSpec(0.24, 7.8, 0.28, 0.86, 0.86, 1.35, 2.4, 0.90, 0.85),
                BarrierSpec(0.68, 6.2, 0.32, 0.50, 0.50, 1.20, 1.8, 0.90, 0.78),
            ],
        ))
        # Seed 6: 4p/36s 2-Barrier Ultra-Thick Deep-V Ferrite
        s6 = copy.deepcopy(seeds[5])
        s6.name = "S6_Fer_4p36s_2B_ThickV"
        s6.p_pairs, s6.q_spp = 2, 3
        s6.R_ro_mm, s6.h_yoke_mm, s6.N_turns = 61.0, 20.5, 4.25
        seeds.append(s6)
        # Seed 7: 8p/48s 2-Barrier Thick Ferrite (tests 8-pole directly!)
        seeds.append(MotorDesign(
            name="S7_Fer_8p48s_2B_Thick", magnet_grade="FERRITE_FB12B", p_pairs=4, q_spp=2,
            R_ro_mm=65.5, gap_mm=0.65, k_tooth=0.49, h_yoke_mm=13.0, b_so_mm=2.2, h_tt_mm=1.0, N_turns=3.75,
            barriers=[
                BarrierSpec(0.25, 7.5, 0.32, 0.84, 0.84, 1.30, 2.2, 0.90, 0.84),
                BarrierSpec(0.69, 5.5, 0.30, 0.48, 0.48, 1.15, 1.5, 0.88, 0.72),
            ],
        ))

    return [_repair_design(s) for s in seeds]


def _mutate_design(parent: MotorDesign, rng: random.Random, scale: float, mode: str,
                   track_name: str, child_name: str) -> MotorDesign:
    """Physics-targeted mutation operator across all 4 co-design layers."""
    d = copy.deepcopy(parent)
    d.name = child_name

    if mode == "topology_jump":
        # Jump pole-slot topology while scaling N_turns & h_yoke physically to preserve V_bus!
        old_p, old_q = d.p_pairs, d.q_spp
        choices = VALID_TOPOLOGIES if track_name != "MAGNETFREE_SYNRM_OPT" else [(2, 3), (2, 4), (3, 2), (3, 3)]
        new_p, new_q = rng.choice(choices)
        if (new_p, new_q) != (old_p, old_q):
            # V_ph ~ p_pairs * q_spp * N_turns * flux_per_pole (~ 1/p_pairs) -> V ~ q_spp * N_turns!
            d.N_turns = d.N_turns * (old_q / new_q) * ((old_p / new_p) ** 0.15)
            d.p_pairs, d.q_spp = new_p, new_q
            if new_p == 2:
                d.h_yoke_mm = rng.uniform(19.0, 23.0)
                d.R_ro_mm = min(d.R_ro_mm, 62.5)
            elif new_p == 3:
                d.h_yoke_mm = rng.uniform(14.5, 18.0)
                d.R_ro_mm = float(np.clip(d.R_ro_mm, 60.5, 65.5))
            else:
                d.h_yoke_mm = rng.uniform(11.5, 14.5)
                d.R_ro_mm = max(d.R_ro_mm, 63.0)
        # Optionally change barrier count (2 <-> 3 <-> 4)
        if rng.random() < 0.45:
            min_b, max_b = (2, 3) if track_name == "NDFEB_OPT" else ((3, 4) if track_name == "MAGNETFREE_SYNRM_OPT" else (2, 4))
            n_cur = len(d.barriers)
            if n_cur < max_b and (n_cur == min_b or rng.random() < 0.5):
                # Add an outer air-shield or inner barrier
                outer = copy.deepcopy(d.barriers[-1])
                outer.r_frac = min(0.88, outer.r_frac + 0.18)
                outer.thick_mm = 2.8
                outer.beta_L_frac = max(0.22, outer.beta_L_frac - 0.18)
                outer.beta_R_frac = max(0.22, outer.beta_R_frac - 0.18)
                outer.mag_center_fill = 0.0
                outer.mag_wing_fill = 0.0
                d.barriers.append(outer)
            elif n_cur > min_b:
                d.barriers.pop(-1)

    elif mode == "ripple_asym_harmonic":
        # Dedicated NVH / torque-ripple crusher: asymmetric wing angles + notches + 5th/7th harmonic injection!
        for b in d.barriers:
            # Independent left/right wing angle perturbation (breaks symmetry cleanly!)
            b.beta_L_frac += rng.gauss(0.0, 0.055 * scale)
            b.beta_R_frac += rng.gauss(0.0, 0.055 * scale)
        d.b_so_mm += rng.gauss(0.0, 0.35 * scale)
        if rng.random() < 0.55:
            d.notch_depth_mm = max(0.0, d.notch_depth_mm + rng.gauss(0.0, 0.35 * scale))
            d.notch_angle_frac += rng.gauss(0.0, 0.06 * scale)
            d.notch_width_deg += rng.gauss(0.0, 0.45 * scale)
        # Co-evolve 5th and 7th harmonic current amplitude & phase!
        if rng.random() < 0.75:
            d.k_h5 = max(0.0, d.k_h5 + rng.gauss(0.0, 0.025 * scale))
            d.phi_h5_rad += rng.gauss(0.0, 1.0 * scale)
            d.k_h7 = max(0.0, d.k_h7 + rng.gauss(0.0, 0.018 * scale))
            d.phi_h7_rad += rng.gauss(0.0, 1.0 * scale)

    elif mode == "magnet_demag_stress":
        # Structural ribs/bridges + barrier thickness + magnet fill distribution (solves demag & SF!)
        for idx, b in enumerate(d.barriers):
            b.thick_mm += rng.gauss(0.0, 0.85 * scale)
            b.bridge_mm += rng.gauss(0.0, 0.22 * scale)
            if rng.random() < 0.6:
                b.rib_mm = max(0.0, b.rib_mm + rng.gauss(0.0, 0.35 * scale))
            if d.magnet_grade != "NONE":
                if rng.random() < 0.30 and idx == len(d.barriers) - 1:
                    # Toggle outer barrier air-shield vs magnet
                    b.mag_center_fill = 0.0 if b.mag_center_fill > 0.1 else 0.80
                    b.mag_wing_fill = 0.0
                else:
                    b.mag_center_fill += rng.gauss(0.0, 0.18 * scale)
                    b.mag_wing_fill += rng.gauss(0.0, 0.18 * scale)
        d.N_turns += rng.choice([-0.5, -0.25, 0.0, 0.25, 0.5])

    elif mode == "stator_bus_opt":
        # Stator split ratio, tooth/yoke iron-vs-copper optimization & 400V bus turns matching!
        d.R_ro_mm += rng.gauss(0.0, 1.6 * scale)
        d.gap_mm += rng.gauss(0.0, 0.08 * scale)
        d.k_tooth += rng.gauss(0.0, 0.035 * scale)
        d.h_yoke_mm += rng.gauss(0.0, 1.5 * scale)
        d.b_so_mm += rng.gauss(0.0, 0.30 * scale)
        d.h_tt_mm += rng.gauss(0.0, 0.18 * scale)
        d.N_turns += rng.gauss(0.0, 0.55 * scale)

    else:  # "full_polish"
        d.R_ro_mm += rng.gauss(0.0, 0.7 * scale)
        d.k_tooth += rng.gauss(0.0, 0.02 * scale)
        d.h_yoke_mm += rng.gauss(0.0, 0.7 * scale)
        d.N_turns += rng.gauss(0.0, 0.25 * scale)
        for b in d.barriers:
            b.r_frac += rng.gauss(0.0, 0.035 * scale)
            b.thick_mm += rng.gauss(0.0, 0.45 * scale)
            b.half_width_frac += rng.gauss(0.0, 0.05 * scale)
            b.beta_L_frac += rng.gauss(0.0, 0.035 * scale)
            b.beta_R_frac += rng.gauss(0.0, 0.035 * scale)
            b.bridge_mm += rng.gauss(0.0, 0.12 * scale)
            b.rib_mm = max(0.0, b.rib_mm + rng.gauss(0.0, 0.18 * scale))
            if d.magnet_grade != "NONE":
                b.mag_center_fill += rng.gauss(0.0, 0.08 * scale)
                b.mag_wing_fill += rng.gauss(0.0, 0.08 * scale)
        if d.k_h5 > 0.0 or rng.random() < 0.4:
            d.k_h5 = max(0.0, d.k_h5 + rng.gauss(0.0, 0.012 * scale))
            d.phi_h5_rad += rng.gauss(0.0, 0.45 * scale)
            d.k_h7 = max(0.0, d.k_h7 + rng.gauss(0.0, 0.008 * scale))
            d.phi_h7_rad += rng.gauss(0.0, 0.45 * scale)

    return _repair_design(d)


def run_track(
    pool: ProcessPoolExecutor,
    track_name: str,
    n_gens: int,
    pop_per_gen: int,
    log_f,
    seed_offset: int,
) -> dict:
    rng = random.Random(202611 + seed_offset)
    t_track0 = time.perf_counter()
    print(f"\n==================== TRACK: {track_name} ({n_gens} gens x {pop_per_gen} = {n_gens*pop_per_gen} FEA evals) ====================", flush=True)

    # Evaluate diverse seeds first (pad to pop_per_gen with mutations)
    seeds = _seed_track(track_name)
    init_pop = list(seeds)
    while len(init_pop) < pop_per_gen:
        base = rng.choice(seeds)
        init_pop.append(_mutate_design(base, rng, 1.0, rng.choice(["magnet_demag_stress", "stator_bus_opt", "ripple_asym_harmonic"]),
                                       track_name, f"{track_name}_G0_C{len(init_pop)}"))

    archive: List[Tuple[float, bool, MotorDesign, dict]] = []
    best_feas: Tuple[float, MotorDesign, dict] | None = None
    best_any: Tuple[float, MotorDesign, dict] | None = None
    total_evals = 0

    modes = ["ripple_asym_harmonic", "magnet_demag_stress", "stator_bus_opt", "full_polish", "topology_jump"]
    mode_weights_early = [0.24, 0.24, 0.20, 0.18, 0.14]
    mode_weights_late = [0.34, 0.18, 0.16, 0.30, 0.02]

    for gen in range(n_gens):
        if gen == 0:
            batch = init_pop
        else:
            # Select top-5 diverse elites (preferring feasible designs + topology diversity)
            archive.sort(key=lambda x: x[0], reverse=True)
            elites: List[MotorDesign] = []
            seen_topos: Dict[Tuple[int, int, int], int] = {}
            for sc, feas, d_el, r_el in archive:
                key = (d_el.p_pairs, d_el.q_spp, len(d_el.barriers))
                if seen_topos.get(key, 0) < 2:
                    elites.append(d_el)
                    seen_topos[key] = seen_topos.get(key, 0) + 1
                if len(elites) >= 5:
                    break
            if not elites:
                elites = [archive[0][2]]

            prog = gen / max(n_gens - 1, 1)
            scale = 1.05 * (1.0 - 0.72 * prog)  # anneal 1.05 -> 0.29
            mw = mode_weights_early if prog < 0.55 else mode_weights_late

            batch = []
            for c_i in range(pop_per_gen):
                # Focus 60% of children on top-2 elites, 40% on diverse elites
                parent = elites[0] if c_i < pop_per_gen // 3 else (
                    elites[min(1, len(elites) - 1)] if c_i < (2 * pop_per_gen) // 3 else rng.choice(elites)
                )
                mode = rng.choices(modes, weights=mw, k=1)[0]
                child = _mutate_design(parent, rng, scale, mode, track_name, f"{track_name}_G{gen}_C{c_i}")
                batch.append(child)

        # Parallel exact nonlinear FEA evaluation across 16 cores!
        tasks = [(design_to_dict(d), 1.25, 6) for d in batch]
        results = list(pool.map(_worker_eval, tasks))

        gen_new_best = False
        for d_dict, res in results:
            total_evals += 1
            d_obj = dict_to_design(d_dict)
            sc = res.get("score", -1e6)
            feas = bool(res.get("feasible", False))
            if sc > -5e5:
                archive.append((sc, feas, d_obj, res))
                if best_any is None or sc > best_any[0]:
                    best_any = (sc, d_obj, res)
                if feas and (best_feas is None or sc > best_feas[0]):
                    best_feas = (sc, d_obj, res)
                    gen_new_best = True

            # Append TSV row
            log_f.write(
                f"{track_name}\t{total_evals}\t{gen}\t{d_obj.name}\t{res.get('topology','?')}\t"
                f"{len(d_obj.barriers)}\t{int(feas)}\t{sc:.3f}\t"
                f"{res.get('eta_sys_wltc_pct',0):.3f}\t{res.get('cos_phi_wltc',0):.4f}\t"
                f"{res.get('ripple_wltc_pct',99):.2f}\t{res.get('T_peak_Nm',0):.1f}\t"
                f"{res.get('S_inv_kVA',0):.1f}\t{res.get('total_sys_cost_usd',0):.1f}\t"
                f"{res.get('sf_stress',0):.2f}\t{res.get('demag_margin_T',-9):.3f}\n"
            )
        log_f.flush()

        # Trim archive to top 80
        archive.sort(key=lambda x: x[0], reverse=True)
        archive = archive[:80]

        cur = best_feas if best_feas is not None else best_any
        if cur is not None:
            sc_b, d_b, r_b = cur
            tag = "FEAS" if r_b["feasible"] else "INFEAS"
            star = " *** NEW BEST!" if gen_new_best else ""
            elap = time.perf_counter() - t_track0
            print(
                f"  [{track_name} G{gen:02d} | {total_evals:3d} evals | {elap:5.1f}s] {tag:6s} "
                f"{r_b['topology']} {len(d_b.barriers)}B | score={sc_b:6.2f} | "
                f"eta_sys={r_b['eta_sys_wltc_pct']:5.2f}% | cos={r_b['cos_phi_wltc']:.3f} | "
                f"Rip={r_b['ripple_wltc_pct']:4.1f}% | Tpk={r_b['T_peak_Nm']:5.1f}Nm | "
                f"Inv={r_b['S_inv_kVA']:5.1f}kVA | Cost=${r_b['total_sys_cost_usd']:6.0f} | "
                f"SF={r_b['sf_stress']:.2f} | Demag={r_b['demag_margin_T']:+.2f}T{star}",
                flush=True
            )

    # Final winner (re-evaluated at standard mesh_scale=1.0, n_steps=12 for exact reporting!)
    winner_tuple = best_feas if best_feas is not None else best_any
    assert winner_tuple is not None
    _, win_design, win_res_fast = winner_tuple
    win_res_ref = evaluate_motor_design(win_design, mesh_scale=1.0, n_steps=12, run_full_demag=True)
    print(
        f"\n>>> TRACK {track_name} WINNER (@ mesh=1.0, n_steps=12): {win_design.name} ({win_res_ref['topology']} {len(win_design.barriers)}B)\n"
        f"    feasible={win_res_ref['feasible']} ({win_res_ref['violations']}) | "
        f"eta_sys={win_res_ref['eta_sys_wltc_pct']:.2f}% (mot={win_res_ref['eta_motor_wltc_pct']:.2f}%) | "
        f"cos={win_res_ref['cos_phi_wltc']:.3f} | Rip={win_res_ref['ripple_wltc_pct']:.2f}% | "
        f"Tpk={win_res_ref['T_peak_Nm']:.1f}Nm | Inv={win_res_ref['S_inv_kVA']:.1f}kVA | "
        f"Cost=${win_res_ref['total_sys_cost_usd']:.0f} (Mag=${win_res_ref['bom']['bom_mag_usd']:.0f}, Inv=${win_res_ref['inv_cost_usd']:.0f}) | "
        f"SF={win_res_ref['sf_stress']:.2f} | Demag={win_res_ref['demag_margin_T']:+.3f}T",
        flush=True
    )
    return {
        "track": track_name,
        "total_evals": total_evals,
        "design": design_to_dict(win_design),
        "result_stage4": win_res_fast,
        "result_n12": win_res_ref,
    }


def main():
    t0 = time.perf_counter()
    n_workers = min(16, os.cpu_count() or 8)
    print(f"=== STAGE-4 MOTOR + INVERTER-kVA AUTORESEARCH ({n_workers} CPU workers) ===", flush=True)

    winners = {}
    with open("motor_autoresearch_log.tsv", "w") as log_f:
        log_f.write(
            "track\teval_id\tgen\tname\ttopology\tn_bar\tfeasible\tscore\t"
            "eta_sys_pct\tcos_phi\tripple_pct\tT_peak_Nm\tS_inv_kVA\tsys_cost_usd\tsf_stress\tdemag_margin_T\n"
        )
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            # Track 1: Rare-Earth NdFeB Reference Ceiling (14 x 16 = 224 evals)
            winners["NDFEB_OPT"] = run_track(pool, "NDFEB_OPT", n_gens=14, pop_per_gen=16, log_f=log_f, seed_offset=101)
            # Track 2: Magnet-Free Pure SynRM Ceiling (18 x 16 = 288 evals)
            winners["MAGNETFREE_SYNRM_OPT"] = run_track(pool, "MAGNETFREE_SYNRM_OPT", n_gens=18, pop_per_gen=16, log_f=log_f, seed_offset=202)
            # Track 3: Rare-Earth-Free Ferrite Co-Design Breakthrough (32 x 16 = 512 evals)
            winners["FERRITE_CO_DESIGN_OPT"] = run_track(pool, "FERRITE_CO_DESIGN_OPT", n_gens=32, pop_per_gen=16, log_f=log_f, seed_offset=303)

    with open("motor_stage4_winners.json", "w") as f:
        json.dump(winners, f, indent=2)

    dt_min = (time.perf_counter() - t0) / 60.0
    tot = sum(w["total_evals"] for w in winners.values())
    print(f"\n=== ALL 3 TRACKS COMPLETE: {tot} exact nonlinear FEA motor evaluations in {dt_min:.2f} min ===", flush=True)


if __name__ == "__main__":
    main()
