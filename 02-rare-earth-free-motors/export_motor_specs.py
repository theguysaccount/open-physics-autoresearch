"""
export_motor_specs.py — Export Manufacturable CAD Geometry, Magnet Dimensions,
Winding & Inverter Specs, and 5-Point Harmonic Control Map for Motor Autoresearch Winners.
"""

from __future__ import annotations
import json
import math
from motor_geometry import MotorDesign, BarrierSpec, R_SO, R_SH, L_STK, K_FILL_MAX


def _export_one(entry: dict, verif_tr: dict, fp_tr: dict, mc_tr: dict) -> dict:
    d_dct = entry["design"]
    r = entry["result_n12"]
    th_c_deg = 90.0 / d_dct["p_pairs"]
    r_avail_min = R_SH * 1000.0 + 3.5
    r_avail_max = d_dct["R_ro_mm"] - 2.6
    span_r = r_avail_max - r_avail_min

    barriers_cad = []
    for i, b in enumerate(d_dct["barriers"]):
        u_mm = r_avail_min + b["r_frac"] * span_r
        barriers_cad.append({
            "barrier_index": i + 1,
            "radial_center_u_mm": round(u_mm, 2),
            "thickness_mm": round(b["thick_mm"], 2),
            "central_half_width_frac": round(b["half_width_frac"], 3),
            "wing_tip_left_deg": round(b["beta_L_frac"] * th_c_deg, 2),
            "wing_tip_right_deg": round(b["beta_R_frac"] * th_c_deg, 2),
            "asymmetry_deg": round((b["beta_R_frac"] - b["beta_L_frac"]) * th_c_deg, 2),
            "outer_bridge_mm": round(b["bridge_mm"], 2),
            "center_rib_mm": round(b["rib_mm"], 2),
            "mag_center_fill": round(b["mag_center_fill"], 3),
            "mag_wing_fill": round(b["mag_wing_fill"], 3),
        })

    return {
        "track": entry["track"],
        "design_name": d_dct["name"],
        "topology": r["topology"],
        "magnet_grade": d_dct["magnet_grade"],
        "rare_earth_free": d_dct["magnet_grade"] != "NDFEB_N42SH",
        "envelope_and_stator_CAD": {
            "stator_OD_mm": round(2.0 * R_SO * 1000.0, 1),
            "stack_length_mm": round(L_STK * 1000.0, 1),
            "shaft_OD_mm": round(2.0 * R_SH * 1000.0, 1),
            "rotor_OD_mm": round(2.0 * d_dct["R_ro_mm"], 2),
            "mechanical_air_gap_mm": round(d_dct["gap_mm"], 3),
            "tooth_width_ratio": round(d_dct["k_tooth"], 3),
            "stator_yoke_mm": round(d_dct["h_yoke_mm"], 2),
            "slot_opening_mm": round(d_dct["b_so_mm"], 2),
            "tooth_tip_height_mm": round(d_dct["h_tt_mm"], 2),
            "effective_series_turns_per_slot": d_dct["N_turns"],
            "slot_fill_factor": K_FILL_MAX,
            "R_ph_120C_mOhm": round(r["R_ph_mOhm"], 2),
            "L_end_winding_uH": round(r["L_ew_uH"], 1),
        },
        "rotor_barriers_CAD": barriers_cad,
        "rotor_surface_notch": {
            "depth_mm": round(d_dct["notch_depth_mm"], 2),
            "angle_from_pole_center_deg": round(d_dct["notch_angle_frac"] * th_c_deg, 2),
            "angular_width_deg": round(d_dct["notch_width_deg"], 2),
        },
        "harmonic_current_injection": {
            "k_h5_frac": round(d_dct["k_h5"], 4),
            "phi_h5_deg": round(math.degrees(d_dct["phi_h5_rad"]), 1),
            "k_h7_frac": round(d_dct["k_h7"], 4),
            "phi_h7_deg": round(math.degrees(d_dct["phi_h7_rad"]), 1),
        },
        "verified_performance_fine_n24": verif_tr["fine_m065_n24"],
        "frozen_permeability_OP2": fp_tr,
        "monte_carlo_tolerance_32draws": mc_tr,
        "bom_and_inverter_economics": {
            "m_steel_net_kg": round(r["bom"]["m_steel_net_kg"], 2),
            "m_cu_kg": round(r["bom"]["m_cu_kg"], 2),
            "m_mag_kg": round(r["bom"]["m_mag_kg"], 2),
            "bom_motor_usd": round(r["bom"]["bom_motor_usd"], 1),
            "sic_inverter_kVA": round(r["S_inv_kVA"], 1),
            "sic_inverter_cost_usd": round(r["inv_cost_usd"], 1),
            "total_drive_system_cost_usd": round(r["total_sys_cost_usd"], 1),
        },
        "wltc_operating_schedule": r["ops"],
        "peak_launch_2500rpm": r["op5_peak"],
    }


def main():
    winners = json.load(open("motor_stage4_winners.json"))
    verif = json.load(open("motor_verification_report.json"))
    specs = {}
    for tr in ("NDFEB_OPT", "MAGNETFREE_SYNRM_OPT", "FERRITE_CO_DESIGN_OPT"):
        specs[tr] = _export_one(
            winners[tr],
            verif["mesh_step_convergence"][tr],
            verif["frozen_permeability"][tr],
            verif["monte_carlo_tolerance"][tr],
        )
    with open("motor_engineering_specs.json", "w") as f:
        json.dump(specs, f, indent=2)
    print("Exported motor_engineering_specs.json cleanly!")


if __name__ == "__main__":
    main()
