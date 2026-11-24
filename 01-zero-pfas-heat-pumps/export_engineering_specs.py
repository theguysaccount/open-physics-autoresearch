"""
Export complete manufacturable engineering & prototype test specifications
for BOTH code-verified breakthroughs:
  - Regime A: 7 kW Dual-Climate Reversible Split AC / Heat Pump (Air-to-Air MCHX)
  - Regime B: 100 kW_th Zero-PFAS Industrial HTHP Skids (90 C, 110 C, 135 C BPHE)
"""

from __future__ import annotations

import json
from typing import Dict

from brazed_plate_hx import BPHEGeometry
from cycle_evaluator import STANDARD_BINS, evaluate_system
from hthp_evaluator import INDUSTRIAL_DUTIES, evaluate_hthp
from stage5_hthp_verification import fluid_pretty, hthp_from_dict
from stage5_verification_gauntlet import design_from_dict, fluid_label


def main():
    with open("stage4_winners.json") as f:
        s4_air = json.load(f)
    with open("verification_report.json") as f:
        v5_air = json.load(f)
    with open("hthp_stage4_winners.json") as f:
        s4_hthp = json.load(f)
    with open("hthp_verification_report.json") as f:
        v5_hthp = json.load(f)

    # Regime A Winning Prototype Spec
    d_air = design_from_dict(s4_air["ZEOTROPE_FULL_OPT"]["design"])
    ev_air = evaluate_system(d_air, enforce_natural_gwp=True, max_glide_K=10.0, n_seg=60)
    ic, oc = d_air.indoor_coil, d_air.outdoor_coil

    regime_a_spec = {
        "system_name": "7 kW Dual-Climate Reversible Split AC / Cold-Climate Heat Pump (Zero-PFAS Natural)",
        "refrigerant": {
            "blend": fluid_label(d_air),
            "components": list(d_air.components),
            "mole_fractions": [round(x, 4) for x in d_air.mole_fractions],
            "GWP100": round(ev_air.gwp, 2),
            "safety_class": "A3 (IEC 60335-2-40 Ed.7 compliant, <50% of 988g indoor limit)",
            "rated_charge_g": round(ev_air.flammable_charge_g, 1),
            "phase_glide_K": round(ev_air.max_glide_K, 2),
        },
        "cycle_architecture": {
            "topology": d_air.topology,
            "ihx_effectiveness": round(d_air.ihx_eff, 3),
            "compressor_displacement_cc_per_rev": round(d_air.comp_disp_cc, 2),
            "speed_range_Hz": [20.0, 110.0],
            "incremental_bom_usd": round(ev_air.inc_bom_usd, 1),
        },
        "indoor_microchannel_coil_CAD": {
            "face_W_x_H_mm": [round(ic.width * 1e3, 1), round(ic.height * 1e3, 1)],
            "n_rows": ic.n_rows,
            "slab_depth_mm": round(ic.slab_depth * 1e3, 2),
            "total_core_depth_mm": round(ic.total_depth * 1e3, 2),
            "tube_pitch_mm": round(ic.tube_pitch * 1e3, 2),
            "tube_minor_height_mm": round(ic.tube_height * 1e3, 2),
            "port_hydraulic_diameter_mm": round(ic.port_dh * 1e3, 2),
            "ports_per_tube": ic.n_ports_per_tube,
            "tubes_per_row": ic.n_tubes_per_row,
            "louver_fin_pitch_mm": round(ic.fin_pitch * 1e3, 2),
            "pass_fractions": list(ic.pass_fractions),
            "reversible_counterflow_bridge": ic.counterflow_bridge,
            "aluminum_mass_kg": round(ic.aluminum_mass, 2),
            "internal_volume_mL": round(ic.internal_volume * 1e6, 1),
        },
        "outdoor_microchannel_coil_CAD": {
            "face_W_x_H_mm": [round(oc.width * 1e3, 1), round(oc.height * 1e3, 1)],
            "n_rows": oc.n_rows,
            "slab_depth_mm": round(oc.slab_depth * 1e3, 2),
            "total_core_depth_mm": round(oc.total_depth * 1e3, 2),
            "tube_pitch_mm": round(oc.tube_pitch * 1e3, 2),
            "tube_minor_height_mm": round(oc.tube_height * 1e3, 2),
            "port_hydraulic_diameter_mm": round(oc.port_dh * 1e3, 2),
            "ports_per_tube": oc.n_ports_per_tube,
            "tubes_per_row": oc.n_tubes_per_row,
            "louver_fin_pitch_mm": round(oc.fin_pitch * 1e3, 2),
            "pass_fractions": list(oc.pass_fractions),
            "reversible_counterflow_bridge": oc.counterflow_bridge,
            "aluminum_mass_kg": round(oc.aluminum_mass, 2),
            "internal_volume_mL": round(oc.internal_volume * 1e6, 1),
        },
        "bin_control_schedule": [
            {
                "bin": op.name, "T_out_C": op.T_out_C, "T_in_C": op.T_in_C, "Q_load_W": op.Q_load_W,
                "v_air_in_m3s": round(d_air.v_air_in_m3s[i], 3),
                "v_air_out_m3s": round(d_air.v_air_out_m3s[i], 3),
                "superheat_K": round(d_air.superheat_K[i], 1),
                "subcool_K": round(d_air.subcool_K[i], 1),
                "COP_sys": round(ev_air.bins[i].COP_sys, 3),
                "P_evap_bar": round(ev_air.bins[i].P_evap_bar, 2),
                "P_cond_bar": round(ev_air.bins[i].P_cond_bar, 2),
                "T_dis_C": round(ev_air.bins[i].T_dis_C, 1),
                "f_comp_Hz": round(ev_air.bins[i].f_comp_Hz, 1),
            }
            for i, op in enumerate(STANDARD_BINS)
        ],
        "verified_performance": {
            "ASPF_n60": round(ev_air.ASPF, 4),
            "SEER_cool": round(ev_air.SEER_cool, 3),
            "SCOP_heat": round(ev_air.SCOP_heat, 3),
            "COP_46C_T3": round(ev_air.COP_46C_T3, 3),
            "COP_m7C_cold": round(ev_air.COP_m7C_cold, 3),
            "rayleigh_leak_20pct_delta_pct": round(v5_air["rayleigh_leak_20pct"]["T3_ZEOTROPE_OPT"]["delta_pct"], 2),
            "monte_carlo_mean_std": [
                round(v5_air["monte_carlo_robustness"]["T3_ZEOTROPE_OPT"]["mean"], 4),
                round(v5_air["monte_carlo_robustness"]["T3_ZEOTROPE_OPT"]["std"], 4),
            ],
        },
    }

    # Regime B Winning Industrial HTHP Specs (all 3 duties, both Max-COP & Min-LCOH)
    regime_b_specs: Dict[str, Dict] = {}
    for du in ("DAIRY_90C", "DISTRICT_110C", "CHEMICAL_135C"):
        regime_b_specs[du] = {}
        for tk in ("NAT_ZEO_COP", "NAT_ZEO_LCOH"):
            dh = hthp_from_dict(s4_hthp[du][tk]["design"])
            evh = evaluate_hthp(dh, enforce_natural_zero_pfas=True, n_seg=80)
            rep = v5_hthp[du][tk]
            regime_b_specs[du][tk] = {
                "duty": {
                    "sink_C": [INDUSTRIAL_DUTIES[du].T_sink_in_C, INDUSTRIAL_DUTIES[du].T_sink_out_C],
                    "source_C": [INDUSTRIAL_DUTIES[du].T_src_in_C, INDUSTRIAL_DUTIES[du].T_src_out_C],
                    "Q_cond_kW": 100.0,
                },
                "refrigerant": {
                    "blend": fluid_pretty(dh),
                    "components": list(dh.components),
                    "mole_fractions": [round(x, 4) for x in dh.mole_fractions],
                    "GWP100": round(evh.gwp, 2),
                    "PFAS_free": True,
                    "cond_glide_K": round(evh.glide_cond_K, 1),
                    "evap_glide_K": round(evh.glide_evap_K, 1),
                },
                "cycle_and_compressor": {
                    "topology": dh.topology,
                    "ihx_effectiveness": round(dh.ihx_eff, 3),
                    "superheat_K": round(dh.superheat_K, 1),
                    "subcool_K": round(dh.subcool_K, 1),
                    "P_evap_out_bar": round(evh.P_evap_out_bar, 2),
                    "P_cond_bar": round(evh.P_cond_bar, 2),
                    "PR": round(evh.PR, 2),
                    "T_dis_C": round(evh.T_dis_C, 1),
                    "dis_superheat_K": round(evh.dis_SH_K, 1),
                    "VHC_kJ_m3": round(evh.VHC_kJ_m3, 0),
                    "required_suction_flow_m3_h": round(evh.V_suc_m3_h, 1),
                    "W_comp_kW": round(evh.W_comp_kW, 2),
                    "W_pumps_kW": round(evh.W_pumps_kW, 2),
                },
                "condenser_BPHE_CAD": {
                    "area_m2": round(dh.cond_hx.area_m2, 2),
                    "L_x_W_mm": [round(dh.cond_hx.plate_length_m * 1e3, 0), round(dh.cond_hx.plate_width_m * 1e3, 0)],
                    "chevron_deg": round(dh.cond_hx.chevron_deg, 1),
                    "channel_gap_mm": round(dh.cond_hx.channel_gap_mm, 2),
                    "channels_per_side": dh.cond_hx.n_channels_per_side,
                    "pinch_min_K": round(evh.pinch_cond_K, 2),
                },
                "evaporator_BPHE_CAD": {
                    "area_m2": round(dh.evap_hx.area_m2, 2),
                    "L_x_W_mm": [round(dh.evap_hx.plate_length_m * 1e3, 0), round(dh.evap_hx.plate_width_m * 1e3, 0)],
                    "chevron_deg": round(dh.evap_hx.chevron_deg, 1),
                    "channel_gap_mm": round(dh.evap_hx.channel_gap_mm, 2),
                    "channels_per_side": dh.evap_hx.n_channels_per_side,
                    "pinch_min_K": round(evh.pinch_evap_K, 2),
                },
                "verified_metrics": {
                    "COP_sys_n80": round(evh.COP_sys, 3),
                    "eta_lorenz_pct": round(100.0 * evh.eta_lorenz, 1),
                    "LCOH_usd_MWh": round(evh.LCOH_usd_MWh, 2),
                    "CAPEX_usd": round(evh.CAPEX_usd, 0),
                    "dual_eos_pcsaft_dCOP_pct": rep["dual_eos_pcsaft"].get("cop_diff_pct"),
                    "rayleigh_leak_20pct": rep["rayleigh_leak_20pct"],
                    "mc_COP_mean_std": [rep["mc_COP_mean"], rep["mc_COP_std"]],
                },
            }

    payload = {
        "regime_A_air_to_air_7kW": regime_a_spec,
        "regime_B_industrial_HTHP_100kW": regime_b_specs,
    }
    with open("engineering_spec_sheets.json", "w") as f:
        json.dump(payload, f, indent=2)
    print("Exported complete manufacturable engineering specs -> engineering_spec_sheets.json!")


if __name__ == "__main__":
    main()
