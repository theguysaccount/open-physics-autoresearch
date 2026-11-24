"""
Stage 3A: Benchmark real industry baselines (R-410A, R-32, pure R-290, and
best textbook pure R-290 + FTVI + IHX) on the locked 6-bin dual-climate harness
and print exact Gouy-Stodola Watt-by-Watt exergy destruction profiles.
"""

from __future__ import annotations

import time
from typing import List

from cycle_evaluator import (
    INDOOR_H,
    INDOOR_W,
    OUTDOOR_H,
    OUTDOOR_W,
    STANDARD_BINS,
    SystemDesign,
    SystemEvaluation,
    evaluate_system,
)
from microchannel_hx import CoilGeometry


def standard_indoor_coil(counterflow_bridge: bool = False) -> CoilGeometry:
    # Realistic 7 kW indoor MCHX slab (single-row 20 mm slab, 2 passes)
    return CoilGeometry(
        name="indoor_std",
        width=INDOOR_W,
        height=INDOOR_H,
        n_rows=1,
        slab_depth=0.020,
        tube_pitch=0.0095,
        tube_height=0.0018,
        port_dh=0.0010,
        fin_pitch=0.0014,
        pass_fractions=(0.55, 0.45),
        counterflow_bridge=counterflow_bridge,
    )


def standard_outdoor_coil(counterflow_bridge: bool = False) -> CoilGeometry:
    # Realistic 7 kW outdoor MCHX slab (single-row 22 mm slab, 2 passes)
    return CoilGeometry(
        name="outdoor_std",
        width=OUTDOOR_W,
        height=OUTDOOR_H,
        n_rows=1,
        slab_depth=0.022,
        tube_pitch=0.0095,
        tube_height=0.0018,
        port_dh=0.0011,
        fin_pitch=0.0015,
        pass_fractions=(0.60, 0.40),
        counterflow_bridge=counterflow_bridge,
    )


# Standard airflows for a 7 kW split AC/HP (~0.32 m3/s indoor max, ~0.85 m3/s outdoor max)
STD_V_AIR_IN = (0.32, 0.32, 0.24, 0.22, 0.28, 0.34)
STD_V_AIR_OUT = (0.90, 0.85, 0.60, 0.55, 0.72, 0.88)
STD_SH = (5.0, 5.0, 5.0, 5.0, 5.0, 5.0)
STD_SC = (5.0, 5.0, 5.0, 5.0, 5.0, 5.0)
ZERO_HOLD = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def build_baselines() -> List[SystemDesign]:
    return [
        SystemDesign(
            label="Baseline 1: R-410A (1-stage TXV, GWP=2088)",
            topology="TXV_BASIC", ihx_eff=0.0, ft_pr_bias=1.0, comp_shift_enabled=False,
            components=("R410A",), mole_fractions=(1.0,),
            indoor_coil=standard_indoor_coil(), outdoor_coil=standard_outdoor_coil(),
            comp_disp_cc=22.0,
            v_air_in_m3s=STD_V_AIR_IN, v_air_out_m3s=STD_V_AIR_OUT,
            superheat_K=STD_SH, subcool_K=STD_SC, liq_hold_frac=ZERO_HOLD,
        ),
        SystemDesign(
            label="Baseline 2: R-32 (1-stage TXV, GWP=675)",
            topology="TXV_BASIC", ihx_eff=0.0, ft_pr_bias=1.0, comp_shift_enabled=False,
            components=("R32",), mole_fractions=(1.0,),
            indoor_coil=standard_indoor_coil(), outdoor_coil=standard_outdoor_coil(),
            comp_disp_cc=20.0,
            v_air_in_m3s=STD_V_AIR_IN, v_air_out_m3s=STD_V_AIR_OUT,
            superheat_K=STD_SH, subcool_K=STD_SC, liq_hold_frac=ZERO_HOLD,
        ),
        SystemDesign(
            label="Baseline 3: Pure R-290 (1-stage TXV, GWP=3)",
            topology="TXV_BASIC", ihx_eff=0.0, ft_pr_bias=1.0, comp_shift_enabled=False,
            components=("Propane",), mole_fractions=(1.0,),
            indoor_coil=standard_indoor_coil(), outdoor_coil=standard_outdoor_coil(),
            comp_disp_cc=32.0,
            v_air_in_m3s=STD_V_AIR_IN, v_air_out_m3s=STD_V_AIR_OUT,
            superheat_K=STD_SH, subcool_K=STD_SC, liq_hold_frac=ZERO_HOLD,
        ),
        SystemDesign(
            label="Baseline 4: Pure R-290 + FTVI + IHX (Best Textbook Pure)",
            topology="FTVI_IHX", ihx_eff=0.55, ft_pr_bias=1.0, comp_shift_enabled=False,
            components=("Propane",), mole_fractions=(1.0,),
            indoor_coil=standard_indoor_coil(), outdoor_coil=standard_outdoor_coil(),
            comp_disp_cc=32.0,
            v_air_in_m3s=STD_V_AIR_IN, v_air_out_m3s=STD_V_AIR_OUT,
            superheat_K=STD_SH, subcool_K=STD_SC, liq_hold_frac=ZERO_HOLD,
        ),
    ]


def print_summary(ev: SystemEvaluation, dt_ms: float) -> None:
    status = "PASS" if ev.feasible else f"FAIL ({'; '.join(ev.violations)})"
    print(f"\n[{ev.label}]  ({dt_ms:.0f} ms)  ->  {status}")
    print(
        f"  ASPF={ev.ASPF:.3f} | Cool SEER={ev.SEER_cool:.3f} | Heat SCOP={ev.SCOP_heat:.3f} | "
        f"46C COP={ev.COP_46C_T3:.3f} | -7C COP={ev.COP_m7C_cold:.3f}"
    )
    print(
        f"  Charge={ev.max_charge_g:.0f}g (flam={ev.flammable_charge_g:.0f}g) | "
        f"MaxTdis={ev.max_T_dis_C:.1f}C | MaxP={ev.max_P_high_bar:.1f}bar | "
        f"Al=({ev.al_mass_in_kg:.2f}+{ev.al_mass_out_kg:.2f})kg | +BOM=${ev.inc_bom_usd:.0f}"
    )
    print("  Bin breakdown:")
    for b in ev.bins:
        tot_I = b.I_comp_W + b.I_cond_W + b.I_evap_W + b.I_exp_W + b.I_ihx_ft_W + b.I_fans_W + b.W_backup_W
        print(
            f"    {b.bin_name:13s}: COP={b.COP_sys:5.2f} ({100*b.COP_sys/max(b.COP_carnot_ref,0.1):4.1f}% Carnot) | "
            f"f={b.f_comp_Hz:5.1f}Hz P=({b.P_evap_bar:4.1f}->{b.P_cond_bar:4.1f}bar) Tdis={b.T_dis_C:5.1f}C "
            f"Back={b.W_backup_W:4.0f}W R_air=({100*b.r_air_evap:2.0f}/{100*b.r_air_cond:2.0f}%) | "
            f"Exergy[W]: Comp={b.I_comp_W:4.0f} Cond={b.I_cond_W:4.0f} Evap={b.I_evap_W:4.0f} "
            f"Exp={b.I_exp_W:3.0f} Fan={b.I_fans_W:3.0f} (Tot={tot_I:4.0f}W)"
        )


if __name__ == "__main__":
    print("=== LOCKED DUAL-CLIMATE BENCHMARK: INDUSTRY BASELINES & EXERGY PROFILE ===")
    for d in build_baselines():
        t0 = time.perf_counter()
        ev = evaluate_system(d, enforce_natural_gwp=False)
        dt_ms = (time.perf_counter() - t0) * 1000.0
        print_summary(ev, dt_ms)
