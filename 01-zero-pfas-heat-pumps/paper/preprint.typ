#set page(
  paper: "us-letter",
  margin: (x: 0.68in, y: 0.72in),
  numbering: "1",
  footer: context [
    #set text(size: 8pt, fill: rgb("#5a6778"))
    Autonomous Multi-Layer Co-Design of Zero-PFAS Heat Pumps #h(1fr) Page #counter(page).display() of #counter(page).final().first()
  ]
)
#set text(font: ("New Computer Modern", "Libertinus Serif", "Georgia", "Times New Roman"), size: 9.2pt, fill: rgb("#111827"))
#set par(justify: true, leading: 0.56em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => {
  v(0.75em, weak: true)
  text(size: 10.8pt, weight: "bold", fill: rgb("#123460"), it)
  v(0.28em, weak: true)
}
#show heading.where(level: 2): it => {
  v(0.55em, weak: true)
  text(size: 9.6pt, weight: "bold", fill: rgb("#1e3a5f"), it)
  v(0.18em, weak: true)
}

#let navy = rgb("#123460")
#let accent = rgb("#e8eef6")
#let slate = rgb("#4b5563")

// Title Banner
#align(center)[
  #text(size: 15.5pt, weight: "bold", fill: navy)[
    Autonomous Multi-Layer Co-Design of Zero-PFAS Heat Pumps Across Residential Cooling and Industrial Process Heat (\<150 °C)
  ]
  #v(0.45em)
  #text(size: 9.8pt, weight: "semibold")[Jack Jay] \
  #text(size: 8.3pt, fill: slate)[
    Independent Researcher (`jackjay.io`) · Assisted by Autonomous Multi-Agent Physics Verification Pipeline \
    Open-Source Code, CAD Specs & Dual-EOS Ledgers: #link("https://github.com/theguysaccount/open-physics-autoresearch")[`github.com/theguysaccount/open-physics-autoresearch`] · Preprint v1.0 (Nov 2026)
  ]
]

#v(0.5em)
#rect(fill: accent, stroke: (left: 2.5pt + navy), inset: 9pt, radius: 2pt, width: 100%)[
  #text(weight: "bold", fill: navy)[Abstract — ]
  Vapor-compression cooling, space heating, and low-temperature (\<200 °C) industrial process heat account for ≈20,000 TWh/yr of global thermal demand and 7–8% of anthropogenic greenhouse-gas emissions, yet commercial heat pumps achieve only 30–37% of ideal Lorenz efficiency while remaining tied to high-GWP HFCs (`R-410A`, `R-32`) or persistent PFAS-class HFOs (`R-1233zd(E)`). Applying a four-stage autonomous research pipeline—132-year first-principles literature synthesis, four-way adversarial screening, profiling micro-experiments, and a 16-core multi-layer evolutionary autoresearch loop—we jointly co-optimize discrete cycle topology, natural working-fluid composition ($"GWP"_(100) <= 3$, zero PFAS), finite-volume heat-exchanger CAD, and operating control policies under strict safety, material, and techno-economic locks. Across *8,526 exact simulations*, autoresearch uncovers a sharp *two-regime physical phase boundary*:
  #v(0.25em)
  *1. Regime A (7 kW Dual-Climate Reversible Split AC/HP, 6-bin ISO 5151 T1/T3 46 °C + EN 14825 −7 °C, ≤500 g A3 charge):* High-glide blends and two-stage flash-tank vapor injection (`FTVI`) lose under Silver–Bell–Ghaly mass-transfer resistance, suction pressure drop, wet-coil dehumidification, and the 500 g charge cap. Instead, an asymmetric 3-row open-fin single-pass indoor coil and 2-row tight-port (0.70 mm) outdoor coil with passive reversible counterflow bridge headers (+\$24 BOM) and `IHX_SUBCOOL` (+\$18 BOM) charged with `R-1270 14.8% / R-290 85.2%` (GWP = 2.86, 479 g) reaches *ASPF = 5.168* (*+40.3%* over `R-410A`, *+32.1%* over `R-32`, *+24.9%* over textbook `R-290 FTVI+IHX`), with bridge headers alone delivering *+15.4%* heating SCOP and *+18.4%* −7 °C COP, and a 20% vapor leak shifting ASPF by only *−0.24%*. \
  *2. Regime B (100 kW#sub[th] Industrial Waste-Heat Recovery at 90 / 110 / 135 °C, fixed 18.0 m² AISI 316L BPHE, $P_"suc" >= 1.05$ bar ATEX guard):* Under counterflow sensible glides (55–75 K sink, 25–35 K source), single-stage natural zeotropes shatter autoresearch-optimized two-stage PFAS `R-1233zd(E)` ceilings. Factory-sealed hermetic skids using `CO2 / R-600` (19–23% / 77–81%) exploit *isenthalpic flash-gas glide skipping* to reach *COP = 4.731 (+32.6%)* at 110 °C and *4.420 (+30.3%)* at 135 °C (46.4–46.7% Lorenz efficiency), slashing LCOH to *\$32.02–\$32.49/MWh#sub[th]* (*−28%* vs. PFAS `R-1233zd(E)`). Field-serviceable hydrocarbon zeotropes (`R-290/R-601a 48/52` at 90 °C; `R-290/R-601 55/45` at 110 °C; `R-290/R-601 72/28` at 135 °C) achieve *COP = 5.327 (+26.1%)*, *4.076 (+14.2%)*, and *4.140 (+22.1%)*, with *\<1%* agreement against SINTEF `ThermoPack PC-SAFT`, *100%* win rates across 64 Monte Carlo draws, and *−2.8% to +0.6%* COP shift after a 20% Rayleigh vapor leak while staying strictly above atmospheric pressure (2.25–7.16 bar).
]
#v(0.55em)

#show: rest => columns(2, gutter: 0.24in, rest)

= Introduction and Historical Motivation

Classical vapor-compression systems have remained anchored for nearly a century to *isothermal pure-fluid cycles* popularized after Thomas Midgley's 1930 chlorofluorocarbon synthesis, stranding four foundational thermodynamic concepts in separate silos: Hans Lorenz's (1894) sliding-temperature cycle [1], Lorenz & Meutzner's (1975) internal-HX fractionation cycles [2], Tondeur & Kvaalen's (1987) Equipartition of Entropy Production [5], and Gustav Lorentzen's (1990) natural-refrigerant revival [3].

With Kigali HFC phasedowns and global PFAS restrictions converging, both residential reversible split heat pumps (bound by IEC 60335-2-40 flammable charge caps [15]) and $<150 degree "C"$ industrial waste-heat recovery heat pumps (trapped between low-capacity PFAS HFOs and oversized pure hydrocarbons) require simultaneous structural co-design. Prior computational studies either optimize scalar operating points inside a fixed schematic or omit Silver–Bell–Ghaly (1973) zeotropic mass-transfer resistance [4], two-phase pressure drop [7], cubic fan/pump parasitics, Rouhani–Axelsson void-fraction charge [8], ATEX suction limits, compressor volumetric CAPEX, and vapor-leak fractionation. Here we enforce every penalty inside locked multiphysics evaluators and apply a 16-core evolutionary autoresearch loop across four coupled layers—topology, fluid blend, finite-volume exchanger CAD, and control policy—to map the true physical phase boundary.

= Locked Multi-Physics Models and Governing Equations

== Dual-EOS Thermodynamics & Rayleigh Leak
Primary evaluations use `CoolProp 8.1.0` multiparameter Helmholtz mixture equations of state (`HEOS`) [13]. Independent verification re-evaluates converged operating states using SINTEF `ThermoPack 2.2.3` `PC-SAFT` [14]. Field vapor-leak fractionation is computed by integrating the exact Rayleigh differential equation over 40 steps at $T_"leak" = 25 degree "C"$ until 20% of molar inventory escapes as vapor:
$ (d x_i) / (d ln N_"tot") = y_i (T_"leak", bold(x)) - x_i $

== Regime A: 7 kW Air-to-Air Reversible Split AC/HP
Six standardized climate bins combine ISO 5151 T1/T3 cooling and EN 14825 heating [16]: T3 heatwave ($46 -> 29 degree "C"$, 6.50 kW, $w=0.12$), T1 rated ($35 -> 27 degree "C"$, 7.00 kW, $w=0.23$), part-load cooling ($29 -> 27 degree "C"$, 3.85 kW, $w=0.25$), mild heating ($+12 -> 20 degree "C"$, 2.80 kW, $w=0.14$), rated heating ($+7 -> 20 degree "C"$, 4.55 kW, $w=0.16$), and cold heating ($-7 -> 20 degree "C"$, 6.20 kW, $w=0.10$, with exact $"COP"=1.0$ electric resistance backup for any heating capacity deficit):
$ "ASPF" = (sum_(b=1)^6 w_b dot(Q)_("load",b)) / (sum_(b=1)^6 w_b (dot(W)_("comp",b) + dot(W)_("fan,in",b) + dot(W)_("fan,out",b) + dot(Q)_("backup",b))) $

Indoor ($780 times 300$ mm, $<= 44$ mm depth, $<= 3.20$ kg Al) and outdoor ($850 times 620$ mm, $<= 52$ mm depth, $<= 7.50$ kg Al) microchannel coils (MCHX) are discretized into $N_"seg"$ finite-volume slabs ($N_"seg"=10$ evolution; $24, 60$ verification) using Kim & Bullard (2002) air-side $j$/$f$ correlations [6], Braun/Threlkeld wet-coil dehumidification, cubic fan power ($eta_"fan" = 0.32$), Shah (1979) condensation [9], Cooper/Chen boiling, Müller-Steinhagen & Heck (1986) two-phase $Delta P$ [7], and mandatory Silver–Bell–Ghaly (1973) mixture mass-transfer degradation [4]:
$ 1 / alpha_"mix" = 1 / alpha_"id" + Z / alpha_v, quad Z = x c_(p,v) ((partial T) / (partial h))_(P,"2ph") $
Charge is integrated with Rouhani–Axelsson (1970) void fraction [8] under strict caps: $M_"chg" <= 500$ g (A3), $T_"dis" <= 115 degree "C"$, $P_"high" <= 42$ bar, and $Delta"BOM" <=$ \$120.

== Regime B: 100 kW#sub[th] Industrial HTHP Skid
Three pressurized-water waste-heat recovery duties are locked at $100 " kW"_"th"$ and total AISI 316L brazed-plate area $A_"cond" + A_"evap" <= 18.0 " m"^2$:
- `DAIRY_90C`: Sink $35 -> 90 degree "C"$ (55 K), Source $45 -> 20 degree "C"$ (25 K).
- `DISTRICT_110C`: Sink $40 -> 110 degree "C"$ (70 K), Source $55 -> 25 degree "C"$ (30 K).
- `CHEMICAL_135C`: Sink $60 -> 135 degree "C"$ (75 K), Source $75 -> 40 degree "C"$ (35 K).

Each BPHE uses exact equal-enthalpy counterflow integration ($N_"seg"=20$ evolution; $40, 80$ verification):
$ A_"req" = sum_(k=1)^(N_"seg") (Delta dot(Q)_k) / (U_k Delta T_("lm",k)), quad 1 / U_k = 1 / alpha_("ref,mix",k) + t_w / k_"316L" + R_f + 1 / alpha_(w,k) $
with Muley & Manglik (1999) water-side HTC/$Delta P$ and pumping power [10], Han–Lee–Kim (2003) boiling [11], Yan–Lio–Lin (1999) condensation [12], Silver–Bell–Ghaly degradation, and $R_f = 5 times 10^(-5) " m"^2 "K/W"$. Safety locks enforce $P_"evap,out" >= 1.05$ bar (ATEX air-ingress guard), $P_"cond" <= 38$ bar, $T_"dis" <= 155 degree "C"$, and dry compression ($"SH"_"dis" >= 4$ K). Levelized Cost of Heat (LCOH) explicitly prices compressor swept volume:
$ "CAPEX"_"comp" = "USD 9,500" + "USD 285" dot(V)_("suc,m"^3"/h")^0.82 $
amortized over 15 yr at 7% WACC, 6,000 h/yr, 2.5% O&M, and \$85/MWh#sub[e] electricity. Exact post-IHX Gouy–Stodola exergy destruction $dot(I)_j = T_0 dot(S)_("gen",j)$ audits every component.

= Regime A: 7 kW Split AC/HP Breakthrough

#place(top, scope: "parent", float: true)[
  #figure(
    image("figures/lorenz_autoresearch_proof.png", width: 100%),
    caption: [*Regime A (7 kW Dual-Climate Reversible Split AC/HP) 4-Stage Autoresearch Proof.* (A) Verified $N_"seg"=60$ seasonal performance (ASPF, 46 °C T3, −7 °C). (B) Post-IHX Gouy–Stodola exergy audit (−35.9% vs. `R-410A`, 974 → 624 W). (C) Evolutionary ASPF trajectory across 2,382 exact HEOS simulations. (D) Stage-3 Dual-EOS (`CoolProp HEOS` vs. SINTEF `ThermoPack PC-SAFT`) glide-matching proof showing Lorenz +35% COP.]
  )
]

Stage-3 profiling micro-experiments (78 targeted trials) falsified two common assumptions: wide-glide blends (`R-290/R-600a 85/15`) degraded ASPF due to wet-coil latent flattening and mass-transfer resistance, while `FTVI_IHX` lost to `IHX_SUBCOOL` under the 500 g A3 charge cap because the flash vessel consumes liquid inventory. Running 2,382 exact `HEOS` evaluations across 16 cores produced the $N_"seg"=60$ winners in Table 1.

#place(bottom, scope: "parent", float: true)[
  #set text(size: 7.7pt)
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
      inset: (x: 4pt, y: 3.2pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else if y == 7 { rgb("#eef6ee") } else { none },
      align: (left, left, left, right, right, right, right, right, right, right),
      [*ID*], [*Fluid (GWP#sub[100])*], [*Topology*], [*ASPF*], [*SEER*], [*SCOP*], [*46 °C*], [*−7 °C*], [*Chg (g)*], [*Exergy*],
      [`B1`], [`R-410A` (2088)], [`TXV_BASIC`], [3.683], [3.560], [3.888], [2.62], [3.12], [833 (A1)], [974 W],
      [`B2`], [`R-32` (675)], [`TXV_BASIC`], [3.914], [3.827], [4.054], [2.88], [3.25], [772 (A2L)], [892 W],
      [`B3`], [`R-290` (3.0)], [`TXV_BASIC`], [3.953], [3.876], [4.076], [2.93], [3.18], [374 (A3)], [896 W],
      [`B4`], [`R-290` (3.0)], [`FTVI_IHX`], [4.137], [4.071], [4.243], [3.20], [3.45], [413 (A3)], [842 W],
      [`T1`], [`R-32` (675) AutoOpt], [`FTVI_IHX`], [5.109], [5.181], [5.001], [4.04], [3.95], [1120 (A2L)], [625 W],
      [`T2`], [`R-1270` (2.0) AutoOpt], [`IHX_SUBCOOL`], [4.983], [5.004], [4.952], [3.83], [3.88], [498 (A3)], [651 W],
      [*`T3`*], [*`R-1270/290` 15/85 (2.86)*], [*`IHX_SUBCOOL`*], [*5.168*], [*5.253*], [*5.043*], [*3.97*], [*3.92*], [*479 (A3)*], [*625 W*],
    ),
    caption: [Regime A (7 kW Dual-Climate Split AC/HP) Verified Performance at $N_"seg"=60$.]
  )
]

== Same-Hardware Fairness Ablation
Testing alternative fluids and disabling the counterflow headers on the *identical* `T3` evolved coils ($N_"seg"=24$, Table 2) isolates three genuine physical drivers:
1. *Asymmetric MCHX Geometry:* A 3-row open-fin single-pass indoor coil ($F_p=1.93$ mm, $D_h=0.97$ mm, 3.17 kg Al) cuts fan power and suction $Delta P$, while a 2-row tight-port outdoor coil ($D_h=0.70$ mm, 2 passes, 7.26 kg Al) collapses condenser liquid charge to keep total inventory at 478.8 g.
2. *Passive Counterflow Bridge Headers (+\$24 BOM):* Restoring row-to-row counterflow upon 4-way valve reversal delivers *+15.4% heating SCOP* ($4.351 -> 5.023$) and *+18.4% −7 °C COP* ($3.31 -> 3.92$) on zero added aluminum.
3. *Charge-Density-Capped Doping:* Pure `R-1270` slightly edges `R-290` thermally (5.185 vs. 5.150) but breaches the 500 g cap (513 g); 14.8% `R-1270` in `R-290` captures peak efficiency inside 479 g with only −0.24% shift after a 20% vapor leak.

#figure(
  {
    set text(size: 7.4pt)
    table(
      columns: (1fr, auto, auto, auto, auto),
      inset: (x: 3.5pt, y: 2.8pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else { none },
      align: (left, right, right, right, right),
      [*Variant on Identical `T3` Coils*], [*ASPF*], [*SCOP*], [*−7 °C*], [*Chg*],
      [*`T3` (`R-1270/290`, Bridge ON)*], [*5.156*], [*5.023*], [*3.92*], [*479 g*],
      [Pure `R-290` (Bridge ON)], [5.150], [5.016], [3.91], [472 g],
      [Pure `R-1270` (Bridge ON)], [5.185], [5.049], [3.94], [513 g\*],
      [`R-290/R-600a` 85/15 (6 K)], [5.101], [4.960], [3.85], [505 g\*],
      [`T3` Fluid, *Bridge OFF*], [4.851], [4.351], [3.31], [479 g],
    )
  },
  caption: [Same-hardware fairness ablation ($N_"seg"=24$). \*Violates 500 g A3 limit.]
)

= Regime B: Industrial \<150 °C HTHP Breakthrough

#place(top, scope: "parent", float: true)[
  #figure(
    image("figures/hthp_autoresearch_proof.png", width: 100%),
    caption: [*Regime B (100 kW#sub[th] Industrial HTHP) 4-Stage Autoresearch Proof.* (A) Converged $N_"seg"=80$ system COP beating Pure/PFAS ceilings by +26.1% to +32.6%. (B) True industrial LCOH slashing heat costs to \$32.0–\$34.3/MWh#sub[th] (31–36% below a \$50/MWh gas boiler). (C) Gouy–Stodola exergy destruction audit eliminating 33.8–37.1% of irreversibilities. (D) Dual-EOS cross-check (`CoolProp HEOS` vs. SINTEF `ThermoPack PC-SAFT`) across all six industrial zeotropic winners.]
  )
]

Across 6,144 exact equal-enthalpy BPHE + `HEOS` evaluations spanning Dairy (90 °C), District (110 °C), and Chemical (135 °C), single-stage natural zeotropes decisively beat autoresearch-optimized two-stage PFAS `R-1233zd(E)` ceilings on every thermodynamic and economic metric (Table 3).

#place(bottom, scope: "parent", float: true)[
  #set text(size: 7.3pt)
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
      inset: (x: 3.3pt, y: 2.8pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else if y in (2, 4, 5, 7, 8) { rgb("#eef6ee") } else { none },
      align: (left, left, left, right, right, right, right, right, right, right),
      [*Duty*], [*Class*], [*Fluid & Topology*], [*COP*], [*$eta_L$*], [*LCOH*], [*VHC*], [*$P_"suc"$*], [*PC-SAFT*], [*20% Leak*],
      [Dairy 90 °C], [Pure Opt], [`R-290` (`IHX`)], [4.224], [37.1%], [\$33.97], [5,295], [8.33], [−0.04%], [Pure fluid],
      [*Dairy 90 °C*], [*Best/Safe*], [*`R-290/601a 48/52` (`1S`)*], [*5.327*], [*46.8%*], [*\$32.02*], [2,399], [*2.72*], [*−0.40%*], [*−1.4% (2.25)*],
      [Dist 110 °C], [Pure PFAS], [`R-1233zd(E)` (`FTVI`)], [3.569], [35.0%], [\$44.59], [1,314], [1.28], [N/A], [Pure fluid],
      [*Dist 110 °C*], [*Hermetic*], [*`CO2/R-600 19/81` (`1S`)*], [*4.728*], [*46.4%*], [*\$32.02*], [*4,333*], [*5.13*], [−2.46%], [Hermetic only],
      [*Dist 110 °C*], [*Leak-Safe*], [*`R-290/601 55/45` (`1S`)*], [*4.076*], [*40.0%*], [*\$36.91*], [2,399], [*3.11*], [*−0.02%*], [*−2.8% (2.45)*],
      [Chem 135 °C], [Pure PFAS], [`R-1233zd(E)` (`FTVI`)], [3.391], [35.8%], [\$43.30], [2,012], [2.14], [N/A], [Pure fluid],
      [*Chem 135 °C*], [*Hermetic*], [*`CO2/R-600 23/77` (`IHX`)*], [*4.420*], [*46.7%*], [*\$32.49*], [*6,604*], [*8.70*], [*−0.88%*], [−17.0% (5.57)],
      [*Chem 135 °C*], [*Leak-Safe*], [*`R-290/601 72/28` (`IHX`)*], [*4.140*], [*43.8%*], [*\$34.33*], [*5,399*], [*8.24*], [*+0.81%*], [*+0.6% (7.16)*],
    ),
    caption: [Regime B (100 kW#sub[th] Industrial HTHP, Fixed 18.0 m² BPHE) Dual-EOS Verified Results ($N_"seg"=80$). LCOH in \$/MWh#sub[th]; VHC in kJ/m³; $P_"suc"$ in bar.]
  )
]

== Two Discovered Physical Mechanisms

*1. Isenthalpic Flash-Gas Glide Skipping:* Industrial sink/source glides exhibit ≈2:1 asymmetry (55–75 K vs. 25–35 K). Rather than requiring a costly flash-tank compositional separator, single-stage throttling drops high-glide zeotropes into the evaporator at vapor quality $x_"in" approx 0.30$–$0.40$ (producing the sharp condenser/evaporator exergy drop in Fig. 2C), automatically skipping the bottom third of the boiling dome so active evaporator glide shrinks ≈40% to match the waste-heat source while preserving the full condenser glide—achieving 46.4–46.8% Lorenz efficiency in a simple 1-stage skid.

*2. Hermetic vs. Field-Serviceable Leak-Immunity Law:* Our 20% Rayleigh leak audit reveals a sharp operational split. Ultra-high-volatility `CO2/R-600` achieves maximum COP (4.42–4.73, +30–33%), 3.3× smaller compressor displacement, and lowest LCOH (\$32.0–\$32.5/MWh) for *factory-welded hermetic skids*, whereas matched-volatility hydrocarbon zeotropes (`R-290/R-601a`, `R-290/R-601`) trade only 6–12% COP while surviving catastrophic 20% vapor leaks with *−2.8% to +0.6% COP change* and strictly above-atmospheric ATEX suction ($2.25$–$7.16$ bar).

= Five-Test Independent Verification Gauntlet

Every winner passed five adversarial audits (Table 4): (1) finite-volume mesh convergence ($N_"seg": 10 -> 60$ and $20 -> 80$, max drift $0.44\%$); (2) dual-EOS verification against SINTEF `ThermoPack PC-SAFT` ($<1\%$ error for all hydrocarbon blends and `CO2/R-600` at 135 °C); (3) 40-step Rayleigh 20% vapor-leak simulation; (4) 64-draw correlated Monte Carlo perturbations (HTC $plus.minus 25\%$, SBG $0.75$–$1.50 times$, $Delta P$ $0.80$–$1.30 times$, $eta_"comp" plus.minus 2$ pt, fouling $0.5$–$2.0 times$) yielding *100.0% win rates*; and (5) exact post-IHX Gouy–Stodola exergy closure confirming 34–37% total irreversibility reduction.

#figure(
  {
    set text(size: 7.4pt)
    table(
      columns: (1fr, auto, auto, auto, auto),
      inset: (x: 3.3pt, y: 2.8pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else { none },
      align: (left, right, right, right, right),
      [*Winning Design*], [*Mesh*], [*PC-SAFT*], [*MC Mean $plus.minus sigma$*], [*Win*],
      [Reg. A `T3` `R-1270/290`], [−0.08%], [N/A\*], [$5.165 plus.minus 0.037$], [100%],
      [Dairy `R-290/601a`], [−0.28%], [−0.40%], [$5.325 plus.minus 0.043$], [100%],
      [Dist `CO2/R-600 19/81`], [−0.44%], [−2.46%], [$4.711 plus.minus 0.037$], [100%],
      [Chem `CO2/R-600 23/77`], [+0.01%], [−0.88%], [$4.416 plus.minus 0.034$], [100%],
      [Chem `R-290/601 72/28`], [−0.30%], [+0.81%], [$4.134 plus.minus 0.016$], [100%],
    )
  },
  caption: [Five-test verification summary. \*Pure R-290 agrees within 0.14% in PC-SAFT across all 6 bins.]
)

= Limitations & Hardware Prototype Roadmap

Three second-order effects require bench validation: two-phase header maldistribution at high glide, compressor POE/PAG oil solubility shifts affecting local boiling HTC, and −7 °C transient defrost cycling. Both regimes are directly buildable from Appendix A: (1) a *\<\$120k Regime A calorimeter build* (4–6 mo) brazing the 3R-open / 2R-tight MCHX slabs with 4-check-valve headers; and (2) a *\<\$180k Regime B 100 kW#sub[th] industrial skid* (6–9 mo) pairing commercial AISI 316L BPHEs with an ATEX semi-hermetic compressor running `R-290/R-601a (48/52)` and `R-290/R-601 (72/28)`.

#heading(numbering: none)[Appendix A: Manufacturable CAD Specs]

#figure(
  {
    set text(size: 7.3pt)
    table(
      columns: (1fr, auto, auto),
      inset: (x: 3.5pt, y: 2.6pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else { none },
      align: (left, right, right),
      [*Regime A `T3` MCHX Parameter*], [*Indoor*], [*Outdoor*],
      [Face $W times H$ (mm)], [$780 times 300$], [$850 times 620$],
      [Rows $times$ Slab Depth (mm)], [$3 times 11.65$], [$2 times 14.34$],
      [Total Core Depth (mm)], [39.96], [31.18],
      [Tube Pitch $P_t$ / Height $H_t$ (mm)], [10.16 / 1.73], [8.11 / 1.62],
      [Port $D_h$ (mm) · Ports $times$ Tubes], [0.97 · $10 times 30$], [0.70 · $17 times 76$],
      [Louver Fin Pitch $F_p$ (mm)], [1.93], [1.44],
      [Passes · Counterflow Bridge], [`[1.0]` · *Yes*], [`[0.5, 0.5]` · *Yes*],
      [Al Mass / Internal Volume], [3.17 kg / 975 mL], [7.26 kg / 1,664 mL],
    )
  },
  caption: [Regime A `T3` coil CAD (`R-1270/R-290 14.8/85.2`, 478.8 g, `IHX_SUBCOOL` $epsilon=0.738$, compressor 27.16 cc/rev, +\$59.7 BOM).]
)

#figure(
  {
    set text(size: 7.1pt)
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto),
      inset: (x: 2.8pt, y: 2.6pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else { none },
      align: (left, right, right, right, right, right, right),
      [*Bin*], [*$dot(Q)$*], [*SH/SC*], [*$P_(e"/"c)$*], [*$T_"dis"$*], [*$f$*], [*COP*],
      [`COOL_46C`], [6.50], [5.0/2.5], [8.1/19.3], [88.1], [56.7], [3.971],
      [`COOL_35C`], [7.00], [5.1/6.5], [7.5/15.4], [68.5], [57.8], [5.168],
      [`COOL_29C`], [3.85], [3.6/5.3], [8.5/13.1], [51.1], [26.1], [8.110],
      [`HEAT_12C`], [2.80], [4.4/3.2], [6.2/10.7], [48.8], [23.1], [7.281],
      [`HEAT_7C`], [4.55], [5.6/2.2], [5.1/10.8], [57.4], [45.7], [5.799],
      [`HEAT_-7C`], [6.20], [4.4/11.2], [3.3/11.8], [73.8], [88.4], [3.923],
    )
  },
  caption: [Regime A `T3` six-bin operating control schedule ($dot(Q)$ in kW, $P$ in bar, $T$ in °C, $f$ in Hz).]
)

#place(bottom, scope: "parent", float: true)[
  #set text(size: 7.2pt)
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto),
      inset: (x: 3.2pt, y: 2.7pt),
      stroke: (x: none, y: 0.4pt + luma(180)),
      fill: (_, y) => if y == 0 { accent } else { none },
      align: (left, left, right, right, right, right, right, right, right),
      [*Industrial Skid*], [*Fluid & Topology*], [*Cond $A$*], [*Cond $L times W$*], [*Cond $beta$/gap/ch*], [*Evap $A$*], [*Evap $L times W$*], [*Evap $beta$/gap/ch*], [*$dot(V)_"suc"$ / CAPEX*],
      [Dairy 90 °C Best], [`R-290/601a 48/52` (`1S`)], [9.56 m²], [665×249], [65°/1.80/24], [8.44 m²], [792×280], [65°/2.19/16], [150.1 m³/h · \$71.5k],
      [Dist 110 °C Herm.], [`CO2/R-600 19/81` (`1S`)], [8.33 m²], [938×227], [65°/1.82/16], [9.67 m²], [948×177], [65°/2.18/24], [83.1 m³/h · \$62.5k],
      [Chem 135 °C Herm.], [`CO2/R-600 23/77` (`IHX`)], [10.03 m²], [949×338], [65°/2.01/13], [7.97 m²], [808×140], [65°/1.80/30], [54.5 m³/h · \$59.0k],
      [Chem 135 °C Safe], [`R-290/601 72/28` (`IHX`)], [10.01 m²], [723×328], [62°/1.80/18], [7.99 m²], [700×170], [65°/1.80/28], [66.7 m³/h · \$61.4k],
    ),
    caption: [Appendix A (cont.): Regime B 100 kW#sub[th] AISI 316L BPHE & compressor specifications ($L times W$ and gap in mm).]
  )
]

#heading(numbering: none)[AI-Assisted Methodology Disclosure]
This computational study was directed, reviewed, and published by the author using an autonomous multi-agent research pipeline to synthesize literature, implement multiphysics solvers, execute parallel evolutionary co-design, and conduct independent verification audits. All source code, CAD specs, and Dual-EOS ledgers are openly published at `github.com/theguysaccount/open-physics-autoresearch`.

#heading(numbering: none)[References]
#set text(size: 7.4pt)
#set par(leading: 0.48em, spacing: 0.42em)
[1] H. Lorenz, *Z. Gesamte Kälte-Ind.* 1, 62–69 (1894). \
[2] A. Lorenz & K. Meutzner, *Proc. XIV Int. Cong. Refrigeration*, Moscow (1975). \
[3] G. Lorentzen, *Int. Patent WO 90/07683* (1990). \
[4] K. J. Bell & M. A. Ghaly, *AIChE Symp. Ser.* 69(131), 72–79 (1973). \
[5] D. Tondeur & E. Kvaalen, *Ind. Eng. Chem. Res.* 26(1), 50–56 (1987). \
[6] M.-H. Kim & C. W. Bullard, *Int. J. Refrig.* 25(3), 390–400 (2002). \
[7] H. Müller-Steinhagen & K. Heck, *Chem. Eng. Process.* 20(6), 297–308 (1986). \
[8] S. Z. Rouhani & E. Axelsson, *Int. J. Heat Mass Transfer* 13(2), 383–393 (1970). \
[9] M. M. Shah, *Int. J. Heat Mass Transfer* 22(4), 547–556 (1979). \
[10] A. Muley & R. M. Manglik, *J. Heat Transfer* 121(1), 110–117 (1999). \
[11] D.-H. Han, K.-J. Lee & Y.-H. Kim, *Appl. Therm. Eng.* 23(10), 1209–1225 (2003). \
[12] Y.-Y. Yan, H.-C. Lio & T.-F. Lin, *Int. J. Heat Mass Transfer* 42(6), 993–1006 (1999). \
[13] I. H. Bell et al., *Ind. Eng. Chem. Res.* 53(6), 2498–2508 (2014). \
[14] Ø. Wilhelmsen et al., *Ind. Eng. Chem. Res.* 56(13), 3503–3515 (2017). \
[15] IEC 60335-2-40:2022 (Ed. 7), *Household and Similar Electrical Appliances — Heat Pumps*. \
[16] EN 14825:2022, *Air Conditioners and Heat Pumps — Seasonal Performance Rating*.
