# Preprint & DOI Submission Guide

Both papers are already live publicly on GitHub with downloadable PDFs and arXiv LaTeX bundles in Release `v1.0.0`. Use the copy-paste packets below to mint an instant **Zenodo DOI** (no endorsement needed, 60 seconds) and/or upload to **arXiv** / **TechRxiv**.

---

## 1. Instant Citable DOI via Zenodo (Recommended First — 60 seconds)

1. Sign into [zenodo.org](https://zenodo.org) with your GitHub account (`theguysaccount`).
2. Visit **Settings → GitHub** (`https://zenodo.org/account/settings/github/`), toggle ON `theguysaccount/open-physics-autoresearch`.
3. Create a quick tag or re-save release `v1.0.0` (or upload the two PDFs directly via `https://zenodo.org/uploads/new`) — Zenodo automatically reads `CITATION.cff` and mints a permanent, Google-Scholar-indexed DOI (`10.5281/zenodo.xxxxxxx`).

---

## 2. Paper 1 — arXiv (`physics.app-ph`, cross-list `eess.SY`) / *Applied Energy*

- **Upload Bundle:** `01-zero-pfas-heat-pumps/paper/lorenz_arxiv_submission.tar.gz` (contains `main.tex`, `references.bib`, and `figures/*.png`)
- **Compiled PDF (for SSRN / Zenodo / ResearchGate):** `01-zero-pfas-heat-pumps/paper/lorenz_autoresearch_preprint.pdf`
- **Title:** `Autonomous Multi-Layer Co-Design of Zero-PFAS Heat Pumps Across Residential Cooling and Industrial Process Heat (<150 C)`
- **Comments:** `5 pages, 2 multi-panel figures, 7 tables. Code, CAD specs, and 8,526-evaluation HEOS/PC-SAFT ledgers at https://github.com/theguysaccount/open-physics-autoresearch`
- **Abstract:**
  > Vapor-compression cooling, space heating, and low-temperature (<200 C) industrial process heat account for ~20,000 TWh/yr of global thermal demand, yet commercial heat pumps achieve only 30-37% of ideal Lorenz efficiency while remaining dependent on high-GWP HFCs or persistent PFAS-class HFOs. Using a locked four-stage autonomous discovery pipeline (8,526 exact CoolProp HEOS finite-volume heat-exchanger simulations cross-checked against SINTEF ThermoPack PC-SAFT, 60-80 segment mesh convergence, 20% Rayleigh vapor-leak fractionation, and 64-draw correlated Monte Carlo), we uncover a sharp two-regime physical phase boundary. In Regime A (7 kW reversible split AC/HP, 6-bin ISO 5151 46 C + EN 14825 -7 C, <=500 g A3 charge), high-glide blends and flash-tank vapor injection lose under Silver-Bell-Ghaly resistance, wet-coil latent flattening, and the 500 g charge cap; instead, an asymmetric 3-row open-fin indoor / 2-row tight-port outdoor microchannel architecture with passive reversible counterflow bridge headers (+15.4% SCOP) and IHX_SUBCOOL charged with R-1270/R-290 (14.8/85.2, 479 g) reaches ASPF = 5.168 (+40.3% over R-410A, +32.1% over R-32, +24.9% over textbook R-290 FTVI+IHX, -0.24% 20% leak shift). In Regime B (100 kW_th industrial waste-heat recovery at 90/110/135 C, fixed 18 m^2 AISI 316L BPHE, P_suc >= 1.05 bar ATEX guard), single-stage natural zeotropes shatter autoresearch-optimized two-stage PFAS R-1233zd(E) ceilings via isenthalpic flash-gas glide skipping: CO2/R-600 reaches COP = 4.728 (+32.5%) at 110 C and 4.420 (+30.3%) at 135 C (LCOH $32.02-$32.49/MWh_th) for factory-sealed hermetic skids, while matched-volatility hydrocarbons (R-290/R-601a 48/52 at 90 C; R-290/R-601 55/45 at 110 C and 72/28 at 135 C) reach COP = 5.327 (+26.1%), 4.076 (+14.2%), and 4.140 (+22.1%) with <1% PC-SAFT agreement and -2.8% to +0.6% shift under a 20% vapor leak.

---

## 3. Paper 2 — arXiv (`eess.SY`, cross-list `physics.app-ph`) / TechRxiv / *IEEE Trans. Energy Conversion*

- **Upload Bundle:** `02-rare-earth-free-motors/paper/motor_arxiv_submission.tar.gz` (contains `main.tex`, `references.bib`, and `figures/motor_autoresearch_proof.png`)
- **Compiled PDF (for TechRxiv / Zenodo / ResearchGate):** `02-rare-earth-free-motors/paper/motor_autoresearch_preprint.pdf`
- **Title:** `Autonomous Multiphysics Co-Design of Rare-Earth-Free Traction Motors Coupled to SiC Inverter kVA, Centrifugal Elasticity, and Cold-Soak Demagnetization`
- **Comments:** `3 pages, 1 multi-panel FEA figure, 3 tables. Gmsh/Numba nonlinear 2D FEA code, CAD specs, and 1,296-evaluation ledgers at https://github.com/theguysaccount/open-physics-autoresearch`
- **Abstract:**
  > Electric vehicle and industrial traction drives remain dependent on Dy-doped sintered NdFeB magnets (~$115/kg) because rare-earth-free alternatives frequently fail under three coupled constraints omitted from simplified studies: low power factor (cos phi) bloating 400 V SiC inverter kVA cost, 2D centrifugal overspeed fracture at 14,400 rpm (1.2x max speed), and cold-soak (-30 C) irreversible demagnetization of ceramic ferrite. Executing 1,296 exact nonlinear 2D Maxwell FEA (Gmsh + Numba Newton-Raphson with node-matching sliding-band rotation) and 2D plane-stress centrifugal elasticity evaluations across 16 CPU cores inside a locked 4.40 L envelope (200 mm OD x 140 mm stack, 400 V bus, 0.5 mm cut-edge lamination degradation, AC proximity copper loss, and 3D end-winding R/L), we map the true rare-earth-free phase boundary (verified at 2x fine mesh, 24 rotation steps, <=1.44% Arkkio-vs-Virtual-Work error, 0.00% frozen-permeability closure, and 32-draw tolerance Monte Carlo). Pure magnet-free SynRM (4p/48s 4B) hits an unbreakable peak-torque and inverter-kVA wall (162.1 Nm, cos phi = 0.652, 178.5 kVA, $1,747 total cost). Conversely, autonomous co-design discovers a 6p/54s 2-barrier enlarged-rotor (137.05 mm OD) asymmetric thick-barrier FB12B La-Co ceramic ferrite traction motor ($12 magnet BOM) deriving 73.0% reluctance and 27.0% ferrite torque that achieves 93.10% WLTC system efficiency, cos phi = 0.809, 99.8 kVA inverter rating, 9.46% WLTC ripple (6.3% rated OP2, with rotor asymmetry alone cutting ripple by 51.5%), SF = 1.52 at 14,400 rpm, +0.019 T cold-soak (-30 C) demagnetization margin, and $1,153 total drive cost (-15.4% below commercial 8p NdFeB $1,363, and within $25 / +2.2% of autoresearch-optimized 6p/54s Dy-NdFeB $1,128—crossing over below $6.46/kVA SiC or above $135/kg NdFeB).
