# Project Status & Tracking: Space NEP Reactor

**Current Date / Milestone**: October 2, 2026  
**Current Reference Design**: **Opt-16** (Optimized Hexagonal Architecture C with Dual-Pass Counter-Current Gas Cooling, 10 High-Authority Control Drums, **$540.44\text{ kg}$**, $163.96\text{ L}$, $\text{SDM} = 2832\text{ pcm}$, $\text{SDM}_\text{stk} = 2394\text{ pcm}$, Drum Worth $= 4741\text{ pcm}$, **Strictly Negative Void Reactivity at $2\sigma$**: $\Delta\rho_\text{void} = -87.1 \pm 37.2\text{ pcm}$, **$28.0\text{ yr}$ Full-Power Lifetime**)  
**Audited Benchmark Designs**: **Opt-15** (Trimmed Hexagonal Study), **Opt-14** (Hexagonal Prototype), **Opt-13** (Square Lattice Baseline), **Opt-12** (Canned Hydride Benchmark)  
**System Rating**: **$400\text{ kWth}$** Reactor Thermal Power / **$100\text{ kWe}$** Electrical Bus Power (~25% Closed Brayton efficiency)  

---

## 1. Multi-Generational Key Performance Indicators (Opt-13 through Opt-16)

| Metric | Target / Requirement | Opt-13 Benchmark Value | Opt-14 Hex Prototype | Opt-15 Trimmed Candidate | Opt-16 Reference Candidate | Opt-17 No-Regression Audit | Compliance & Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Architecture Topology** | Non-square Dual-Pass | Square | Hexagonal (Arch C) | Hexagonal (Arch C Opt) | **Hexagonal (Arch C Opt)** | **OPTIMIZED HEX** |
| **Manufacturing Feasibility** | 100% Fabricable Materials | Canned $\text{YH}_2$ | 100% Fabricable | 100% Fabricable | **100% FABRICABLE** | **FULL INTEGRITY** |
| **Total Assembly Mass** | Space envelope ($< 650\text{ kg}$) | $565.10\text{ kg}$ | $580.29\text{ kg}$ | $530.30\text{ kg}$ | **$540.44\text{ kg}$** | **CONFORMING ($+109.56\text{ kg}$ margin)** |
| **Total Assembly Volume**| Compact ($< 170\text{ L}$) | $164.50\text{ L}$ | $163.96\text{ L}$ | $161.22\text{ L}$ | **$163.96\text{ L}$** | **COMPACT FLIGHT ENVELOPE** |
| **Active Pins in Core** | Integer symmetric grid | $377$ | $217$ (9 rings) | $217$ (9 rings) | **$217$ (9 rings)** | **OPTIMIZED HEX LATTICE** |
| **Cooling Topology** | Counter-Current Dual-Pass | Single pass | Dual-Pass Counter-Current | Dual-Pass Counter-Current | **Dual-Pass Counter-Current** | **ISOTHERMAL MODERATOR** |
| **Axial Streaming Ducts**| Penetrating BeO reflector | $0.820\text{ cm}$ | $0.820\text{ cm} + 6\text{ aux/pin}$ | $0.835\text{ cm} + 6\text{ aux/pin}$ | **$0.820\text{ cm} + 6\text{ aux/pin}$** | **ENHANCED STREAMING** |
| **Control System Architecture**| Drums in Radial Reflector | 10 drums | 10 drums ($140^\circ$) | 10 drums ($142^\circ$) | **10 drums ($145^\circ$, 90% B-10)** | **HIGH AUTHORITY** |
| **BOL Multiplication ($k_\text{op}$)** | $\ge 1.005$ at 400 kWth | $1.01501$ | $1.02007 \pm 0.00051$ | $1.01167 \pm 0.00034$ | **$1.01994 \pm 0.00027$** | **CONFORMING ($+1909.1\text{ pcm}$)** |
| **Shutdown Multiplication ($k_\text{in}$)**| $\le 0.980$ (SDM $\ge 2000\text{ pcm}$) | $0.95552$ | $0.97303 \pm 0.00070$ | $0.96301 \pm 0.00082$ | **$0.97246 \pm 0.00063$** | **CONFORMING ($\text{SDM} = 2831.8\text{ pcm}$)** |
| **1-Stuck Drum Margin ($\text{SDM}_\text{stk}$)**| $\ge 1000\text{ pcm}$ | $3859\text{ pcm}$ | $2425.2\text{ pcm}$ | $3244.0\text{ pcm}$ | **$2393.7\text{ pcm}$** | **CONFORMING ($>2.39\times$ margin)** |
| **Total Control Drum Worth** | $\ge 4000\text{ pcm}$ | $6133\text{ pcm}$ | $4739.0\text{ pcm}$ | $5025.6\text{ pcm}$ | **$4740.9\text{ pcm}$** | **CONFORMING ($>4700\text{ pcm}$)** |
| **Fuel Doppler Feedback ($\alpha_\text{fuel}$)** | Prompt Negative ($3\sigma < 0$) | $-1.56\text{ pcm/K}$ | $-1.08 \pm 0.18\text{ pcm/K}$ | $-0.88 \pm 0.17\text{ pcm/K}$ | **$-1.11 \pm 0.17\text{ pcm/K}$** | **ROBUSTLY NEGATIVE ($3\sigma \le -0.61\text{ pcm/K} < 0$)** |
| **Moderator Feedback ($\alpha_\text{mod}$)**| Strictly Negative ($3\sigma < 0$) | $-0.94\text{ pcm/K}$ | $-1.47 \pm 0.18\text{ pcm/K}$ | $-1.42 \pm 0.21\text{ pcm/K}$ | **$-1.60 \pm 0.19\text{ pcm/K}$** | **ROBUSTLY NEGATIVE ($3\sigma \le -1.02\text{ pcm/K} < 0$)** |
| **Coolant Feedback ($\alpha_\text{cool}$)** | Non-positive ($1\sigma < 0$) | $-0.033\text{ pcm/K}$ | $-0.061 \pm 0.042\text{ pcm/K}$ | $-0.014 \pm 0.045\text{ pcm/K}$ | **$-0.083 \pm 0.035\text{ pcm/K}$** | **ROBUSTLY NEGATIVE ($2\sigma \le -0.012\text{ pcm/K} < 0$)** |
| **Total Temperature Coeff ($\alpha_\text{total}$)**| Negative ($3\sigma < 0$) | $-2.54\text{ pcm/K}$ | $-2.61 \pm 0.26\text{ pcm/K}$ | $-2.32 \pm 0.28\text{ pcm/K}$ | **$-2.79 \pm 0.26\text{ pcm/K}$** | **ROBUSTLY NEGATIVE ($3\sigma \le -2.02\text{ pcm/K} < 0$)** |
| **Coolant Void Reactivity ($\Delta\rho_\text{void}$)** | Strictly Negative ($1\sigma, 2\sigma < 0$)| $-34.3\text{ pcm}$ | $-64.5 \pm 43.8\text{ pcm}$ | $-15.0 \pm 47.0\text{ pcm}$ | **$-87.1 \pm 37.2\text{ pcm}$** | **STRICTLY NEGATIVE AT $2\sigma$ ($2\sigma \le -12.6\text{ pcm} < 0$)** |
| **Achievable Full-Power Lifetime**| $\ge 10.0\text{ years}$ verified | $20.4\text{ yr}$ | $28.8\text{ yr}$ | $17.4\text{ yr}$ | **$28.0\text{ YEARS}$** | **CONFORMING ($>2.8\times$ margin)** |

---

## 2. Principal Geometric & Material Parameters (Candidate Opt-13 Reference)

* **Core Radius ($R_\text{core}$)**: $19.45\text{ cm}$
* **Active Core Height ($H_\text{core}$)**: $42.20\text{ cm}$ ($H/R = 2.170$, $H/D = 1.085$)
* **Fuel Element**: Sintered UN pellet ($r_\text{fuel} = 0.580\text{ cm}$, $19.75\text{ wt\%}$ $^{235}\text{U}$), Nb-1Zr cladding ($t_\text{clad} = 0.050\text{ cm}$)
* **Coolant Channel Gap**: $t_\text{cool} = 0.075\text{ cm}$ ($0.75\text{ mm}$ annular gas passage, pure Helium)
* **Axial Streaming Ducts**: Open non-overlapping ducts of radius $r_\text{duct} = 0.820\text{ cm}$ penetrating top and bottom BeO reflectors
* **Pin Lattice**: 377 fuel pins on square lattice, moderator pitch $P = 1.700\text{ cm}$ ($V_\text{mod}/V_\text{fuel} = 1.505$)
* **Coolant Working Fluid**: Pure Helium ($^4\text{He}$) at $P = 2.5\text{ MPa}$, $\bar{T} = 1050\text{ K}$, $\rho = 0.00115\text{ g/cm}^3$ (Zero parasitic Xenon poison)
* **Moderator Block**: Sintered stoichiometric $\text{YH}_{2.0}$ blocks ($\rho = 4.30\text{ g/cm}^3$) hermetically canned in refractory Nb-1Zr with internal CVD $\text{Al}_2\text{O}_3$ hydrogen permeation barrier ($P_{\text{H}_2}(800\text{ K}) \approx 5\text{ mbar} \ll 25\text{ bar}$)
* **Reflector Sizing**: Pure sintered BeO blocks: radial thickness $9.80\text{ cm}$ ($R_\text{outer} = 29.25\text{ cm}$); axial BeO blocks $9.50\text{ cm}$ each ($H_\text{total} = 61.20\text{ cm}$)
* **Control Drums**: 10 rotating BeO drums ($R_\text{drum} = 3.95\text{ cm}$) with $138.0^\circ$ sector of $90\%$ $^{10}\text{B}$-enriched $\text{B}_4\text{C}$ ($t_\text{abs} = 1.40\text{ cm}$)
* **Total Reactor Assembly Mass**: **$565.10\text{ kg}$** ($+84.90\text{ kg}$ safety margin below $650\text{ kg}$ limit)
* **Total Reactor Assembly Volume**: **$164.50\text{ L}$**
* **Depletion Longevity**: **$20.39\text{ years}$** continuous full power at $400\text{ kWth}$ ($k_\text{eff}(10\text{ yr}) = 1.00759 > 1.000$)

---

## 3. Archival Snapshot & Independent Audit Verification
* **Archival Directory**: `best_configurations/opt13/` (complete recursive snapshot, 30 files, 100% SHA-256 match).
* **Audit Directory**: `best_configurations/opt13_audit/` (includes `OPT13_INDEPENDENT_AUDIT.md`, `OPT13_CORRECTED_SUMMARY.txt`, `OPT13_FILE_MANIFEST.csv`, `OPT13_VALIDATION_RESULTS.csv`).
* **Audit Confirmation**: All reported geometry, mass, control margins, Doppler and total TMC coefficients independently re-computed and verified. Statistical limitation on void reactivity ($2\sigma, 3\sigma$ crossing zero) explicitly audited and documented.

---

## 4. Iteration Opt-17: Strict No-Regression Investigation & Pareto Frontier Mapping
* **Directive Objective**: Optimize Opt-16 under a strict no-regression rule where EVERY metric must be no worse than Opt-16 with at least one meaningful improvement.
* **Candidates Evaluated**: 11 configurations spanning geometric fuel zoning, radial enrichment grading, control drum sector widening, axial streaming duct enlargement, and active core aspect ratio adjustments.
* **Key Findings**:
  - Radial power flattening ($F_q$ reduced from $3.428$ to $3.283$) and enhanced void reactivity (down to $-144.9\text{ pcm}$ and $-296.1\text{ pcm}$) were demonstrated.
  - However, every candidate suffered a statistically confirmed regression in operational lifetime ($20.5\text{--}27.2\text{ yr} < 28.0\text{ yr}$) due to central fissile reactivity loss, parasitic reflection degradation at $0^\circ$, or increased continuous operational leakage.
* **Final Program Determination**: Under the strict no-regression rule, **Candidate Opt-16 is reaffirmed and formally retained as the reference flight design**.
* **Deliverables Generated**:
  - Comparison CSV: `results_second_architecture/opt17/candidate_comparison.csv`
  - No-Regression Audit: `results_second_architecture/opt17/no_regression_check.csv`
  - 6 High-Resolution 300 DPI Figures: `results_second_architecture/opt17/figures/`
  - Comprehensive Report: `results_second_architecture/opt17/opt17_comprehensive_report.md`
  - Plain-Text Report: `results_second_architecture/opt17/reports/opt17_comprehensive_report.txt` and `results/report.txt`

