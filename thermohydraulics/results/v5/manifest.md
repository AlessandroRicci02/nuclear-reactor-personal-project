# V5 MULTIPHYSICS OPTIMIZATION CAMPAIGN MANIFEST

**Project Directory:** `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Results Directory:** `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5/`  
**Execution Timestamp:** October 3, 2026  
**Status:** COMPLETE, AUDITED, AND REPRODUCIBLE  

---

## 1. Directory Structure

```text
results/v5/
├── V5_FINAL_REPORT.md                 # Primary comprehensive engineering report
├── v5_summary.json                    # Machine-readable high-level campaign summary
├── manifest.md                        # Complete file index and verification manifest
├── v5_baseline_design_vector.json     # Reconstructed Opt-17 baseline parameters
├── v5_local_sensitivities.csv         # 36 perturbation points around baseline
├── v5_derivatives_summary.csv         # First-order local and normalized derivatives
├── v5_screening_population.csv        # 90 LHS & targeted DOE candidate evaluations
├── v5_pareto_front.csv                # 26 Rank-1 non-dominated Pareto designs
├── v5_selected_finalists_summary.csv  # Trade-off selection table for the 5 finalists
├── v5_finalists_comparison.csv        # High-fidelity OpenMC + THM comparison table
├── v5_robustness_analysis.csv         # 78 off-normal operational perturbation cases
├── candidates/                        # JSON records for each DOE candidate (V5_CAND_*.json)
├── finalists/                         # JSON records for each finalist (BASELINE, MIN_MASS, etc.)
│   ├── BASELINE_OPT17.json
│   ├── V5_CAND_MIN_MASS.json
│   ├── V5_CAND_MIN_TEMP.json
│   ├── V5_CAND_MIN_DP.json
│   ├── V5_CAND_BALANCED.json
│   └── V5_CAND_MAX_MARGIN.json
├── pareto/                            # Pareto sorting algorithm exports
│   └── pareto_summary.json
└── figures/                           # 10 High-resolution (300 DPI) publication figures
    ├── v5_01_pareto_front_mass_vs_tclad.png
    ├── v5_02_pareto_front_mass_vs_deltap.png
    ├── v5_03_pareto_front_tclad_vs_deltap.png
    ├── v5_04_orificing_exit_temp_comparison.png
    ├── v5_05_finalists_radar_tradeoff.png
    ├── v5_06_local_sensitivity_tornado.png
    ├── v5_07_robustness_thermal_margins.png
    ├── v5_08_core_mass_breakdown_finalists.png
    ├── v5_09_lifetime_burnup_trajectory.png
    └── v5_10_design_space_correlation_matrix.png
```

---

## 2. Optimization Software Modules (`optimization/`)

1. **`optimization/design_space.py`**
   - Formulates the 10 bounded continuous/discrete design variables.
   - Formalizes the 9 hard engineering constraints and their authoritative literature sources.
   - Defines the 5 objective functions ($F_1$ mass, $F_2$ peak clad temp, $F_3$ core pressure drop, $F_4$ 3D power peaking, $F_5$ hot-channel elevation).
   - Reconstructs the exact Opt-17 Cand 17-Clean baseline design vector.

2. **`inputs/build_opt17_v5.py`**
   - Parametric OpenMC model generator with variable $r_f$, $t_{\mathrm{clad}}$, $t_{\mathrm{gap}}$, pitch, active height, reflector thicknesses, enrichment zoning, and control drum angles.

3. **`optimization/orificing.py`**
   - Solves multi-channel compressible flow equalization across 217 channels in 9 concentric hexagonal rings.
   - Implements 4 flow distribution regimes: Uniform, Ideal Power-Proportional, 3-Zone, and Optimal 9-Ring.
   - Computes orifice loss coefficients $K_p$, orifice plate hole sizing $d_{\mathrm{orif}}$, and pumping power.

4. **`optimization/evaluator.py`**
   - Multi-fidelity candidate evaluator (Level 1 calibrated surrogate, Level 2 intermediate OpenMC, Level 3 high-statistics OpenMC).
   - Exact CAD component mass accounting (fuel, clad, moderator, BeO reflector, $B_4C$ absorber, vessel).
   - Sourced engineering constraints audit.

5. **`optimization/local_sensitivity.py`**
   - Evaluates $\pm 5\%$ and $\pm 10\%$ perturbations around the Opt-17 baseline.
   - Calculates local central-difference derivatives and non-dimensional normalized sensitivities $(dY/Y)/(dX/X)$.

6. **`optimization/generator.py`**
   - Latin Hypercube Sampling (LHS) design-of-experiments generator ($N=90$ candidates).
   - Executes Level 1 screening and filters feasible vs infeasible candidates against the 9 constraints.

7. **`optimization/pareto.py`**
   - Implements fast non-dominated sorting (NSGA-II concept) across the 5 objectives.
   - Extracts the Rank-1 non-dominated Pareto front (26 designs).
   - Identifies the 5 flagship engineering trade-off finalists.

8. **`optimization/run_finalists.py`**
   - Executes high-fidelity multi-physics verification on all 6 reactors (the 5 finalists plus Opt-17 baseline).
   - Evaluates operating criticality, 10-drum shutdown margin, stuck-drum single failure criterion, Doppler reactivity feedback, and multi-channel compressible thermal-hydraulics.

9. **`optimization/robustness.py`**
   - Evaluates 78 off-normal perturbation cases ($\pm 10\%$ mass flow, $\pm 20\text{ K}$ inlet temp, $\pm 5\%$ power, $\pm 25\%$ gap conductance, $\pm 15\%$ HTC, $\pm 10\%$ pressure).
   - Verifies that positive thermal margins are strictly preserved under all credible operational off-normal states.

10. **`optimization/plot_v5.py`**
    - Generates 10 publication-quality figures at 300 DPI with publication-ready styling.

---

## 3. High-Fidelity Deliverables Summary

| Role | Subassembly Mass [kg] | Peak Clad Temp [K] | Core $\Delta P$ [kPa] | $\Delta T_{\mathrm{hot}}$ [K] | Operating $k_{\mathrm{eff}}$ | SDM (All-In) [pcm] | SDM (Stuck) [pcm] | Doppler Coeff [pcm/K] | Full Lifetime [yr] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **BASELINE_OPT17** | 539.66 | 1273.0 | 24.89 | 65.3 | 1.01985 | 2790.5 | 2391.7 | -1.650 | 28.5 |
| **V5_CAND_MIN_MASS** | **517.97** | 1219.2 | 38.57 | 18.8 | 1.01287 | **3865.9** | **3067.0** | -1.002 | 18.6 |
| **V5_CAND_MIN_TEMP** | 553.99 | **1185.3** | 39.26 | **0.0** | 1.02190 | 2634.7 | 2136.6 | -0.194 | 31.4 |
| **V5_CAND_MIN_DP** | 538.76 | 1287.8 | **13.00** | 64.2 | 1.02461 | 2087.5 | 1561.3 | +1.512 | 35.2 |
| **V5_CAND_BALANCED** | **526.18** | **1220.1** | **29.96** | **18.5** | 1.01022 | **3716.9** | **3297.4** | -0.150 | 14.8 |
| **V5_CAND_MAX_MARGIN** | 612.30 | 1224.8 | 35.39 | 19.8 | **1.02782** | 1931.3 | 1613.5 | **-4.412** | **39.7** |

---

## 4. Verification and Reproducibility

All calculations were performed natively on Linux using:
- **OpenMC 0.16.0**: `/home/aleri/miniforge3/envs/openmc-env/bin/openmc`
- **Python 3.13**: `/home/aleri/miniforge3/envs/openmc-env/bin/python3` (with NumPy 2.2.6, SciPy 1.18.1, Pandas 2.2.3, Matplotlib 3.10.0)
- **Nuclear Data**: Windowed Multipole data at `/home/aleri/wmp_data/cross_sections_wmp.xml`
- **Cardinal/MOOSE**: Verified object registries for FlowChannel1Phase and solid conduction benchmarks.
