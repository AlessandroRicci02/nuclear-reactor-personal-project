# V7 Campaign File & Artifact Manifest

**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/`  
**Generated**: October 2026  
**Auditor**: Antigravity Autonomous Multiphysics Engineering Agent  

---

## 1. Primary Engineering Reports

| File Path | Description | Format |
| :--- | :--- | :---: |
| [`results/v7/V7_FINAL_REPORT.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/V7_FINAL_REPORT.md) | Comprehensive engineering final report for the V7 campaign | Markdown |
| [`results/v7/V7_FINAL_REPORT.txt`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/V7_FINAL_REPORT.txt) | Plain-text copy of the comprehensive final report (per user directive) | Plain Text |
| [`results/v7/V6_REVALIDATION.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/V6_REVALIDATION.md) | Forensic audit of V6 claims, solver provenance, and candidate reclassification | Markdown |
| [`results/v7/ROBUST_DESIGN_SPACE.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/ROBUST_DESIGN_SPACE.md) | Formulation of robust feasibility and targeted second-stage optimization | Markdown |
| [`results/v7/ORIFICING_ROBUSTNESS.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/ORIFICING_ROBUSTNESS.md) | 9-ring orifice subsystem engineering, manufacturing tolerances, off-design matrix | Markdown |
| [`results/v7/CORRELATION_MODEL_UNCERTAINTY.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/CORRELATION_MODEL_UNCERTAINTY.md) | Model-form uncertainty in the transitional Reynolds regime ($Re \in [1300, 2200]$) | Markdown |
| [`results/v7/BURNUP_LIFETIME_AUDIT.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/BURNUP_LIFETIME_AUDIT.md) | Audit of lifetime extrapolation, depletion physics, and 95% confidence bounds | Markdown |
| [`results/v7/PHYSICS_MATURITY_ASSESSMENT.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/PHYSICS_MATURITY_ASSESSMENT.md) | Assessment of 14 omitted phenomena, margin relevance, and required tools | Markdown |
| [`results/v7/MISSION_ARCHITECTURES.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/MISSION_ARCHITECTURES.md) | Mission-specific architectural specifications and trade-off matrices | Markdown |
| [`results/v7/DESIGN_SELECTION_MATRIX.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/DESIGN_SELECTION_MATRIX.md) | Maturity levels, multi-attribute scoring, and final design parameters | Markdown |
| [`results/v7/v7_summary.json`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/v7_summary.json) | Machine-readable JSON summary of findings and finalist parameters | JSON |
| [`results/v7/manifest.md`](file:///home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v7/manifest.md) | Complete file and artifact directory manifest | Markdown |

---

## 2. Generated Data & CSV Artifacts

| Directory | File Name | Description |
| :--- | :--- | :--- |
| `audit/` | `v6_forensic_metrics_audit.csv` | Metric-by-metric audit of V6 results, solver identity, statistical uncertainty |
| `audit/` | `v6_revalidation_feasibility_matrix.csv` | Reclassification of all 10 candidates against all hard constraints |
| `correlations/` | `correlation_sweep_re.csv` | Parameter sweep across 5 friction and 6 Nusselt models ($Re \in [500, 4500]$) |
| `correlations/` | `candidate_ranking_by_correlation.csv` | Candidate temperature and pressure drop ranking across correlations |
| `orificing/` | `robust_orificing_matrix.csv` | Comparison of 4 flow distribution regimes |
| `orificing/` | `robust_orificing_ring_specs.csv` | Concentric 9-ring orifice hole diameters, loss coefficients, and exit temps |
| `orificing/` | `robust_orificing_perturbations.csv` | 14-case off-design and manufacturing perturbation suite |
| `depletion/` | `burnup_depletion_audit.csv` | Fissile consumption, burnup loss rates, and specific power breakdown |
| `depletion/` | `lifetime_uncertainty_bounds.csv` | 95% confidence upper and lower bounds on full and reserve lifetime |
| `optimization/` | `all_candidates_pareto_comparison.csv` | Complete ledger of all 12 evaluated candidates with nominal and perturbed metrics |
| `optimization/` | `nominal_pareto_front.csv` | Non-dominated candidate set under nominal operating conditions |
| `optimization/` | `robust_pareto_front.csv` | Non-dominated candidate set under worst-case perturbed conditions |
| `optimization/` | `robust_candidates_summary.csv` | Summary table of the 3 primary V7 robust finalists |
| `optimization/` | `V7_BALANCED_ROBUST.json` | Detailed parameter ledger for Active Reference Architecture |
| `optimization/` | `V7_LONG_LIFE_ROBUST.json` | Detailed parameter ledger for Deep-Space Flagship Architecture |
| `optimization/` | `V7_BRAYTON_ROBUST.json` | Detailed parameter ledger for Low-Pressure-Drop Architecture |

---

## 3. High-Quality Publication Figures (300 DPI)

All figures reside in `results/v7/figures/` and are replicated in the conversation artifact directory:

1. `v7_01_nominal_vs_robust_pareto.png` — Nominal vs. Robust Pareto Frontier Under Perturbations
2. `v7_02_mass_vs_lifetime.png` — Reactor Mass vs. Operational Depletion Lifetime
3. `v7_03_mass_vs_pressure_drop.png` — Mass vs. Pressure Drop Trade-Off
4. `v7_04_pressure_drop_vs_peak_clad_temp.png` — Core Pressure Drop vs. Peak Cladding Temperature
5. `v7_05_orifice_tolerance_sensitivity.png` — 9-Ring Orifice Plate Tolerance Sensitivity ($\pm 0.05$ mm)
6. `v7_06_thermal_maldistribution_robustness.png` — Radial Exit Temperature Equalization & Maldistribution Robustness
7. `v7_07_correlation_model_form_uncertainty.png` — Transitional Reynolds Number Model-Form Uncertainty Assessment
8. `v7_08_monte_carlo_uncertainty_vs_margins.png` — Monte Carlo Statistical Uncertainty vs. Constraint C07 Margins
9. `v7_09_burnup_lifetime_uncertainty.png` — Long-Term Reactivity Depletion Trajectory & Uncertainty Envelope
10. `v7_10_mission_specific_feasible_regions.png` — Mission-Specific Architectural Feasibility Envelopes
11. `v7_11_final_architecture_comparison_radar.png` — Final Multi-Physics Architecture Trade-Off Radar

---

## 4. Executable Python Scripts in `v7_scripts/`

| Script Path | Primary Function |
| :--- | :--- |
| `v7_scripts/v6_forensic_revalidation.py` | Audits V6 data, produces forensic metrics and feasibility reclassification |
| `v7_scripts/transitional_reynolds_study.py` | Performs Re sweep and candidate ranking invariance checks |
| `v7_scripts/robust_orificing_study.py` | Sizes 9-ring plate, runs 14-case off-design perturbation suite |
| `v7_scripts/audit_lifetime_depletion.py` | Analyzes depletion physics and propagates lifetime uncertainty bounds |
| `v7_scripts/evaluate_v7_robust_candidates.py` | Evaluates coupled performance of the 3 V7 robust finalists |
| `v7_scripts/generate_pareto_fronts.py` | Computes nominal and robust non-dominated Pareto sets |
| `v7_scripts/generate_v7_figures.py` | Generates all 11 publication figures at 300 DPI and mirrors to artifacts |
