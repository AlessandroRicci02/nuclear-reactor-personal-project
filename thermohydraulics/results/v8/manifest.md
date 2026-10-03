# V8 Campaign Manifest & File Deliverables Directory

**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Campaign**: V8 Independent Verification, Fuel Performance & System-Level Engineering  
**Active Reference Architecture**: `V7_CAND_BALANCED_ROBUST`  

---

## 1. Directory Structure

```text
results/v8/
├── audit/
│   ├── design_genealogy_audit.csv
├── verification/
│   ├── independent_revalidation_summary.csv
│   ├── V7_BALANCED_ROBUST_channels.csv
│   ├── V7_LONG_LIFE_ROBUST_channels.csv
│   └── V7_BRAYTON_ROBUST_channels.csv
├── fuel_performance/
│   ├── fuel_performance_trajectory.csv
│   └── fuel_performance_summary.csv
├── yh2_moderator/
│   ├── yh2_permeation_trajectory.csv
│   └── yh2_loss_sensitivity.csv
├── plenum/
│   └── plenum_maldistribution_channels.csv
├── brayton/
│   ├── brayton_cycle_summary.csv
│   └── brayton_state_points.csv
├── radiator/
│   ├── radiator_summary.csv
│   └── radiator_trade_curves.csv
├── vvuq/
│   └── joint_uncertainty_quantiles.csv
├── figures/
│   ├── v8_01_simulation_space_registry.png
│   ├── v8_02_design_genealogy.png
│   ├── v8_03_revalidation_profiles.png
│   ├── v8_04_fuel_performance_pcmi.png
│   ├── v8_05_yh2_moderator_permeation.png
│   ├── v8_06_plenum_cfd_maldistribution.png
│   ├── v8_07_brayton_cycle_ts_diagram.png
│   ├── v8_08_radiator_sizing_tradeoffs.png
│   ├── v8_09_off_design_transients.png
│   ├── v8_10_joint_uncertainty_distributions.png
│   └── v8_11_mature_architectures_radar.png
├── COMPLETE_SIMULATION_REGISTRY.csv
├── COMPLETE_SIMULATION_REGISTRY.md
├── FUEL_PERFORMANCE_REPORT.md
├── YH2_MODERATOR_REPORT.md
├── PLENUM_CFD_REPORT.md
├── BRAYTON_SYSTEM_REPORT.md
├── RADIATOR_SYSTEM_REPORT.md
├── STRUCTURAL_THERMAL_REPORT.md
├── OFF_DESIGN_REPORT.md
├── VVUQ_MATRIX.md
├── DESIGN_MATURITY_ASSESSMENT.md
├── v8_summary.json
├── manifest.md
├── V8_FINAL_REPORT.md
└── V8_FINAL_REPORT.txt
```

---

## 2. Deliverables Inventory

| File Path | Description | Epistemic Status |
| :--- | :--- | :---: |
| `COMPLETE_SIMULATION_REGISTRY.csv` | Full audit of all 129 runs across V4, V5, V6, and V7 campaigns | `VERIFIED` |
| `COMPLETE_SIMULATION_REGISTRY.md` | Human-readable markdown audit of all runs | `VERIFIED` |
| `audit/design_genealogy_audit.csv` | Multi-generational design evolution table (Opt-16 to V8) | `DIRECTLY_CALCULATED` |
| `verification/independent_revalidation_summary.csv` | Independent 217-channel thermal-hydraulic re-evaluation | `DIRECTLY_CALCULATED` |
| `fuel_performance/fuel_performance_trajectory.csv` | BISON 0-40 yr swelling, gap closure, and creep strain | `DIRECTLY_CALCULATED` |
| `yh2_moderator/yh2_loss_sensitivity.csv` | Hydrogen dissociation, permeation loss, and reactivity | `DIRECTLY_CALCULATED` |
| `plenum/plenum_maldistribution_channels.csv` | 3D header static pressure and channel flow maldistribution | `MODEL-DEPENDENT` |
| `brayton/brayton_cycle_summary.csv` | Closed Brayton cycle component state points and power | `DIRECTLY_CALCULATED` |
| `radiator/radiator_summary.csv` | Deployable radiator area, mass, and minimum-alpha sweep | `DIRECTLY_CALCULATED` |
| `vvuq/joint_uncertainty_quantiles.csv` | Latin Hypercube (N=10,000) joint uncertainty quantiles | `DIRECTLY_CALCULATED` |
| `V8_FINAL_REPORT.md` | Comprehensive 25-section master engineering report | `SYNTHESIZED` |
| `V8_FINAL_REPORT.txt` | Exact plain-text copy of master engineering report | `SYNTHESIZED` |
