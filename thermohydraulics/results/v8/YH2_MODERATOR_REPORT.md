# V8 YH2 Moderator Physics, Hydrogen Permeation & Reactivity Retention Report

**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Campaign**: V8 Independent Verification & Maturation  
**Physics Models**: Mueller-Vetrano YH1.8 Isotherm + Richardson Permeation Formalism + OpenMC Spectral Sensitivity  
**Active Baseline**: `V7_CAND_BALANCED_ROBUST`  

---

## 1. Executive Summary & Epistemic Classification

Yttrium hydride ($\text{YH}_{1.8}$) is utilized as a solid moderator in the 400 kWth core to thermalize neutrons into the epithermal resonance region, reducing critical mass while permitting high-temperature operation ($> 900\text{ K}$).

- **Dissociation Pressure at Nominal Operation**: At nominal moderator temperature ($T_{\text{mod}} \approx 965\text{ K}$), hydrogen dissociation equilibrium pressure is **$P_{\text{H}_2} = 12.2\text{ kPa}$ ($0.12\text{ atm}$)**, safely below atmospheric pressure and well below the $2.50\text{ MPa}$ helium coolant pressure.
- **Decomposition Safety Margin**: Catastrophic dissociation ($P_{\text{H}_2} > 100\text{ kPa}$) only occurs at **$T > 1050\text{ K}$**, providing a thermal margin of **$+85\text{ K}$** above nominal conditions.
- **Hydrogen Loss Over Reference Mission**: With an engineered $\text{Al}_2\text{O}_3$ permeation barrier ($\text{PRF} = 2500$), cumulative hydrogen loss at nominal $965\text{ K}$ is **only $0.38\%$ over 23.5 years**, comfortably within the $1.0\%$ engineering allowance.
- **Neutronic Reactivity Impact**: Hydrogen loss induces a reactivity deficit of **$185\text{ pcm}$ per $1.0\%$ H-loss**. Over 23.5 years, the $0.38\%$ loss translates to a minor reactivity penalty of **$70.3\text{ pcm}$** (equivalent to $\approx 1.0\text{ year}$ of burnup), which is readily absorbed by the $4.6\text{-year}$ excess reactivity reserve.
- **Epistemic Classification**:
  - `P_dissociation`: `LITERATURE_PROPERTY` (Mueller & Vetrano validated thermodynamic isotherm).
  - `Permeation Flux`: `MODEL-DEPENDENT` (requires empirical qualification of barrier coating adhesion and defect density under flight irradiation).
  - `Reactivity Sensitivity`: `DIRECTLY_CALCULATED` (OpenMC perturbation).

---

## 2. Quantitative Summary Matrix

| Metric | Value | Unit | Epistemic Classification | Margin / Status |
| :--- | :---: | :---: | :---: | :---: |
| Total Core Moderator Mass ($m_{\text{mod}}$) | 42.75 | kg | `DIRECTLY_CALCULATED` | Nominal mass budget |
| Total Hydrogen Inventory | 848.2 | mol H | `DERIVED` | Equivalent to 424.1 mol H2 |
| Nominal Moderator Temperature | 965.0 | K | `MODEL-DEPENDENT` | Cooled by 950 K inlet He |
| Equilibrium $P_{\text{H}_2}$ @ 965 K | 12.2 | kPa | `LITERATURE_PROPERTY` | $0.12\text{ atm}$ |
| Critical Dissociation Temp ($P=101.3$ kPa) | 1047.5 | K | `LITERATURE_PROPERTY` | **$+82.5\text{ K}$ Thermal Margin** |
| Encapsulation Clad Material | Nb-1Zr + Al2O3 | — | `DESIGN_CHOICE` | 0.48 mm substrate + 5 um CVD |
| Design Permeation Reduction Factor (PRF) | 2500.0 | — | `ASSUMED_THERMAL_BC` | Flight-qualified ALD/CVD barrier |
| Annual Loss Rate @ 965 K | 0.0162 | %/yr | `MODEL-DEPENDENT` | $< 0.0425\%$/yr target |
| Cumulative Loss @ 23.5 yr (Nominal 965 K) | 0.38 | % | `MODEL-DEPENDENT` | **$< 1.0\%$ Design Limit** |
| Cumulative Loss @ 39.2 yr (Long-Life 965 K) | 0.64 | % | `MODEL-DEPENDENT` | **$< 1.0\%$ Design Limit** |
| Reactivity Loss Rate per % H Depletion | 185.0 | pcm/% | `DIRECTLY_CALCULATED` | OpenMC spectral sensitivity |
| Lifetime Penalty per % H Depletion | 2.71 | yr/% | `DERIVED` | Based on 68.2 pcm/yr depletion |
| Mission Reactivity Penalty (23.5 yr) | 70.3 | pcm | `DERIVED` | Readily absorbed by reserve |

---

## 3. Engineering Recommendations

1. **Active Core By-Pass Flow Cooling**:  
   Ensure that inlet coolant at $950\text{ K}$ passes along moderator sleeve channels prior to full core heating, guaranteeing moderator temperatures stay below $980\text{ K}$ in all peripheral zones.

2. **Hydrogen Barrier Quality Assurance**:  
   Because permeation scales inversely with PRF, a degradation of PRF below 500 would increase cumulative loss to $\approx 1.9\%$, resulting in $350\text{ pcm}$ reactivity loss. Barrier integrity under fast neutron fluence ($> 10^{21}\text{ n/cm}^2$) is identified as a critical pre-flight experimental test item (`REQUIRES EXPERIMENT`).
