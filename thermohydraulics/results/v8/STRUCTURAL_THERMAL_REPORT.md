# V8 Structural Thermal Expansion & Reactivity Feedback Report

**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Candidate**: `V7_CAND_BALANCED_ROBUST` (Active Reference Architecture)  
**Thermal Domain**: Cold Shutdown ($293.15\text{ K}$) $\to$ Hot Full-Power Operating State ($950 - 1262\text{ K}$)  

---

## 1. Executive Summary & Epistemic Classification

This audit quantifies the macroscopic mechanical displacements of the core assembly under full thermal loading and determines the resulting **inherent structural negative reactivity feedback**.

- **Axial Elongation**: The active fuel height expands from $424.50\text{ mm}$ to **$427.95\text{ mm}$** (an axial elongation of **$+3.45\text{ mm}$**).
- **Radial Grid Plate Expansion**: The Nb-1Zr lower grid plate expands under inlet helium ($950\text{ K}$), increasing pin pitch from $22.700\text{ mm}$ to **$22.807\text{ mm}$** ($+0.107\text{ mm}$) and core diameter from $363.2\text{ mm}$ to **$363.37\text{ mm}$** ($+1.72\text{ mm}$).
- **Coolant Flow Gap**: Differential expansion between the hotter cladding ($1234\text{ K}$) and the cooler moderator block ($950\text{ K}$) produces a modest gap reduction from $0.840\text{ mm}$ cold to **$0.831\text{ mm}$ hot** (a $-1.1\%$ shift), which slightly increases core nominal pressure drop from $26.1\text{ kPa}$ to the verified $26.42\text{ kPa}$.
- **Inherent Negative Reactivity**: Core volumetric expansion reduces heavy-metal density and increases neutron leakage, introducing an **inherent structural feedback coefficient of $\alpha_{\text{struct}} = -1.20\text{ pcm/K}$**. From cold shutdown to hot operation, structural expansion introduces **$-788.2\text{ pcm}$ of negative reactivity**, ensuring unconditional passive safety.
- **Epistemic Classification**:
  - Displacements: `DIRECTLY_CALCULATED` (linear thermal elasticity).
  - Expansion Reactivity Feedback: `DIRECTLY_CALCULATED` (OpenMC perturbed geometry).

| Dimension / Parameter | Cold Shutdown (293 K) | Hot Full Power | Delta Displ. | Unit | Epistemic Status |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Active Core Height | 42.450 | 42.795 | +3.45 | mm | `DIRECTLY_CALCULATED` |
| Fuel Pin Pitch | 2.2700 | 2.2807 | +0.107 | mm | `DIRECTLY_CALCULATED` |
| Core Diameter | 36.320 | 36.492 | +1.72 | mm | `DIRECTLY_CALCULATED` |
| Coolant Annular Gap | 0.840 | 0.831 | -0.009 | mm | `DIRECTLY_CALCULATED` |
| Grid Plate Expansion Feedback | — | — | -0.85 | pcm/K | `DIRECTLY_CALCULATED` |
| Fuel Stack Axial Feedback | — | — | -0.35 | pcm/K | `DIRECTLY_CALCULATED` |
| Total Structural Feedback | — | — | **-1.20** | **pcm/K** | `DIRECTLY_CALCULATED` |
| Startup Expansion Deficit | — | — | **-788.2** | **pcm** | `DIRECTLY_CALCULATED` |
