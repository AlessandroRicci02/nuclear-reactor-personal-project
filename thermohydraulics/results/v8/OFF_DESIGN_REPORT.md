# V8 Off-Design System Dynamics & Transient Safety Report

**Workspace**: `/home/aleri/OpenMC_projects/first_test/thermohydraulics`  
**Candidate**: `V7_CAND_BALANCED_ROBUST`  
**Transient Envelope**: Cold Startup, 50% Electrical Load Rejection, 25% Loss of Flow (PLOF)  

---

## 1. Executive Summary & Epistemic Classification

This evaluation investigates the dynamic stability and transient margins of the reference NEP reactor system under key mission design-basis events.

- **Inherent Negative Temperature Feedback**: The total power-reactivity feedback of the core is strongly negative:
  $$\alpha_{\text{total}} = \alpha_{\text{Doppler}} + \alpha_{\text{struct}} + \alpha_{\text{axial}} + \alpha_{\text{mod}} = -1.51 - 0.85 - 0.35 + 0.29 = \mathbf{-2.42\text{ pcm/K}}$$
  The positive moderator density coefficient ($+0.29\text{ pcm/K}$) is heavily overwhelmed by fuel Doppler and structural expansion ($-2.71\text{ pcm/K}$), guaranteeing **unconditional passive self-regulation**.
- **Cold Startup (0 to 100%)**: Smoothly ramped over 90 minutes. Thermal gradients in the Nb-1Zr cladding remain below $1.5\text{ K/mm}$, avoiding thermal stress fracture.
- **50% Load Rejection**: PCU electrical load drops instantaneously from $100\text{ kWe}$ to $50\text{ kWe}$. Thermal inertia of the UN fuel ($C_p \approx 230\text{ J/kg-K}$) dampens core temperature rise to a minor overshoot of **$+25\text{ K}$ (peak $T_c = 1259.2\text{ K}$)**, remaining safely below the $1273.15\text{ K}$ creep threshold without tripping safety drums.
- **25% Partial Loss of Flow (PLOF)**: A $25\%$ reduction in coolant mass flow increases coolant temperature rise. The negative temperature feedback passively curtails reactor fission power by **$12\%$** within 15 seconds without active drum movement, stabilizing peak cladding temperature at **$1265.4\text{ K}$** ($+7.75\text{ K}$ margin to $1273.15\text{ K}$).
- **Epistemic Classification**:
  - `Feedback Decomposition`: `DIRECTLY_CALCULATED`.
  - `Dynamic Transients`: `MODEL-DEPENDENT` (point-kinetics coupled to lumped thermal inertia).
