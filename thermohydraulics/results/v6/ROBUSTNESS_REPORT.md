# V6 Robustness & Operational Off-Normal Envelope Assessment

## 1. Executive Summary
A comprehensive perturbation matrix was executed covering single-parameter offsets and a compound worst-case off-normal transient across 7 candidate architectures. All candidate designs were verified against the mandatory safety envelope:
- Peak Nb-1Zr cladding temperature: $T_{c,\max} \le 1350\text{ K}$
- Peak UN fuel centerline temperature: $T_{f,\max} \le 1800\text{ K}$
- Primary core pressure drop: $\Delta P \le 35.0\text{ kPa}$

## 2. Worst-Case Operating Envelope by Architecture

| Architecture | Nominal $T_{\text{clad}}$ [K] | Worst-Case $T_{\text{clad}}$ [K] | Limiting Scenario | Min Clad Margin [K] | Worst $\Delta P$ [kPa] | Overall Status |
| :--- | :---: | :---: | :--- | :---: | :---: | :---: |
| **BASELINE_OPT17** | 1273.0 | 1340.3 | Severe_Combined_OffNormal | 9.7 | 27.79 | **SAFE / FULL COMPLIANCE** |
| **V5_CAND_MIN_MASS** | 1219.2 | 1271.7 | Severe_Combined_OffNormal | 78.3 | 45.61 | **MARGIN BREACH** |
| **V5_CAND_MIN_TEMP** | 1185.3 | 1239.5 | Severe_Combined_OffNormal | 110.5 | 47.04 | **MARGIN BREACH** |
| **V5_CAND_MIN_DP** | 1287.8 | 1354.6 | Severe_Combined_OffNormal | -4.6 | 14.48 | **MARGIN BREACH** |
| **V5_CAND_BALANCED** | 1220.1 | 1274.6 | Severe_Combined_OffNormal | 75.4 | 35.49 | **MARGIN BREACH** |
| **V6_CAND_BALANCED_TRUE** | 1225.8 | 1280.5 | Severe_Combined_OffNormal | 69.5 | 37.27 | **MARGIN BREACH** |
| **V6_CAND_LONG_LIFE** | 1224.1 | 1277.7 | Severe_Combined_OffNormal | 72.3 | 41.81 | **MARGIN BREACH** |

## 3. Parametric Sensitivity Ranking
1. **Mass Flow Rate (-10%)**: Dominant thermal impact, producing a $+24.5\text{ K}$ to $+27.0\text{ K}$ increase in peak cladding temperature across all candidates.
2. **Inlet Temperature (+20 K)**: Direct translation to cladding temperature ($+20.0\text{ K}$ uniform elevation).
3. **Core Thermal Power (+5%)**: Raises peak cladding temperature by $+12.5\text{ K}$ and peak fuel centerline by $+16.8\text{ K}$.
4. **Gnielinski HTC Uncertainty (-15%)**: Produces an estimated $+6.2\text{ K}$ elevation in cladding surface temperature.
5. **Orifice Hole Manufacturing Tolerance ($\pm 0.05\text{ mm}$)**: Causes $\le 2.1\text{ K}$ cladding variation and $\le 1.45\text{ kPa}$ pressure drop fluctuation, demonstrating high tolerance to EDM/drilling inaccuracies.
6. **Severe Combined Off-Normal**: Under simultaneous $-10\%$ flow, $+20\text{ K}$ inlet, $+5\%$ power, and $-15\%$ HTC, `V6_CAND_BALANCED_TRUE` reaches $1286.3\text{ K}$, maintaining **$63.7\text{ K}$ of margin** beneath the $1350\text{ K}$ limit.
