"""
Consolidated Design Comparison and Iteration Candidate Table Updater
====================================================================
Maintains multi-generational design tracking across:
- Baseline 0 (Initial UO2 lattice)
- Baseline 1 (Initial UN/YH2/BeO/Nb-1Zr NEP model)
- Opt-2 (5.5-Year Balanced Reference)
- Opt-3 (First 10-Year Full-Power 2.0 MWth Design, 670.5 kg)
- Opt-4 (Mass-Optimized 10-Year Optimum, 622.7 kg)
- Opt-5 (Ultra-Compact 400 kWth / 100 kWe Space Reactor, ~514 kg)
"""

import os
import csv

def update_comparison_tables():
    comp_csv = "results/design_comparison.csv"
    with open(comp_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Generation", "Designation", "Thermal Power [kWth]", "Electrical Power [kWe]",
            "Pins", "Core Radius [cm]", "Core Height [cm]", "Total Mass [kg]", "Total Volume [L]",
            "k_out (BOL)", "k_in (BOL)", "SDM [pcm]", "Stuck SDM [pcm]", "Total TMC [pcm/K]",
            "Void drho [pcm]", "10-Yr Status", "Achievable Lifetime [yr]", "Key Innovation"
        ])
        writer.writerow([
            "Gen 0", "Baseline 0", "N/A", "N/A",
            "1", "0.40", "N/A", "N/A", "N/A",
            "1.3421", "N/A", "N/A", "N/A", "N/A",
            "N/A", "N/A", "N/A", "Infinite UO2 pin-cell validation"
        ])
        writer.writerow([
            "Gen 1", "Baseline 1", "2000.0", "N/A",
            "397", "18.00", "48.00", "638.30", "154.20",
            "1.0242", "0.9671", "3400", "N/A", "-2.30",
            "+120", "Failed (<2 yr)", "< 2.0", "First 3D UN/YH2/BeO space reactor concept"
        ])
        writer.writerow([
            "Gen 2", "Opt-2", "2000.0", "N/A",
            "349", "17.00", "45.00", "619.50", "148.10",
            "1.0315", "0.9680", "3305", "1250", "-2.10",
            "+110", "Failed (5.5 yr)", "5.5", "6-drum rotation & enhanced B4C sectors"
        ])
        writer.writerow([
            "Gen 3", "Opt-3", "2000.0", "N/A",
            "349", "17.70", "46.20", "670.53", "162.39",
            "1.0623", "0.9749", "2575", "1624", "-2.33",
            "+116", "MET (10.0 yr)", "11.2", "First successful 10-year full-power core"
        ])
        writer.writerow([
            "Gen 4", "Opt-4", "400.0", "100.0",
            "325", "17.35", "45.60", "622.68", "145.23",
            "1.0574", "0.9718", "2902", "932", "-2.85",
            "-56", "MET (10.0 yr)", "> 30.0", "Mass-optimized reflector albedo & 0.5mm coolant gap"
        ])
        writer.writerow([
            "Gen 5", "Opt-5", "400.0", "100.0",
            "277", "16.00", "43.50", "514.08", "120.86",
            "1.0205", "0.9329", "7195", "5835", "-2.48",
            "+165", "MET (10.0 yr)", "20.0", "Ultra-compact core tailored for 400 kWth NEP system (-108.6 kg, -17.4%)"
        ])
    print(f"Updated: {comp_csv}")

if __name__ == "__main__":
    update_comparison_tables()
