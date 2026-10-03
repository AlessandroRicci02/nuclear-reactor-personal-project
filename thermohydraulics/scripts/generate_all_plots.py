"""
PUBLICATION PLOT GENERATOR FOR OPT-17 THERMOHYDRAULICS COUPLED CAMPAIGN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/generate_all_plots.py

Generates all 13 required figures at 300 DPI:
1. keff_vs_iteration.png
2. reactivity_vs_iteration.png
3. peak_power_vs_iteration.png
4. power_convergence.png
5. fuel_temperature_vs_iteration.png
6. clad_temperature_vs_iteration.png
7. coolant_temperature_vs_iteration.png
8. density_feedback_vs_iteration.png
9. axial_power_profile_final.png
10. radial_power_distribution_final.png
11. fuel_temperature_axial_profile.png
12. coolant_temperature_axial_profile.png
13. pressure_drop.png
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_DIR = os.path.join(WORKSPACE_ROOT, "results")
NEUT_DIR = os.path.join(RESULTS_DIR, "neutronics")
THERM_DIR = os.path.join(RESULTS_DIR, "thermal")
COUP_DIR = os.path.join(RESULTS_DIR, "coupling")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")

plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "axes.grid": True,
    "grid.alpha": 0.4,
    "grid.linestyle": "--"
})


def generate_all_13_plots():
    print("=" * 70)
    print("  GENERATING ALL 13 PUBLICATION FIGURES AT 300 DPI")
    print("=" * 70)
    os.makedirs(FIG_DIR, exist_ok=True)
    
    # 1. Load dataframes
    hist_csv = os.path.join(COUP_DIR, "coupling_iterations.csv")
    df_hist = pd.read_csv(hist_csv)
    
    pin_csv = os.path.join(NEUT_DIR, "final_pin_power.csv")
    df_pins = pd.read_csv(pin_csv)
    
    axial_csv = os.path.join(NEUT_DIR, "final_axial_power.csv")
    df_axial = pd.read_csv(axial_csv)
    
    therm_csv = os.path.join(THERM_DIR, "final_thermal_field.csv")
    df_therm = pd.read_csv(therm_csv)
    
    iters = df_hist["iteration"].values
    
    # -------------------------------------------------------------------------
    # 1. keff_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(iters, df_hist["k_eff"], yerr=df_hist["sigma_k"], fmt="o-", color="navy",
                capsize=5, elinewidth=1.5, capthick=1.5, label=r"$k_{\mathrm{eff}} \pm \sigma$")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel(r"Eigenvalue $k_{\mathrm{eff}}$")
    ax.set_xticks(iters)
    ax.set_title("Opt-17 Cand 17-Clean: Eigenvalue Convergence with Monte Carlo Uncertainty", pad=12)
    ax.legend(frameon=True, loc="upper right")
    fig.tight_layout()
    p1 = os.path.join(FIG_DIR, "keff_vs_iteration.png")
    plt.savefig(p1)
    plt.close()
    print(f"  [+] Saved: {p1}")

    # -------------------------------------------------------------------------
    # 2. reactivity_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["reactivity"], "s--", color="crimson", linewidth=2.0, label=r"Doppler Shift $\Delta \rho$ [pcm]")
    ax.axhline(0.0, color="gray", linestyle=":")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel(r"Reactivity Shift $\Delta \rho$ [pcm]")
    ax.set_xticks(iters)
    ax.set_title("Net Doppler & Thermal Feedback Reactivity Shift vs. Iteration", pad=12)
    ax.legend(frameon=True, loc="lower left")
    fig.tight_layout()
    p2 = os.path.join(FIG_DIR, "reactivity_vs_iteration.png")
    plt.savefig(p2)
    plt.close()
    print(f"  [+] Saved: {p2}")

    # -------------------------------------------------------------------------
    # 3. peak_power_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(8, 5))
    color = "purple"
    ax1.set_xlabel("Picard Coupling Iteration")
    ax1.set_ylabel("Peak Pin Power [W]", color=color, fontweight="bold")
    ax1.plot(iters, df_hist["peak_pin_power"], "o-", color=color, linewidth=2.0, label="Peak Pin Power [W]")
    ax1.tick_params(axis="y", labelcolor=color)
    ax1.set_xticks(iters)
    
    ax2 = ax1.twinx()
    color2 = "darkgreen"
    ax2.set_ylabel("Peaking Factors ($F_r$, $F_q$)", color=color2, fontweight="bold")
    ax2.plot(iters, df_hist["F_r"], "s--", color="green", linewidth=1.8, label=r"Radial Peaking $F_r$")
    ax2.plot(iters, df_hist["F_q"], "^-.", color="darkolivegreen", linewidth=1.8, label=r"3D Peaking $F_q$")
    ax2.tick_params(axis="y", labelcolor=color2)
    ax2.grid(False)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")
    plt.title("Peak Pin Power and Peaking Factor Trajectory across Coupling", pad=12)
    fig.tight_layout()
    p3 = os.path.join(FIG_DIR, "peak_power_vs_iteration.png")
    plt.savefig(p3)
    plt.close()
    print(f"  [+] Saved: {p3}")

    # -------------------------------------------------------------------------
    # 4. power_convergence.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters[1:], df_hist["power_RMS_change"].values[1:] * 100.0, "o-", color="darkorange",
            linewidth=2.0, label="Power RMS Change [%]")
    ax.plot(iters[1:], df_hist["delta_peak_power"].values[1:] * 100.0, "s--", color="brown",
            linewidth=1.8, label="Peak Pin Power Change [%]")
    ax.axhline(2.5, color="gray", linestyle=":", label="Monte Carlo Statistical Noise Threshold (~2.5%)")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel("Relative Power Variation [%]")
    ax.set_xticks(iters[1:])
    ax.set_title("Power Distribution Convergence vs. Monte Carlo Tally Noise", pad=12)
    ax.legend(frameon=True, loc="upper right")
    fig.tight_layout()
    p4 = os.path.join(FIG_DIR, "power_convergence.png")
    plt.savefig(p4)
    plt.close()
    print(f"  [+] Saved: {p4}")

    # -------------------------------------------------------------------------
    # 5. fuel_temperature_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["Tfuel_max"], "o-", color="firebrick", linewidth=2.2, label="Max Fuel Centerline Temp [K]")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel("Peak Fuel Temperature [K]")
    ax.set_xticks(iters)
    ax.set_title("Fuel Hot Spot Temperature Stability vs. Iteration", pad=12)
    ax.legend(frameon=True, loc="center right")
    fig.tight_layout()
    p5 = os.path.join(FIG_DIR, "fuel_temperature_vs_iteration.png")
    plt.savefig(p5)
    plt.close()
    print(f"  [+] Saved: {p5}")

    # -------------------------------------------------------------------------
    # 6. clad_temperature_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["Tclad_max"], "s-", color="darkgoldenrod", linewidth=2.0, label="Max Clad Outer Temp [K]")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel("Peak Cladding Temperature [K]")
    ax.set_xticks(iters)
    ax.set_title("Cladding Outer Hot Spot Temperature Stability vs. Iteration", pad=12)
    ax.legend(frameon=True, loc="center right")
    fig.tight_layout()
    p6 = os.path.join(FIG_DIR, "clad_temperature_vs_iteration.png")
    plt.savefig(p6)
    plt.close()
    print(f"  [+] Saved: {p6}")

    # -------------------------------------------------------------------------
    # 7. coolant_temperature_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["Tcoolant_out"], "^-", color="teal", linewidth=2.0, label="Mixed Core Outlet Temp [K]")
    ax.plot(iters, df_hist["Tcoolant_max"], "d--", color="steelblue", linewidth=1.8, label="Hot Channel Exit Temp [K]")
    ax.axhline(1150.0, color="gray", linestyle=":", label="Nominal Target Outlet (1150.0 K)")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel("Coolant Temperature [K]")
    ax.set_xticks(iters)
    ax.set_title("Helium Coolant Outlet Temperatures vs. Coupling Iteration", pad=12)
    ax.legend(frameon=True, loc="center right")
    fig.tight_layout()
    p7 = os.path.join(FIG_DIR, "coolant_temperature_vs_iteration.png")
    plt.savefig(p7)
    plt.close()
    print(f"  [+] Saved: {p7}")

    # -------------------------------------------------------------------------
    # 8. density_feedback_vs_iteration.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["He_density_min"], "o-", color="darkcyan", linewidth=2.0, label=r"Min He Density (Hot Outlet) [$\mathrm{kg/m}^3$]")
    ax.plot(iters, df_hist["He_density_max"], "s--", color="midnightblue", linewidth=1.8, label=r"Max He Density (Inlet) [$\mathrm{kg/m}^3$]")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel(r"Helium Coolant Density [$\mathrm{kg/m}^3$]")
    ax.set_xticks(iters)
    ax.set_title("Axial Helium Density Feedback Range vs. Picard Iteration", pad=12)
    ax.legend(frameon=True, loc="center right")
    fig.tight_layout()
    p8 = os.path.join(FIG_DIR, "density_feedback_vs_iteration.png")
    plt.savefig(p8)
    plt.close()
    print(f"  [+] Saved: {p8}")

    # -------------------------------------------------------------------------
    # 9. axial_power_profile_final.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    z = df_axial["z_mid_cm"].values
    ax.plot(z, df_axial["axial_peaking_Fz"], "o-", color="navy", linewidth=2.2, label=r"Coupled $F_z(z)$")
    ax.axhline(1.0, color="gray", linestyle=":", label="Core Average ($F_z = 1.0$)")
    ax.set_xlabel("Axial Elevation z [cm] (Flow: z = 0 -> 42.45 cm)")
    ax.set_ylabel(r"Axial Peaking Factor $F_z(z)$")
    ax.set_title("Final Converged Core Axial Power Profile", pad=12)
    ax.set_xlim(0.0, 42.45)
    ax.legend(frameon=True, loc="upper right")
    fig.tight_layout()
    p9 = os.path.join(FIG_DIR, "axial_power_profile_final.png")
    plt.savefig(p9)
    plt.close()
    print(f"  [+] Saved: {p9}")

    # -------------------------------------------------------------------------
    # 10. radial_power_distribution_final.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 6.8))
    sc = ax.scatter(df_pins["x_cm"], df_pins["y_cm"], c=df_pins["power_W"], cmap="viridis",
                    s=140, edgecolors="k", linewidth=0.6)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pin Power [W]")
    
    p_peak_row = df_pins.loc[df_pins["power_W"].idxmax()]
    ax.plot(p_peak_row["x_cm"], p_peak_row["y_cm"], "r*", markersize=14,
            label=f"Peak Pin #{int(p_peak_row['pin_id'])} ({p_peak_row['power_W']:.1f} W)")
    
    ax.set_xlabel("X Coordinate [cm]")
    ax.set_ylabel("Y Coordinate [cm]")
    ax.set_title("Final Converged 2D Radial Pin Power Distribution", pad=12)
    ax.set_aspect("equal")
    ax.legend(loc="upper left")
    fig.tight_layout()
    p10 = os.path.join(FIG_DIR, "radial_power_distribution_final.png")
    plt.savefig(p10)
    plt.close()
    print(f"  [+] Saved: {p10}")

    # -------------------------------------------------------------------------
    # 11. fuel_temperature_axial_profile.png
    # -------------------------------------------------------------------------
    p_peak_id = int(p_peak_row["pin_id"])
    p_avg_id = int(df_pins.iloc[(df_pins["power_W"] - df_pins["power_W"].mean()).abs().argsort()[:1]]["pin_id"].values[0])
    p_edge_id = int(df_pins.loc[df_pins["ring"] == 8]["power_W"].idxmin())
    
    df_peak_th = df_therm[df_therm["pin_id"] == p_peak_id]
    df_avg_th = df_therm[df_therm["pin_id"] == p_avg_id]
    df_edge_th = df_therm[df_therm["pin_id"] == p_edge_id]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df_peak_th["z_mid_cm"], df_peak_th["T_fuel_center_K"], "r-", linewidth=2.2, label=f"Peak Pin #{p_peak_id} Centerline")
    ax.plot(df_peak_th["z_mid_cm"], df_peak_th["T_fuel_vol_K"], "r--", linewidth=1.8, label=f"Peak Pin #{p_peak_id} Vol-Avg")
    ax.plot(df_avg_th["z_mid_cm"], df_avg_th["T_fuel_center_K"], "g-", linewidth=2.0, label=f"Average Pin #{p_avg_id} Centerline")
    ax.plot(df_avg_th["z_mid_cm"], df_avg_th["T_fuel_vol_K"], "g--", linewidth=1.6, label=f"Average Pin #{p_avg_id} Vol-Avg")
    ax.plot(df_edge_th["z_mid_cm"], df_edge_th["T_fuel_center_K"], "b-", linewidth=1.8, label=f"Edge Pin #{p_edge_id} Centerline")
    
    ax.set_xlabel("Axial Elevation $z$ [cm]")
    ax.set_ylabel("Fuel Temperature [K]")
    ax.set_title("Axial Fuel Temperature Profiles (Centerline and Volume-Average)", pad=12)
    ax.set_xlim(0.0, 42.45)
    ax.legend(frameon=True, loc="upper left")
    fig.tight_layout()
    p11 = os.path.join(FIG_DIR, "fuel_temperature_axial_profile.png")
    plt.savefig(p11)
    plt.close()
    print(f"  [+] Saved: {p11}")

    # -------------------------------------------------------------------------
    # 12. coolant_temperature_axial_profile.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df_peak_th["z_mid_cm"], df_peak_th["T_cool_bulk_K"], "r-", linewidth=2.0, label=f"Hot Channel #{p_peak_id} Bulk He")
    ax.plot(df_avg_th["z_mid_cm"], df_avg_th["T_cool_bulk_K"], "g-", linewidth=2.0, label=f"Average Channel #{p_avg_id} Bulk He")
    ax.plot(df_edge_th["z_mid_cm"], df_edge_th["T_cool_bulk_K"], "b-", linewidth=1.8, label=f"Edge Channel #{p_edge_id} Bulk He")
    ax.plot(df_peak_th["z_mid_cm"], df_peak_th["T_clad_out_K"], "m--", linewidth=1.8, label=f"Hot Pin #{p_peak_id} Clad Surface")
    ax.plot(df_avg_th["z_mid_cm"], df_avg_th["T_clad_out_K"], "c--", linewidth=1.6, label=f"Average Pin #{p_avg_id} Clad Surface")
    
    ax.set_xlabel("Axial Elevation z [cm] (Flow from Inlet z=0 to Outlet z=42.45 cm)")
    ax.set_ylabel("Coolant & Cladding Temperature [K]")
    ax.set_title("Axial Coolant Bulk and Outer Cladding Temperature Profiles", pad=12)
    ax.set_xlim(0.0, 42.45)
    ax.legend(frameon=True, loc="upper left")
    fig.tight_layout()
    p12 = os.path.join(FIG_DIR, "coolant_temperature_axial_profile.png")
    plt.savefig(p12)
    plt.close()
    print(f"  [+] Saved: {p12}")

    # -------------------------------------------------------------------------
    # 13. pressure_drop.png
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_hist["pressure_drop"] / 1000.0, "o-", color="darkslateblue", linewidth=2.0, label=r"Total Channel $\Delta P$ [kPa]")
    ax.set_xlabel("Picard Coupling Iteration")
    ax.set_ylabel(r"Total Coolant Pressure Drop $\Delta P$ [kPa]")
    ax.set_xticks(iters)
    ax.set_title("Core Channel Pressure Drop vs. Coupling Iteration", pad=12)
    ax.legend(frameon=True, loc="center right")
    fig.tight_layout()
    p13 = os.path.join(FIG_DIR, "pressure_drop.png")
    plt.savefig(p13)
    plt.close()
    print(f"  [+] Saved: {p13}")

    print("All 13 publication-quality figures generated successfully!")


if __name__ == "__main__":
    generate_all_13_plots()
