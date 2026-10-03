"""
PLOTTING SUITE FOR V3 ADVANCED COUPLED CAMPAIGN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/generate_v3_plots.py

Generates 9 high-resolution (300 DPI) publication-quality figures:
1. v3_01_convergence_k.png: k_eff vs Picard iteration with MC error bars and convergence band
2. v3_02_convergence_temperature.png: Max fuel, clad, and coolant temperatures vs iteration
3. v3_03_convergence_power.png: Peak pin power and relative RMS power change vs iteration
4. v3_04_axial_power.png: Axial power profile across iterations comparing uncoupled vs converged
5. v3_05_radial_power.png: 2D radial pin power distribution map across 9 hexagonal rings
6. v3_06_fuel_clad_coolant_temp.png: Complete axial temperature profiles (T_fuel, T_clad, T_cool)
7. v3_07_pressure.png: Axial coolant pressure distribution and component pressure drops
8. v3_08_density.png: Helium coolant density rho(z) and flow velocity u_z(z) profiles
9. v3_09_mesh_sensitivity.png: Axial discretization comparison (10 slices vs 20 slices)
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.collections import PatchCollection

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_V3_DIR = os.path.join(WORKSPACE_ROOT, "results", "v3")
FIG_DIR = os.path.join(RESULTS_V3_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Set style
plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "grid.alpha": 0.5,
    "grid.linestyle": "--"
})


def plot_all_v3_figures():
    print("Generating all 9 V3 publication-grade figures...")
    
    # Load dataframes
    df_iters = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_iterations.csv"))
    df_therm = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_thermal_field.csv"))
    df_cool = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_coolant_field.csv"))
    df_pins = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_pin_power.csv"))
    df_mesh = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_mesh_sensitivity.csv"))
    df_comp = pd.read_csv(os.path.join(RESULTS_V3_DIR, "v3_comparison.csv"))
    
    iters = df_iters["iteration"].values
    
    # -------------------------------------------------------------------------
    # 1. Convergence k
    # -------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.errorbar(iters, df_iters["k_eff"], yerr=df_iters["sigma_k"], fmt="o-", color="#1f77b4",
                 ecolor="#aec7e8", elinewidth=2, capsize=4, label=r"$k_{\mathrm{eff}} \pm \sigma$")
    ax1.set_xlabel("Picard Iteration")
    ax1.set_ylabel(r"Effective Multiplication Factor $k_{\mathrm{eff}}$", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.grid(True)
    
    ax2 = ax1.twinx()
    ax2.plot(iters[1:], np.abs(df_iters["delta_k"].values[1:]), "s--", color="#d62728", label=r"$|\Delta k|$")
    ax2.plot(iters[1:], df_iters["sigma_k"].values[1:] * 1.5, "k:", label=r"Convergence Threshold $1.5\sigma_k$")
    ax2.set_ylabel(r"$|\Delta k|$ between iterations", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    ax2.set_yscale("log")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")
    plt.title("V3 Multiphysics Coupling: Eigenvalue Convergence vs Picard Iteration")
    plt.tight_layout()
    p1 = os.path.join(FIG_DIR, "v3_01_convergence_k.png")
    plt.savefig(p1, dpi=300)
    plt.close()
    print(f"  [1/9] Saved {p1}")

    # -------------------------------------------------------------------------
    # 2. Convergence Temperature
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(iters, df_iters["T_fuel_max_K"], "o-", color="#d62728", label=r"Peak Fuel Centerline $T_{f,\max}$")
    ax.plot(iters, df_iters["cardinal_T_fuel_max_K"], "^--", color="#ff7f0e", label=r"Cardinal FE Peak Fuel $T_{f,\max}$")
    ax.plot(iters, df_iters["T_clad_max_K"], "s-", color="#2ca02c", label=r"Peak Cladding Outer $T_{c,\max}$")
    ax.plot(iters, df_iters["T_coolant_max_K"], "d-", color="#1f77b4", label=r"Peak Coolant Bulk $T_{cool,\max}$")
    ax.plot(iters, df_iters["T_coolant_out_K"], "x-", color="#9467bd", label=r"Mixed Coolant Outlet $\bar{T}_{out}$")
    ax.set_xlabel("Picard Iteration")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("V3 Multiphysics Coupling: Temperature Convergence vs Picard Iteration")
    ax.grid(True)
    ax.legend(loc="center right")
    plt.tight_layout()
    p2 = os.path.join(FIG_DIR, "v3_02_convergence_temperature.png")
    plt.savefig(p2, dpi=300)
    plt.close()
    print(f"  [2/9] Saved {p2}")

    # -------------------------------------------------------------------------
    # 3. Convergence Power Distribution
    # -------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(iters, df_iters["peak_pin_power_W"], "o-", color="#1f77b4", label=r"Peak Pin Power $P_{\mathrm{peak}}$")
    ax1.set_xlabel("Picard Iteration")
    ax1.set_ylabel("Peak Pin Power [W]", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.grid(True)
    
    ax2 = ax1.twinx()
    ax2.plot(iters[1:], df_iters["power_RMS_change"].values[1:] * 100.0, "r^--", label=r"Pin Power RMS Change (%)")
    ax2.axhline(1.5, color="k", linestyle=":", label="Convergence Target (1.5%)")
    ax2.set_ylabel("Relative RMS Change [%]", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")
    plt.title("V3 Multiphysics Coupling: Power Distribution Convergence")
    plt.tight_layout()
    p3 = os.path.join(FIG_DIR, "v3_03_convergence_power.png")
    plt.savefig(p3, dpi=300)
    plt.close()
    print(f"  [3/9] Saved {p3}")

    # -------------------------------------------------------------------------
    # 4. Axial Power Distribution
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    z_cm = df_therm["z_cm"].values
    q_vol_peak = df_therm["q_vol_peak_W_m3"].values / 1e6
    q_vol_avg = df_therm["q_vol_avg_W_m3"].values / 1e6
    
    ax.plot(z_cm, q_vol_peak, "r-o", label=f"Peak Pin #{int(df_iters['peak_pin_idx'].iloc[-1])} Volumetric Power")
    ax.plot(z_cm, q_vol_avg, "b--s", label="Core-Average Axial Profile")
    ax.axhline(np.mean(q_vol_avg), color="gray", linestyle=":", label="Core Average Reference")
    ax.set_xlabel("Axial Elevation z [cm]")
    ax.set_ylabel(r"Volumetric Heat Generation $q'''$ [MW/m$^3$]")
    ax.set_title("V3 Final Converged Axial Volumetric Power Profile")
    ax.grid(True)
    ax.legend(loc="upper right")
    plt.tight_layout()
    p4 = os.path.join(FIG_DIR, "v3_04_axial_power.png")
    plt.savefig(p4, dpi=300)
    plt.close()
    print(f"  [4/9] Saved {p4}")

    # -------------------------------------------------------------------------
    # 5. Radial Power Distribution Map
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(df_pins["x_cm"], df_pins["y_cm"], c=df_pins["peaking_Fr"],
                    s=120, cmap="inferno", edgecolors="black", linewidths=0.8)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(r"Radial Peaking Factor $F_r = P_{\mathrm{pin}} / \bar{P}_{\mathrm{pin}}$")
    ax.set_xlabel("Core X [cm]")
    ax.set_ylabel("Core Y [cm]")
    ax.set_title("V3 Opt-17 Converged Radial Pin Power Peaking Map (217 Pins)")
    ax.set_aspect("equal")
    ax.grid(True)
    plt.tight_layout()
    p5 = os.path.join(FIG_DIR, "v3_05_radial_power.png")
    plt.savefig(p5, dpi=300)
    plt.close()
    print(f"  [5/9] Saved {p5}")

    # -------------------------------------------------------------------------
    # 6. Fuel, Clad, Coolant Temperatures
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(z_cm, df_therm["T_fuel_center_peak_K"], "r-o", label=r"Peak Pin Fuel Centerline $T_{f,0}(z)$")
    ax.plot(z_cm, df_therm["T_fuel_vol_peak_K"], "m--^", label=r"Peak Pin Fuel Vol-Avg $T_{f,\mathrm{vol}}(z)$")
    ax.plot(z_cm, df_therm["T_clad_out_peak_K"], "g-s", label=r"Peak Pin Clad Outer Surface $T_{c,o}(z)$")
    ax.plot(z_cm, df_therm["T_cool_bulk_peak_K"], "b-d", label=r"Peak Channel Coolant Bulk $T_b(z)$")
    ax.plot(z_cm, df_therm["T_cool_bulk_core_avg_K"], "c--", label=r"Core-Average Coolant Bulk $\bar{T}_b(z)$")
    ax.axhline(1150.0, color="gray", linestyle=":", label="Nominal Mixed Exit Target (1150 K)")
    ax.set_xlabel("Axial Elevation z [cm]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("V3 Coupled Reactor: Axial Temperature Profiles Along Active Core")
    ax.grid(True)
    ax.legend(loc="upper left")
    plt.tight_layout()
    p6 = os.path.join(FIG_DIR, "v3_06_fuel_clad_coolant_temperature.png")
    plt.savefig(p6, dpi=300)
    plt.close()
    print(f"  [6/9] Saved {p6}")

    # -------------------------------------------------------------------------
    # 7. Pressure Distribution & Drop
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    p_tot_kPa = df_cool["pressure_drop_accum_Pa"].values / 1000.0
    p_abs_MPa = (2.5e6 + df_cool["pressure_drop_accum_Pa"].values[::-1]) / 1e6
    ax.plot(z_cm, p_tot_kPa, "r-o", label=r"Cumulative Channel Pressure Drop $\Delta P(z)$ [kPa]")
    ax.set_xlabel("Axial Elevation z [cm]")
    ax.set_ylabel("Pressure Drop [kPa]", color="r")
    ax.tick_params(axis="y", labelcolor="r")
    ax.grid(True)
    
    ax2 = ax.twinx()
    ax2.plot(z_cm, p_abs_MPa, "b--s", label="Absolute Pressure [MPa]")
    ax2.set_ylabel("Absolute System Pressure [MPa]", color="b")
    ax2.tick_params(axis="y", labelcolor="b")
    
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="center left")
    plt.title("V3 Helium Coolant Flow: Axial Pressure Distribution")
    plt.tight_layout()
    p7 = os.path.join(FIG_DIR, "v3_07_pressure.png")
    plt.savefig(p7, dpi=300)
    plt.close()
    print(f"  [7/9] Saved {p7}")

    # -------------------------------------------------------------------------
    # 8. Coolant Density & Velocity
    # -------------------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(z_cm, df_cool["density_He_kg_m3"], "b-o", label=r"Helium Coolant Density $\rho_{\mathrm{He}}(z)$")
    ax1.set_xlabel("Axial Elevation z [cm]")
    ax1.set_ylabel(r"Helium Density [kg/m$^3$]", color="b")
    ax1.tick_params(axis="y", labelcolor="b")
    ax1.grid(True)
    
    ax2 = ax1.twinx()
    ax2.plot(z_cm, df_cool["velocity_m_s"], "g--s", label=r"Coolant Velocity $u_z(z)$")
    ax2.set_ylabel("Coolant Velocity [m/s]", color="g")
    ax2.tick_params(axis="y", labelcolor="g")
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")
    plt.title("V3 Annular Flow Channel: Helium Coolant Density & Velocity Profiles")
    plt.tight_layout()
    p8 = os.path.join(FIG_DIR, "v3_08_density.png")
    plt.savefig(p8, dpi=300)
    plt.close()
    print(f"  [8/9] Saved {p8}")

    # -------------------------------------------------------------------------
    # 9. Mesh Discretization Sensitivity
    # -------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    
    # Bar chart of key metrics
    metrics_to_plot = ["k_eff", "F_z", "F_q", "F_r"]
    df_sub = df_mesh[df_mesh["metric"].isin(metrics_to_plot)].copy()
    x = np.arange(len(df_sub))
    w = 0.35
    ax1.bar(x - w/2, df_sub["case_10_slices"], w, label="10 Slices (Baseline)", color="#1f77b4")
    ax1.bar(x + w/2, df_sub["case_20_slices"], w, label="20 Slices (Refined)", color="#ff7f0e")
    ax1.set_xticks(x)
    ax1.set_xticklabels(df_sub["metric"])
    ax1.set_title("Neutronic Peaking & Eigenvalue Stability")
    ax1.legend()
    ax1.grid(True)
    
    # Temperature comparison bar
    temp_metrics = ["T_fuel,max (K)", "T_clad,max (K)"]
    df_temp = df_mesh[df_mesh["metric"].isin(temp_metrics)].copy()
    xt = np.arange(len(df_temp))
    ax2.bar(xt - w/2, df_temp["case_10_slices"], w, label="10 Slices (Baseline)", color="#1f77b4")
    ax2.bar(xt + w/2, df_temp["case_20_slices"], w, label="20 Slices (Refined)", color="#ff7f0e")
    ax2.set_xticks(xt)
    ax2.set_xticklabels(df_temp["metric"])
    ax2.set_ylabel("Temperature [K]")
    ax2.set_ylim(1000, 1400)
    ax2.set_title("Thermal Maxima Sensitivity")
    ax2.legend()
    ax2.grid(True)
    
    plt.suptitle("V3 Spatial Discretization Sensitivity: 10 vs 20 Axial Slices (2,170 vs 4,340 Cells)")
    plt.tight_layout()
    p9 = os.path.join(FIG_DIR, "v3_09_mesh_sensitivity.png")
    plt.savefig(p9, dpi=300)
    plt.close()
    print(f"  [9/9] Saved {p9}")
    print("All 9 figures generated successfully at 300 DPI!")


if __name__ == "__main__":
    plot_all_v3_figures()
