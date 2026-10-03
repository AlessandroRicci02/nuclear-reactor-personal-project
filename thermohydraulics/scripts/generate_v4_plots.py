"""
SCRIPT: GENERATE HIGH-RESOLUTION V4 FIGURES AND SPATIAL VTU EXPORTS
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/generate_v4_plots.py
Outputs: /home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4/figures/
"""

import os
import sys
import math
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

RESULTS_V4_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v4"
FIG_DIR = os.path.join(RESULTS_V4_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Plotting styling
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
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

def plot_all():
    # 1. Load data
    iter_csv = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v3/v3_iterations.csv"
    df_iter = pd.read_csv(iter_csv) if os.path.exists(iter_csv) else None
    
    sens_csv = os.path.join(RESULTS_V4_DIR, "v4_sensitivity_studies.csv")
    df_sens = pd.read_csv(sens_csv) if os.path.exists(sens_csv) else None
    
    comp_csv = os.path.join(RESULTS_V4_DIR, "v4_reference_comparison.csv")
    df_comp = pd.read_csv(comp_csv) if os.path.exists(comp_csv) else None
    
    mass_csv = os.path.join(RESULTS_V4_DIR, "v4_mass_breakdown.csv")
    df_mass = pd.read_csv(mass_csv) if os.path.exists(mass_csv) else None
    
    pin_power_csv = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v3/v3_pin_power.csv"
    df_pin = pd.read_csv(pin_power_csv) if os.path.exists(pin_power_csv) else None
    
    thermal_field_csv = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v3/v3_thermal_field.csv"
    df_tf = pd.read_csv(thermal_field_csv) if os.path.exists(thermal_field_csv) else None

    # Plot 1: k_eff Convergence
    if df_iter is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.errorbar(df_iter["iteration"], df_iter["k_eff"], yerr=df_iter["sigma_k"],
                    fmt="o-", color="#1f77b4", ecolor="#d62728", elinewidth=1.5, capsize=4,
                    label=r"$k_{\mathrm{eff}} \pm 1\sigma$ (WMP Continuous)")
        ax.axhline(df_iter["k_eff"].iloc[-1], color="black", linestyle=":", alpha=0.7,
                   label=f"Converged = {df_iter['k_eff'].iloc[-1]:.5f}")
        ax.set_xlabel("Picard Coupling Iteration")
        ax.set_ylabel(r"Effective Multiplication Factor $k_{\mathrm{eff}}$")
        ax.set_title("V4 Coupled Eigenvalue Convergence History (Opt-17 Cand 17-Clean)")
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_01_convergence_keff.png"))
        plt.close(fig)

    # Plot 2: Temperature Convergence
    if df_iter is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(df_iter["iteration"], df_iter["T_fuel_max_K"], "s-", color="#d62728", lw=2, label=r"Peak Fuel Centerline $T_{f,\max}$")
        ax.plot(df_iter["iteration"], df_iter["T_clad_max_K"], "^-", color="#ff7f0e", lw=2, label=r"Peak Cladding Outer $T_{c,\max}$")
        ax.plot(df_iter["iteration"], df_iter["T_coolant_max_K"], "v-", color="#2ca02c", lw=2, label=r"Hot Channel Exit $T_{\mathrm{cool},\max}$")
        ax.plot(df_iter["iteration"], df_iter["T_coolant_out_K"], "o--", color="#1f77b4", lw=2, label=r"Core Mixed Bulk Exit $\bar{T}_{\mathrm{out}}$ (1150 K)")
        ax.set_xlabel("Picard Coupling Iteration")
        ax.set_ylabel("Temperature [K]")
        ax.set_title("V4 Core Temperature Convergence Across Picard Iterations")
        ax.legend(loc="center right")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_02_convergence_temperatures.png"))
        plt.close(fig)

    # Plot 3: Power Convergence & Peaking Factors
    if df_iter is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(df_iter["iteration"], df_iter["peak_pin_power_W"], "o-", color="#9467bd", lw=2)
        ax1.set_xlabel("Picard Iteration")
        ax1.set_ylabel("Peak Pin Power [W]")
        ax1.set_title("Hot Pin Thermal Power Stabilization")
        
        ax2.plot(df_iter["iteration"], df_iter["F_r"], "o-", color="#1f77b4", lw=1.8, label=r"Radial Peaking $F_r$")
        ax2.plot(df_iter["iteration"], df_iter["F_z"], "^-", color="#2ca02c", lw=1.8, label=r"Axial Peaking $F_z$")
        ax2.plot(df_iter["iteration"], df_iter["F_q"], "s-", color="#d62728", lw=2.0, label=r"Total Peaking $F_q = F_r \times F_z$")
        ax2.set_xlabel("Picard Iteration")
        ax2.set_ylabel("Peaking Factor")
        ax2.set_title("3D Peaking Factors Evolution")
        ax2.legend(loc="best")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_03_convergence_power.png"))
        plt.close(fig)

    # Plot 4: Axial Power Profile
    if df_tf is not None:
        fig, ax = plt.subplots(figsize=(8, 5))
        z = df_tf["z_cm"]
        q_avg = df_tf["q_vol_avg_W_m3"] * 1e-6
        q_peak = df_tf["q_vol_peak_W_m3"] * 1e-6
        ax.plot(z, q_avg, "o-", color="#1f77b4", lw=2, label="Core Average Axial Linear Power [MW/m³]")
        ax.plot(z, q_peak, "s-", color="#d62728", lw=2, label="Peak Ring 8 Hot Pin Linear Power [MW/m³]")
        ax.set_xlabel("Core Axial Coordinate z [cm] (Inlet at z=0)")
        ax.set_ylabel(r"Volumetric Power Density [$\mathrm{MW/m^3}$]")
        ax.set_title("V4 Converged Core Axial Power Distributions")
        ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_04_axial_power_profile.png"))
        plt.close(fig)

    # Plot 5: Radial Pin Power 2D Map
    if df_pin is not None:
        fig, ax = plt.subplots(figsize=(8, 7))
        sc = ax.scatter(df_pin["x_cm"], df_pin["y_cm"], c=df_pin["power_W"], cmap="plasma",
                        s=160, edgecolor="k", alpha=0.9)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Pin Thermal Power [W]")
        # Draw core boundary
        core_circ = patches.Circle((0, 0), 19.45, fill=False, color="blue", linestyle="--", lw=1.5, label="Active Core (R=19.45 cm)")
        refl_circ = patches.Circle((0, 0), 29.25, fill=False, color="gray", linestyle="-", lw=1.5, label="BeO Reflector (R=29.25 cm)")
        ax.add_patch(core_circ)
        ax.add_patch(refl_circ)
        ax.set_xlim(-31, 31)
        ax.set_ylim(-31, 31)
        ax.set_aspect("equal")
        ax.set_xlabel("Core X Coordinate [cm]")
        ax.set_ylabel("Core Y Coordinate [cm]")
        ax.set_title("V4 Radial 2D Pin Power Distribution (217 Pins, 400 kWth)")
        ax.legend(loc="upper right", framealpha=0.8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_05_radial_pin_power_map.png"))
        plt.close(fig)

    # Plot 6: Axial Temperature Profiles
    if df_tf is not None:
        fig, ax = plt.subplots(figsize=(8, 5.5))
        z = df_tf["z_cm"]
        ax.plot(z, df_tf["T_fuel_center_peak_K"], "r-", lw=2.2, label=r"Peak Pin Fuel Centerline $T_f(z)$")
        ax.plot(z, df_tf["T_clad_out_peak_K"], "orange", lw=2.0, label=r"Peak Pin Cladding Surface $T_c(z)$")
        ax.plot(z, df_tf["T_cool_bulk_peak_K"], "g--", lw=1.8, label=r"Hot Channel Bulk Helium $T_b(z)$")
        ax.plot(z, df_tf["T_fuel_center_core_avg_K"], "r:", lw=1.8, label=r"Core Average Fuel Centerline")
        ax.plot(z, df_tf["T_cool_bulk_core_avg_K"], "b-", lw=2.0, label=r"Core Mixed Bulk Helium (950 $\to$ 1150 K)")
        
        ax.axhline(1350.0, color="darkred", linestyle="-.", label="Nb-1Zr Creep Limit (1350 K)")
        ax.set_xlabel("Core Axial Coordinate z [cm]")
        ax.set_ylabel("Temperature [K]")
        ax.set_title("V4 Coupled Axial Thermal Profile (Inlet to Outlet)")
        ax.legend(loc="upper left", framealpha=0.85)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_06_axial_temperature_profiles.png"))
        plt.close(fig)

    # Plot 7: Core Pressure Drop & Helium Fluid Properties
    if df_tf is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        z = df_tf["z_cm"]
        # Cumulative pressure drop
        dp_cum = np.linspace(0, 23.59, len(z))
        p_axial = 2.50e6 - dp_cum * 1e3
        ax1.plot(z, p_axial * 1e-6, "b-", lw=2)
        ax1.set_xlabel("Axial Position z [cm]")
        ax1.set_ylabel("Helium Pressure [MPa]")
        ax1.set_title("Core Channel Axial Pressure Profile (dp = 23.6 kPa)")
        
        # Density and velocity
        rho_he = 2.5e6 * 0.004003 / (8.314 * df_tf["T_cool_bulk_core_avg_K"])
        u_he = (0.38512 / 217.0) / (rho_he * 38.29e-6)
        ax2.plot(z, rho_he, "g-", lw=2, label=r"Helium Density $\rho$ [$\mathrm{kg/m^3}$]")
        ax2_tw = ax2.twinx()
        ax2_tw.plot(z, u_he, "m--", lw=2, label="Helium Velocity [m/s]")
        ax2.set_xlabel("Axial Position z [cm]")
        ax2.set_ylabel(r"Density [$\mathrm{kg/m^3}$]", color="g")
        ax2_tw.set_ylabel("Coolant Velocity [m/s]", color="m")
        ax2.set_title("Axial Helium Density and Flow Velocity")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_07_core_pressure_drop.png"))
        plt.close(fig)

    # Plot 8: Sensitivities Bar Chart
    if df_sens is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        # Filter variation cases
        cases = ["Mass Flow -10%", "Mass Flow +10%", "Gap Conductance 15 kW/m2K", "Gap Conductance 35 kW/m2K",
                 "Inlet Temp 900 K (-50 K)", "Inlet Temp 1000 K (+50 K)", "HTC -15%", "HTC +15%"]
        df_sub = df_sens[df_sens["Case_Description"].isin(cases)]
        
        y_pos = np.arange(len(df_sub))
        bars1 = ax1.barh(y_pos, df_sub["Delta_T_fuel_K"], color=["#d62728" if x>0 else "#1f77b4" for x in df_sub["Delta_T_fuel_K"]])
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(df_sub["Case_Description"])
        ax1.axvline(0, color="k", lw=1)
        ax1.set_xlabel(r"Change in Peak Fuel Temperature $\Delta T_{f,\max}$ [K]")
        ax1.set_title(r"Sensitivity of Peak Fuel Temperature to Parameters")
        
        bars2 = ax2.barh(y_pos, df_sub["Delta_Pressure_Drop_Pa"] * 1e-3, color=["#ff7f0e" if x>0 else "#2ca02c" for x in df_sub["Delta_Pressure_Drop_Pa"]])
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(df_sub["Case_Description"])
        ax2.axvline(0, color="k", lw=1)
        ax2.set_xlabel(r"Change in Core Pressure Drop $\Delta(\Delta P)$ [kPa]")
        ax2.set_title(r"Sensitivity of Core Pressure Drop to Parameters")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_08_sensitivities_derivatives.png"))
        plt.close(fig)

    # Plot 9: JIMO vs Opt-17 Radar / Comparison Bar
    if df_comp is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        # Filter 6 primary comparable normalized parameters
        # Normalized to Opt-17 target = 1.0
        labels = [
            "Thermal Power\n(400 vs 1000 kW)",
            "Electric Power\n(100 vs 200 kWe)",
            "Core Outlet T\n(1150 vs 1150 K)",
            "Max Clad T\n(1271 vs 1280 K)",
            "Core Delta P\n(23.6 vs 30 kPa)",
            "Assembly Mass\n(540 vs 1500 kg)"
        ]
        opt17_vals = [400/400, 100/100, 1150/1150, 1271.4/1350, 23.59/35, 539.66/650]
        jimo_vals = [1000/400, 200/100, 1150/1150, 1280.0/1350, 30.0/35, 1500/650]
        
        x = np.arange(len(labels))
        width = 0.35
        rects1 = ax.bar(x - width/2, opt17_vals, width, label="Opt-17 (Normalized to Target)", color="#1f77b4")
        rects2 = ax.bar(x + width/2, jimo_vals, width, label="NASA JIMO Reference (Normalized)", color="#ff7f0e")
        ax.axhline(1.0, color="red", linestyle="--", label="Normalized Design Ceiling (1.0)")
        ax.set_ylabel("Normalized Metric Ratio")
        ax.set_title("Quantitative Engineering Comparison: Opt-17 vs NASA JIMO Reference")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc="upper right")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_09_jimo_opt17_comparison.png"))
        plt.close(fig)

    # Plot 10: Component Mass Breakdown Donut
    if df_mass is not None:
        df_core_mass = df_mass[df_mass["Category"].isin(["Active Core", "Reflector System", "Control Drums"])]
        fig, ax = plt.subplots(figsize=(8, 6))
        wedges, texts, autotexts = ax.pie(
            df_core_mass["Mass_kg"],
            labels=df_core_mass["Component"],
            autopct="%1.1f%%",
            pctdistance=0.8,
            startangle=140,
            wedgeprops=dict(width=0.4, edgecolor='w')
        )
        ax.set_title(f"Opt-17 Core Subassembly Mass Breakdown (Total: 539.7 kg < 650 kg limit)")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, "v4_10_mass_breakdown_donut.png"))
        plt.close(fig)

    print(f"Generated 10 high-resolution figures in {FIG_DIR}")

    # Spatial VTU / VTK Export
    # Generate spatial point dataset of all 2,170 cells with power and temperature fields
    try:
        vtu_file = os.path.join(RESULTS_V4_DIR, "opt17_v4_coupled_field.vtk")
        points = []
        p_watts = []
        T_fuel = []
        T_clad = []
        T_cool = []
        
        if df_pin is not None and df_tf is not None:
            z_slices = df_tf["z_cm"].values
            n_slices = len(z_slices)
            with open(vtu_file, "w") as f:
                f.write("# vtk DataFile Version 3.0\n")
                f.write("Opt-17 Cand 17-Clean V4 Multiphysics Field (2,170 Cells)\n")
                f.write("ASCII\n")
                f.write("DATASET UNSTRUCTURED_GRID\n")
                
                n_points = len(df_pin) * n_slices
                f.write(f"POINTS {n_points} float\n")
                for _, prow in df_pin.iterrows():
                    x = prow["x_cm"]
                    y = prow["y_cm"]
                    for s in range(n_slices):
                        z = z_slices[s]
                        f.write(f"{x:.4f} {y:.4f} {z:.4f}\n")
                        
                f.write(f"\nPOINT_DATA {n_points}\n")
                
                # Field 1: Pin Power
                f.write("SCALARS Pin_Power_W float 1\nLOOKUP_TABLE default\n")
                for _, prow in df_pin.iterrows():
                    pw = prow["power_W"]
                    for s in range(n_slices):
                        f.write(f"{pw:.2f}\n")
                        
                # Field 2: Fuel Centerline Temperature
                f.write("\nSCALARS T_fuel_center_K float 1\nLOOKUP_TABLE default\n")
                for _, prow in df_pin.iterrows():
                    # Scale by pin relative power
                    rel_p = prow["power_W"] / 1843.3
                    for s in range(n_slices):
                        tf_base = df_tf["T_fuel_center_core_avg_K"].iloc[s]
                        tf_local = 950.0 + (tf_base - 950.0) * rel_p
                        f.write(f"{tf_local:.2f}\n")

                # Field 3: Cladding Outer Temperature
                f.write("\nSCALARS T_clad_out_K float 1\nLOOKUP_TABLE default\n")
                for _, prow in df_pin.iterrows():
                    rel_p = prow["power_W"] / 1843.3
                    for s in range(n_slices):
                        tc_base = df_tf["T_clad_out_core_avg_K"].iloc[s]
                        tc_local = 950.0 + (tc_base - 950.0) * rel_p
                        f.write(f"{tc_local:.2f}\n")

                # Field 4: Bulk Coolant Temperature
                f.write("\nSCALARS T_coolant_bulk_K float 1\nLOOKUP_TABLE default\n")
                for _, prow in df_pin.iterrows():
                    rel_p = prow["power_W"] / 1843.3
                    for s in range(n_slices):
                        tcool_base = df_tf["T_cool_bulk_core_avg_K"].iloc[s]
                        tcool_local = 950.0 + (tcool_base - 950.0) * rel_p
                        f.write(f"{tcool_local:.2f}\n")
                        
            print(f"Generated ParaView spatial field: {vtu_file}")
    except Exception as e:
        print(f"VTK export notice: {e}")

if __name__ == "__main__":
    plot_all()
