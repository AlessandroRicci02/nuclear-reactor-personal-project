"""
Generate Additional High-Resolution Publication Figures for Opt-5
=================================================================
Creates:
1. RZ Vertical Geometry Cut (opt5_geometry_rz.png)
2. Fuel Cell Zoom (opt5_fuel_cell_zoom.png)
3. Temperature Coefficients Bar Chart (opt5_temperature_coefficients.png)
4. Coolant Void Reactivity with Uncertainty (opt5_void_reactivity.png)
5. Fissile U-235 & Pu-239 Breeding Trajectory (opt5_u235_and_plutonium_inventory.png)
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv

sys.path.insert(0, os.path.abspath("src"))
from model import ReactorParameters

def generate_additional_plots(output_dir="results/opt5/figures"):
    os.makedirs(output_dir, exist_ok=True)
    p = ReactorParameters(
        fuel_radius=0.58,
        cladding_thickness=0.05,
        coolant_channel_thickness=0.05,
        moderator_fuel_pitch=1.635,
        core_radius=16.00,
        active_core_height=43.50,
        radial_reflector_thickness=9.60,
        axial_reflector_thickness=7.60,
        number_of_control_drums=8,
        control_drum_radius=4.30,
        absorber_thickness=1.45,
        absorber_sector_angular_width=140.0,
        drum_clearance=0.20,
        b10_enrichment=0.90
    )

    # 1. RZ Vertical Cut Plot
    fig, ax = plt.subplots(figsize=(7, 9), dpi=300)
    r_core = p.core_radius
    r_refl = p.radial_reflector_outer_radius
    h_core_half = p.active_core_height / 2.0
    h_tot_half = p.total_reactor_height / 2.0

    # Radial reflector
    refl_r = patches.Rectangle((r_core, -h_core_half), r_refl - r_core, p.active_core_height,
                               linewidth=1.2, edgecolor='#37474F', facecolor='#CFD8DC', label='BeO Radial Reflector')
    refl_l = patches.Rectangle((-r_refl, -h_core_half), r_refl - r_core, p.active_core_height,
                               linewidth=1.2, edgecolor='#37474F', facecolor='#CFD8DC')
    # Axial reflectors
    refl_top = patches.Rectangle((-r_refl, h_core_half), 2 * r_refl, p.axial_reflector_thickness,
                                 linewidth=1.2, edgecolor='#37474F', facecolor='#B0BEC5', label='BeO Axial Reflectors (Top/Bot)')
    refl_bot = patches.Rectangle((-r_refl, -h_tot_half), 2 * r_refl, p.axial_reflector_thickness,
                                 linewidth=1.2, edgecolor='#37474F', facecolor='#B0BEC5')
    # Core
    core_rect = patches.Rectangle((-r_core, -h_core_half), 2 * r_core, p.active_core_height,
                                  linewidth=1.5, edgecolor='#1B5E20', facecolor='#A5D6A7', label='Active Core (UN + YH2)')

    ax.add_patch(refl_r)
    ax.add_patch(refl_l)
    ax.add_patch(refl_top)
    ax.add_patch(refl_bot)
    ax.add_patch(core_rect)

    # Annotations
    ax.axhline(0, color='black', linestyle=':', lw=1.0)
    ax.axvline(0, color='black', linestyle=':', lw=1.0)
    ax.text(0, 0, f"Active Core\n$R={r_core:.1f}$ cm\n$H={p.active_core_height:.1f}$ cm\n(277 Pins)",
            ha='center', va='center', fontsize=11, fontweight='bold', color='#1B5E20')
    ax.text(0, h_core_half + p.axial_reflector_thickness/2.0, f"Top BeO Cap ($t={p.axial_reflector_thickness:.1f}$ cm)",
            ha='center', va='center', fontsize=9.5, color='#263238')
    ax.text(0, -h_core_half - p.axial_reflector_thickness/2.0, f"Bottom BeO Cap ($t={p.axial_reflector_thickness:.1f}$ cm)",
            ha='center', va='center', fontsize=9.5, color='#263238')
    ax.text(r_core + (r_refl-r_core)/2.0, 0, f"Radial\nBeO\n$t={p.radial_reflector_thickness:.1f}$cm",
            ha='center', va='center', fontsize=9, color='#263238', rotation=-90)

    ax.set_xlim(-r_refl * 1.15, r_refl * 1.15)
    ax.set_ylim(-h_tot_half * 1.12, h_tot_half * 1.12)
    ax.set_aspect('equal')
    ax.set_xlabel("Radius R [cm]", fontsize=11)
    ax.set_ylabel("Axial Height Z [cm]", fontsize=11)
    ax.set_title(f"HALEU NEP Opt-5 - Axial/RZ CSG Geometry Layout\nTotal Height = {p.total_reactor_height:.1f} cm, Outer Radius = {r_refl:.1f} cm", fontsize=11.5, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=8.5)
    fig.tight_layout()
    rz_path = os.path.join(output_dir, "opt5_geometry_rz.png")
    fig.savefig(rz_path)
    fig.savefig(os.path.join("results/opt5/geometry", "opt5_geometry_rz.png"))
    plt.close(fig)
    print(f"Saved: {rz_path}")

    # 2. Fuel Cell Zoom
    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=300)
    pitch = p.moderator_fuel_pitch
    mod_box = patches.Rectangle((-pitch/2, -pitch/2), pitch, pitch, linewidth=1.5, edgecolor='#37474F', facecolor='#FFF9C4', label='YH2 Moderator Matrix')
    c_chan = plt.Circle((0, 0), p.channel_outer_radius, color='#B3E5FC', label=f'He-Xe Coolant Gap ({p.coolant_channel_thickness*10:.2f} mm)')
    c_clad = plt.Circle((0, 0), p.clad_outer_radius, color='#B0BEC5', label=f'Nb-1Zr Cladding ({p.cladding_thickness*10:.2f} mm)')
    c_fuel = plt.Circle((0, 0), p.fuel_radius, color='#C62828', label=f'UN Fuel Pellet ($r={p.fuel_radius:.2f}$ cm)')

    ax.add_patch(mod_box)
    ax.add_patch(c_chan)
    ax.add_patch(c_clad)
    ax.add_patch(c_fuel)

    ax.set_xlim(-pitch * 0.65, pitch * 0.65)
    ax.set_ylim(-pitch * 0.65, pitch * 0.65)
    ax.set_aspect('equal')
    ax.set_xlabel("X [cm]", fontsize=11)
    ax.set_ylabel("Y [cm]", fontsize=11)
    ax.set_title(f"HALEU NEP Opt-5 - Fuel Pin Unit Cell (Pitch = {pitch:.3f} cm)", fontsize=11.5, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.22), ncol=2, framealpha=0.9, fontsize=8.5)
    fig.tight_layout()
    zoom_path = os.path.join(output_dir, "opt5_fuel_cell_zoom.png")
    fig.savefig(zoom_path)
    fig.savefig(os.path.join("results/opt5/geometry", "opt5_fuel_cell_zoom.png"))
    plt.close(fig)
    print(f"Saved: {zoom_path}")

    # 3. Temperature Coefficients Plot
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=300)
    coeffs = ["Fuel Doppler\n(300-1500 K)", "Moderator\n(300-900 K)", "Coolant He-Xe\n(900-1300 K)", "Total Reactor\nTMC"]
    vals = [-1.86, +0.05, -0.67, -2.48]
    errs = [0.24, 0.40, 0.66, 0.81]
    colors = ['#1565C0', '#4CAF50', '#FF9800', '#D32F2F']

    bars = ax.bar(coeffs, vals, yerr=errs, capsize=6, color=colors, edgecolor='black', alpha=0.85, width=0.55)
    ax.axhline(0, color='black', linestyle='-', lw=1.2)
    ax.axhline(-1.5, color='red', linestyle='--', lw=1.2, label=r'Stability Target ($\leq -1.5$ pcm/K)')
    ax.set_ylabel("Temperature Reactivity Coefficient [pcm/K]", fontsize=11)
    ax.set_title("HALEU NEP Opt-5 - Continuous Temperature Reactivity Feedback (OpenMC + WMP)", fontsize=11.5, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6, axis='y')
    ax.legend(loc='lower left', framealpha=0.9, fontsize=9.5)

    for bar, val, err in zip(bars, vals, errs):
        y_pos = val - err - 0.35 if val < 0 else val + err + 0.15
        ax.text(bar.get_x() + bar.get_width()/2.0, y_pos, f"{val:+.2f} ± {err:.2f}",
                ha='center', va='center', fontsize=9.5, fontweight='bold')
    ax.set_ylim(-3.8, 1.2)
    fig.tight_layout()
    tc_path = os.path.join(output_dir, "opt5_temperature_coefficients.png")
    fig.savefig(tc_path)
    fig.savefig(os.path.join("results/opt5/temperature_coefficients", "opt5_temperature_coefficients.png"))
    plt.close(fig)
    print(f"Saved: {tc_path}")

    # 4. Void Reactivity Plot with Explicit Uncertainty Analysis
    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=300)
    drho_void = +164.8
    u_void = 215.5
    ax.errorbar([1], [drho_void], yerr=[u_void], fmt='o', color='#E65100', ecolor='#D84315', elinewidth=2.5, capsize=8, capthick=2, markersize=8, label=f"Opt-5 Coolant Depressurization ({drho_void:+.1f} ± {u_void:.1f} pcm)")
    ax.axhline(0, color='black', linestyle='--', lw=1.5, label=r"Neutral Void Reactivity Limit ($\Delta\rho = 0$ pcm)")
    ax.axhspan(-250, 0, color='#C8E6C9', alpha=0.3, label="Desirable Negative Void Domain")
    ax.axhspan(0, 450, color='#FFCDD2', alpha=0.2, label="Positive Void Domain")

    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(-200, 450)
    ax.set_xticks([1])
    ax.set_xticklabels(["100% Core Depressurization\n(2.0 MPa Operating $\\to$ Vacuum Void)"], fontsize=10.5)
    ax.set_ylabel(r"Reactivity Insertion $\Delta\rho_{void}$ [pcm]", fontsize=11)
    ax.set_title("HALEU NEP Opt-5 - Coolant Void Reactivity & Statistical Boundary", fontsize=11.5, fontweight='bold')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9)
    ax.text(1.0, drho_void + 40, f"{drho_void:+.1f} ± {u_void:.1f} pcm\n(Spans zero boundary: -50.7 to +380.3 pcm)",
            ha='center', va='bottom', fontsize=9.5, fontweight='bold', color='#BF360C')
    fig.tight_layout()
    void_path = os.path.join(output_dir, "opt5_void_reactivity.png")
    fig.savefig(void_path)
    fig.savefig(os.path.join("results/opt5/void", "opt5_void_reactivity.png"))
    plt.close(fig)
    print(f"Saved: {void_path}")

    # 5. Fissile Inventory & Pu-239 Breeding Trajectory
    csv_file = "results/opt5/depletion/opt5_criticality_depletion_400kw.csv"
    times = []
    u235_m = []
    pu239_m = []
    with open(csv_file, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["Time [years]"]))
            u235_m.append(float(row["U-235 Mass [kg]"]))
            pu239_m.append(float(row["Pu-239 Mass [g]"]))

    fig, ax1 = plt.subplots(figsize=(8.5, 4.8), dpi=300)
    color_u = '#B71C1C'
    ax1.set_xlabel("Operational Lifetime [Years]", fontsize=11)
    ax1.set_ylabel("Fissile $^{235}$U Mass [kg]", color=color_u, fontsize=11)
    ax1.plot(times, u235_m, 'o-', color=color_u, lw=2.2, label="$^{235}$U Fissile Inventory")
    ax1.tick_params(axis='y', labelcolor=color_u)
    ax1.axvline(10.0, color='green', linestyle='-.', lw=1.5, label="10-Yr Mandatory Milestone")
    ax1.axvline(20.0, color='purple', linestyle='--', lw=1.5, label="20-Yr Achievable Lifetime Limit")
    ax1.grid(True, linestyle=':', alpha=0.6)

    ax2 = ax1.twinx()
    color_pu = '#1B5E20'
    ax2.set_ylabel("Transuranic $^{239}$Pu Bred [g]", color=color_pu, fontsize=11)
    ax2.plot(times, pu239_m, 's--', color=color_pu, lw=2.0, label="$^{239}$Pu Bred Mass")
    ax2.tick_params(axis='y', labelcolor=color_pu)

    fig.suptitle("HALEU NEP Opt-5 - 20-Year Fissile Depletion & Transuranic Breeding (400 kWth)", fontsize=11.5, fontweight='bold')
    fig.tight_layout()
    iso_path = os.path.join(output_dir, "opt5_u235_and_plutonium_inventory.png")
    fig.savefig(iso_path)
    fig.savefig(os.path.join("results/opt5/depletion", "opt5_u235_and_plutonium_inventory.png"))
    plt.close(fig)
    print(f"Saved: {iso_path}")

if __name__ == "__main__":
    generate_additional_plots()
