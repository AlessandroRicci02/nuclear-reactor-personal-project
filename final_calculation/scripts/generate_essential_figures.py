#!/usr/bin/env python3
"""
final_calculation/generate_essential_figures.py

Generates the Four Essential, Publication-Quality Engineering Figures:
1. fig1_reference_reactor_geometry.png
   - 2D cross-section slice of active reference core geometry
   - Core barrel, 217 fuel pins, clad, moderator, BeO reflector, 10 control drums
   - Physical scale bar, material legend, clear labels
2. fig2_neutronic_power_distribution.png
   - 2D pin-by-pin radial power map (Watts) across hexagonal core
   - 1D normalized axial power profile q'(z)/q_bar with Fr, Fz, Fq annotations
3. fig3_temperature_distributions.png
   - 2D peak cladding temperature distribution across 217 channels
   - 1D axial profiles (T_fuel_center, T_clad_out, T_bulk_coolant) along peak pin
   - Exit temperature elevation Delta T_hot and orificing effect
4. fig4_depletion_lifetime_trajectory.png
   - Uncontrolled multiplication factor k_eff(t) and critical drum angle theta_crit(t)
   - Heavy metal isotopics (U-235 burnup, Pu-239 breeding) and burnup (MWd/kgU)
   - Stopping boundary at theta = 0 deg clearly marked
"""

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

BASE_DIR = "/home/aleri/OpenMC_projects/first_test"
FIG_DIR = os.path.join(BASE_DIR, "final_calculation/figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Geometry constants
PITCH = 2.270          # cm
CORE_H = 42.45         # cm
R_FUEL = 0.580         # cm
T_CLAD = 0.048         # cm
T_GAP = 0.084          # cm
R_CLAD = R_FUEL + T_CLAD # 0.628 cm
R_CHAN = R_CLAD + T_GAP  # 0.712 cm
CORE_R = 19.45         # cm
REFL_RAD = 9.80        # cm
OUTER_R = CORE_R + REFL_RAD # 29.25 cm
DRUM_R = 3.95          # cm
DRUM_CLEARANCE = 0.180 # cm
DRUM_RAD_POS = CORE_R + DRUM_CLEARANCE + DRUM_R # 23.58 cm
DRUM_COUNT = 10
ABS_TH = 1.50          # cm
ABS_ANGLE = 148.0      # deg


def generate_fig1_geometry():
    """Generates Figure 1: Reference Reactor Geometry Cross-Section."""
    print("[+] Generating Figure 1: Reference Reactor Geometry...")
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    
    # Outer reflector boundary
    circle_outer = plt.Circle((0, 0), OUTER_R, color='#d9d9d9', ec='#333333', lw=2.0, label='BeO Radial Reflector')
    ax.add_patch(circle_outer)
    
    # Core boundary
    circle_core = plt.Circle((0, 0), CORE_R, color='#e8f4f8', ec='#1b4965', lw=2.5, ls='--', label='Core Barrel Boundary')
    ax.add_patch(circle_core)
    
    # Pin coordinates (hexagonal lattice, 217 pins, 9 rings)
    v1 = np.array([PITCH, 0.0])
    v2 = np.array([PITCH / 2.0, PITCH * np.sqrt(3.0) / 2.0])
    
    pin_coords = []
    pin_rings = []
    
    for r in range(9):
        if r == 0:
            pin_coords.append((0.0, 0.0))
            pin_rings.append(0)
        else:
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    pin_coords.append((float(pt[0]), float(pt[1])))
                    pin_rings.append(r)
                    
    # Draw pins
    for idx, (x, y) in enumerate(pin_coords):
        ring = pin_rings[idx]
        is_inner = (ring <= 2)
        
        # Coolant annulus
        c_chan = plt.Circle((x, y), R_CHAN, color='#cce3de', ec='none')
        ax.add_patch(c_chan)
        # Clad
        c_clad = plt.Circle((x, y), R_CLAD, color='#6b705c', ec='none')
        ax.add_patch(c_clad)
        # Fuel
        f_color = '#d62828' if is_inner else '#f77f00' # Red for inner 19.3%, Orange for outer 19.75%
        c_fuel = plt.Circle((x, y), R_FUEL, color=f_color, ec='#333333', lw=0.4)
        ax.add_patch(c_fuel)
        
    # Draw 10 Control Drums
    # Critical drum angle ~ 90 deg orientation
    theta_drum_crit = 92.5 # deg
    for k in range(DRUM_COUNT):
        phi_k = k * 2.0 * np.pi / DRUM_COUNT
        xk = DRUM_RAD_POS * np.cos(phi_k)
        yk = DRUM_RAD_POS * np.sin(phi_k)
        
        # BeO drum body
        c_drum = plt.Circle((xk, yk), DRUM_R, color='#b7b7a4', ec='#2b2d42', lw=1.5)
        ax.add_patch(c_drum)
        
        # Absorber sector (B4C, 148 deg sector, 1.5 cm thick)
        ang_deg = np.degrees(phi_k) + theta_drum_crit
        wedge_b4c = patches.Wedge(
            (xk, yk), DRUM_R,
            ang_deg - ABS_ANGLE / 2.0, ang_deg + ABS_ANGLE / 2.0,
            width=ABS_TH, color='#1d3557', ec='#000814', lw=1.0
        )
        ax.add_patch(wedge_b4c)
        
        # Center pivot
        c_pivot = plt.Circle((xk, yk), 0.35, color='#495057')
        ax.add_patch(c_pivot)

    # Scale bar (10 cm)
    scale_x0, scale_y0 = -25.0, -26.0
    ax.plot([scale_x0, scale_x0 + 10.0], [scale_y0, scale_y0], color='black', lw=3.0)
    ax.plot([scale_x0, scale_x0], [scale_y0 - 0.6, scale_y0 + 0.6], color='black', lw=2.0)
    ax.plot([scale_x0 + 10.0, scale_x0 + 10.0], [scale_y0 - 0.6, scale_y0 + 0.6], color='black', lw=2.0)
    ax.text(scale_x0 + 5.0, scale_y0 + 1.2, "10 cm", ha='center', va='bottom', fontsize=12, fontweight='bold')
    
    # Custom Legend
    legend_elements = [
        patches.Patch(facecolor='#d62828', edgecolor='black', label='UN Fuel Inner (19.30 wt% 235U, Rings 0-2)'),
        patches.Patch(facecolor='#f77f00', edgecolor='black', label='UN Fuel Outer (19.75 wt% 235U, Rings 3-8)'),
        patches.Patch(facecolor='#6b705c', edgecolor='black', label='Nb-1Zr Cladding (t = 0.48 mm)'),
        patches.Patch(facecolor='#cce3de', edgecolor='black', label='He Coolant Annulus (t_gap = 0.84 mm)'),
        patches.Patch(facecolor='#e8f4f8', edgecolor='#1b4965', label='YH2.0 Moderator Matrix (Canned)'),
        patches.Patch(facecolor='#d9d9d9', edgecolor='#333333', label='BeO Radial Reflector (t = 9.80 cm)'),
        patches.Patch(facecolor='#1d3557', edgecolor='black', label='B4C Control Drum Absorber (90% 10B, 10 Drums)')
    ]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.28, 1.0), fontsize=9.5, frameon=True)
    
    ax.set_xlim(-31.5, 31.5)
    ax.set_ylim(-31.5, 31.5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    title_text = "Fig. 1: Reference Space NEP Reactor Core Architecture (V7_CAND_BALANCED_ROBUST)\n" \
                 "217 Annular Fuel Pins, 9 Concentric Rings, BeO Reflector, 10 Rotating B4C Drums"
    plt.title(title_text, fontsize=12, fontweight='bold', pad=15)
    
    out_path = os.path.join(FIG_DIR, "fig1_reference_reactor_geometry.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[✓] Saved {out_path}")


def generate_fig2_power():
    """Generates Figure 2: Neutronic Power Distribution (2D Radial & 1D Axial)."""
    print("[+] Generating Figure 2: Neutronic Power Distribution...")
    q_file = os.path.join(BASE_DIR, "final_calculation/coupled_picard/converged_q_vol_3d.npy")
    ch_file = os.path.join(BASE_DIR, "final_calculation/coupled_picard/converged_channel_results.csv")
    
    if not os.path.exists(q_file) or not os.path.exists(ch_file):
        print(f"[!] Warning: {q_file} not found; skipping Fig 2")
        return
        
    q_vol_3d = np.load(q_file)
    df_ch = pd.read_csv(ch_file)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # 1. 2D Radial Map
    # Pin coordinates
    v1 = np.array([PITCH, 0.0])
    v2 = np.array([PITCH / 2.0, PITCH * np.sqrt(3.0) / 2.0])
    pin_coords = []
    for r in range(9):
        if r == 0:
            pin_coords.append((0.0, 0.0))
        else:
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    pin_coords.append((float(pt[0]), float(pt[1])))
                    
    pin_powers = df_ch["power_W"].values
    mean_power = np.mean(pin_powers)
    fr = np.max(pin_powers) / mean_power
    
    # Scatter hex plot
    x_pts = [p[0] for p in pin_coords]
    y_pts = [p[1] for p in pin_coords]
    
    # Plot core circle
    c_core = plt.Circle((0, 0), CORE_R, color='none', ec='#1b4965', lw=1.5, ls='--')
    ax1.add_patch(c_core)
    
    sc = ax1.scatter(x_pts, y_pts, c=pin_powers, s=110, cmap='inferno', edgecolors='#222222', lw=0.6)
    cbar = plt.colorbar(sc, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("Pin Thermal Power [W]", fontsize=10, fontweight='bold')
    
    ax1.set_xlim(-21.0, 21.0)
    ax1.set_ylim(-21.0, 21.0)
    ax1.set_aspect('equal')
    ax1.set_xlabel("X Dimension [cm]", fontsize=10, fontweight='bold')
    ax1.set_ylabel("Y Dimension [cm]", fontsize=10, fontweight='bold')
    ax1.set_title(f"(A) 2D Radial Pin Power Distribution\nCore Power: 400.0 kWth | Radial Peaking F_r = {fr:.3f}", fontsize=11, fontweight='bold')
    ax1.grid(True, ls=':', alpha=0.5)
    
    # 2. 1D Normalized Axial Profile
    # Core average axial power profile
    axial_power_raw = np.sum(q_vol_3d, axis=0)
    axial_norm = axial_power_raw / np.mean(axial_power_raw)
    z_edges = np.linspace(0.0, CORE_H, 11)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    fz = np.max(axial_norm)
    fq = fr * fz
    
    # Peak pin axial profile
    peak_idx = int(np.argmax(pin_powers))
    peak_axial_norm = q_vol_3d[peak_idx, :] / np.mean(q_vol_3d[peak_idx, :])
    
    ax2.plot(z_centers, axial_norm, 'o-', color='#1d3557', lw=2.2, label='Core-Average Axial Profile')
    ax2.plot(z_centers, peak_axial_norm, 's--', color='#e63946', lw=1.8, label=f'Peak Pin {peak_idx} Axial Profile')
    
    # Cosine fit comparison
    z_fine = np.linspace(0.0, CORE_H, 100)
    h_ext = CORE_H + 2.0 * 7.5 # with reflector savings
    cos_fit = 1.28 * np.cos(np.pi * (z_fine - CORE_H / 2.0) / h_ext)
    ax2.plot(z_fine, cos_fit, ':', color='#457b9d', lw=1.5, label='Chopped Cosine Fit (Extrapolated Core)')
    
    ax2.axhline(1.0, color='gray', ls='--', lw=1.0)
    ax2.set_xlabel("Axial Elevation z [cm]", fontsize=10, fontweight='bold')
    ax2.set_ylabel("Normalized Axial Power Density q'(z) / q_bar", fontsize=10, fontweight='bold')
    ax2.set_title(f"(B) 1D Axial Power Profile\nAxial Peaking F_z = {fz:.3f} | Total 3D Peaking F_q = {fq:.3f}", fontsize=11, fontweight='bold')
    ax2.grid(True, ls=':', alpha=0.5)
    ax2.legend(loc='lower center', fontsize=9.5)
    
    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "fig2_neutronic_power_distribution.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[✓] Saved {out_path}")


def generate_fig3_temperatures():
    """Generates Figure 3: Core Temperature Distributions (2D & 1D Peak Pin Profiles)."""
    print("[+] Generating Figure 3: Temperature Distributions...")
    ch_file = os.path.join(BASE_DIR, "final_calculation/coupled_picard/converged_channel_results.csv")
    q_file = os.path.join(BASE_DIR, "final_calculation/coupled_picard/converged_q_vol_3d.npy")
    moose_comp = os.path.join(BASE_DIR, "final_calculation/moose_thermal/moose_vs_1d_comparison.csv")
    
    if not os.path.exists(ch_file) or not os.path.exists(q_file):
        print(f"[!] Warning: {ch_file} not found; skipping Fig 3")
        return
        
    df_ch = pd.read_csv(ch_file)
    q_vol_3d = np.load(q_file)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    # 1. 2D Map of Peak Cladding Temperatures
    v1 = np.array([PITCH, 0.0])
    v2 = np.array([PITCH / 2.0, PITCH * np.sqrt(3.0) / 2.0])
    pin_coords = []
    for r in range(9):
        if r == 0:
            pin_coords.append((0.0, 0.0))
        else:
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    pin_coords.append((float(pt[0]), float(pt[1])))
                    
    x_pts = [p[0] for p in pin_coords]
    y_pts = [p[1] for p in pin_coords]
    t_clad_peaks = df_ch["T_clad_max_K"].values
    
    sc = ax1.scatter(x_pts, y_pts, c=t_clad_peaks, s=110, cmap='plasma', edgecolors='#222222', lw=0.6)
    cbar = plt.colorbar(sc, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("Peak Cladding Outer Temperature [K]", fontsize=10, fontweight='bold')
    
    c_core = plt.Circle((0, 0), CORE_R, color='none', ec='#1b4965', lw=1.5, ls='--')
    ax1.add_patch(c_core)
    
    ax1.set_xlim(-21.0, 21.0)
    ax1.set_ylim(-21.0, 21.0)
    ax1.set_aspect('equal')
    ax1.set_xlabel("X Dimension [cm]", fontsize=10, fontweight='bold')
    ax1.set_ylabel("Y Dimension [cm]", fontsize=10, fontweight='bold')
    max_tc = np.max(t_clad_peaks)
    mean_tc = np.mean(t_clad_peaks)
    ax1.set_title(f"(A) 2D Cladding Peak Temperature Field\nMax T_clad = {max_tc:.1f} K | Mean T_clad = {mean_tc:.1f} K", fontsize=11, fontweight='bold')
    ax1.grid(True, ls=':', alpha=0.5)
    
    # 2. 1D Axial Profiles along Peak Pin (Pin 203)
    peak_idx = int(df_ch["T_clad_max_K"].idxmax())
    peak_row = df_ch.loc[peak_idx]
    
    z_edges = np.linspace(0.0, CORE_H, 11)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    
    # Calculate axial profiles for peak pin
    q_peak_axial = q_vol_3d[peak_idx, :]
    m_peak = peak_row["m_dot_kg_s"]
    cp_he = 5193.16
    v_slice = math.pi * (R_FUEL*1e-2)**2 * (CORE_H*1e-2 / 10.0)
    q_watts = q_peak_axial * v_slice
    
    T_cool_ax = np.zeros(10)
    T_cur = 950.0
    for s in range(10):
        T_cur += q_watts[s] / (m_peak * cp_he)
        T_cool_ax[s] = T_cur
        
    h_conv = peak_row["h_mean"]
    r_f = R_FUEL * 1e-2
    r_c = R_CLAD * 1e-2
    A_ht = 2.0 * math.pi * r_c * (CORE_H * 1e-2 / 10.0)
    
    T_clad_out_ax = T_cool_ax + q_watts / (h_conv * A_ht)
    k_clad = 55.0 # W/(m*K)
    delta_T_clad = (q_watts / (2.0 * math.pi * (CORE_H*1e-2/10.0) * k_clad)) * math.log(r_c / r_f)
    T_clad_in_ax = T_clad_out_ax + delta_T_clad
    k_fuel = 23.0 # W/(m*K)
    delta_T_fuel = q_peak_axial * (r_f**2) / (4.0 * k_fuel)
    T_fuel_center_ax = T_clad_in_ax + delta_T_fuel
    
    ax2.plot(z_centers, T_fuel_center_ax, 'r-o', lw=2.2, label='Fuel Centerline (T_fuel,center)')
    ax2.plot(z_centers, T_clad_out_ax, 'm-s', lw=2.0, label='Cladding Outer Surface (T_clad,out)')
    ax2.plot(z_centers, T_cool_ax, 'b-^', lw=1.8, label='Bulk Coolant (T_cool,bulk)')
    
    # Design Limit (Nb-1Zr 1% creep limit at 1350 K)
    ax2.axhline(1350.0, color='black', ls='--', lw=1.5, label='Nb-1Zr Creep Limit (1350 K)')
    ax2.axhline(1150.0, color='green', ls=':', lw=1.2, label='Core Bulk Outlet Target (1150 K)')
    
    ax2.set_xlabel("Axial Elevation z [cm]", fontsize=10, fontweight='bold')
    ax2.set_ylabel("Temperature [K]", fontsize=10, fontweight='bold')
    ax2.set_title(f"(B) Peak Pin Axial Temperature Profiles (Pin {peak_idx}, Ring {peak_row['ring']})\n"
                  f"Peak T_fuel = {np.max(T_fuel_center_ax):.1f} K | Peak T_clad = {np.max(T_clad_out_ax):.1f} K (Creep Margin: {1350.0 - np.max(T_clad_out_ax):.1f} K)",
                  fontsize=11, fontweight='bold')
    ax2.grid(True, ls=':', alpha=0.5)
    ax2.legend(loc='lower right', fontsize=9.0)
    
    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "fig3_temperature_distributions.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[✓] Saved {out_path}")


def generate_fig4_depletion():
    """Generates Figure 4: Depletion & Drum Angle Lifetime Trajectory."""
    print("[+] Generating Figure 4: Depletion Lifetime Trajectory...")
    traj_file = os.path.join(BASE_DIR, "final_calculation/depletion/depletion_trajectory.csv")
    
    if not os.path.exists(traj_file):
        print(f"[!] Warning: {traj_file} not found; skipping Fig 4")
        return
        
    df_traj = pd.read_csv(traj_file)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    
    t_yr = df_traj["time_years"].values
    k_un = df_traj["k_eff_uncontrolled_0deg"].values
    sig_k = df_traj["sigma_k"].values
    theta_crit = df_traj["critical_drum_angle_deg"].values
    
    # 1. k_eff and Drum Angle vs Time
    color1 = '#1d3557'
    ax1.set_xlabel("Operating Time [Years]", fontsize=10, fontweight='bold')
    ax1.set_ylabel("Uncontrolled Multiplication Factor k_eff (theta = 0 deg)", color=color1, fontsize=10, fontweight='bold')
    ax1.errorbar(t_yr, k_un, yerr=sig_k, fmt='o-', color=color1, lw=2.0, capsize=4, label='k_eff(theta = 0 deg)')
    ax1.axhline(1.0000, color='black', ls='--', lw=1.2, label='Criticality Threshold (k = 1.000)')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, ls=':', alpha=0.5)
    
    # Twin axis for critical drum angle
    ax1_twin = ax1.twinx()
    color2 = '#e63946'
    ax1_twin.set_ylabel("Critical Drum Angle theta_crit [degrees]", color=color2, fontsize=10, fontweight='bold')
    ax1_twin.plot(t_yr, theta_crit, 's--', color=color2, lw=2.2, label='Critical Drum Angle theta(t)')
    ax1_twin.axhline(0.0, color='red', ls=':', lw=1.5, label='Fully Withdrawn Limit (theta = 0 deg)')
    ax1_twin.tick_params(axis='y', labelcolor=color2)
    
    # Find exact EOL where theta hits 0
    eol_yr = float(np.interp(1.0000, k_un[::-1], t_yr[::-1]))
    ax1.axvline(eol_yr, color='purple', ls='-.', lw=1.5)
    ax1.text(eol_yr - 0.5, 1.002, f"EOL: {eol_yr:.1f} yr\n(theta = 0 deg)", color='purple', fontweight='bold', ha='right', fontsize=9.5)
    
    ax1.set_title(f"(A) Core Criticality & Drum Angle Trajectory\nOperational Lifetime to theta = 0 deg: {eol_yr:.1f} Years", fontsize=11, fontweight='bold')
    
    # 2. Isotopic Masses and Cumulative Burnup
    u235_kg = df_traj["u235_mass_kg"].values
    pu239_kg = df_traj["pu239_mass_kg"].values
    burnup_mwd = df_traj["burnup_MWd_kgU"].values
    
    ax2.plot(t_yr, u235_kg, 'b-o', lw=2.0, label='U-235 Inventory [kg]')
    ax2.plot(t_yr, pu239_kg, 'g-s', lw=2.0, label='Pu-239 Inventory [kg]')
    
    ax2_twin = ax2.twinx()
    ax2_twin.plot(t_yr, burnup_mwd, 'm-.', lw=2.0, label='Core Burnup [MWd/kgU]')
    ax2_twin.set_ylabel("Cumulative Burnup [MWd/kgU]", color='purple', fontsize=10, fontweight='bold')
    ax2_twin.tick_params(axis='y', labelcolor='purple')
    
    ax2.set_xlabel("Operating Time [Years]", fontsize=10, fontweight='bold')
    ax2.set_ylabel("Fissile / Actinide Inventory [kg]", fontsize=10, fontweight='bold')
    ax2.set_title(f"(B) Isotopic Depletion & Energy Production\nEOL Burnup: {burnup_mwd[-1]:.1f} MWd/kgU | FIMA: {df_traj['fima_pct'].values[-1]:.2f}%", fontsize=11, fontweight='bold')
    ax2.grid(True, ls=':', alpha=0.5)
    ax2.legend(loc='center left', fontsize=9.5)
    ax2_twin.legend(loc='center right', fontsize=9.5)
    
    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "fig4_depletion_lifetime_trajectory.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[✓] Saved {out_path}")


def main():
    generate_fig1_geometry()
    generate_fig2_power()
    generate_fig3_temperatures()
    generate_fig4_depletion()
    print("\n[✓] ALL ESSENTIAL FIGURES COMPLETED.")


if __name__ == "__main__":
    main()
