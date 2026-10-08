#!/usr/bin/env python3
"""
generate_high_precision_plots_c01.py

Ultra-High-Precision Technical CAD and CSG Figure Generator (300 DPI)
for Winning Reactor C01 (matrix_dual_pass):
(217 Pins, UN Fuel, YH2.0 Matrix, 10 BeO/B4C Drums, 400 kWth, 100 kWe).

Generates 6 Authoritative High-Resolution Figures:
1. 01_reactor_radial_xy.png: Full core radial cross section (217 pins, 6 aux channels/cell,
   10 rotating drums, BeO reflector, dimensions and custom technical legend).
2. 02_reactor_axial_xz.png: Full core axial cross section (active core, BeO axial reflectors,
   fuel pins, control drums, dimension callouts with arrows).
3. 03_unit_cell_5layer_detailed.png: Close-up 5-layer unit cell (UN, Nb-1Zr, He gap, Al2O3 liner,
   YH2 matrix with 6x aux channels, exact dimensions and callout annotations).
4. 04_control_drum_mechanism.png: Control drum mechanism and rotation trajectory (Drums OUT 0° vs IN 180°,
   155° B4C absorber arc, 1.8 cm thickness, flux direction, shutdown physics).
5. 05_openmc_csg_slice_xy.png: Ultra-high-resolution OpenMC CSG pixel raytrace (3000x3000)
   with customized high-contrast material palette.
6. 06_openmc_csg_slice_xz.png: Ultra-high-resolution OpenMC CSG axial pixel raytrace (2400x3000).
"""

import os
import sys
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches

ROOT_SCRIPT = Path(__file__).resolve()
REACTOR_DIR = ROOT_SCRIPT.parents[1]
FIG_DIR = REACTOR_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Also define artifact directory
ARTIFACT_DIR = Path("/home/aleri/.gemini/antigravity/brain/0c8a66db-b02e-44e7-8d87-eec426335ff7")

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "figure.titlesize": 14
})

# Geometric constants for C01
PITCH = 2.280       # cm
H_CORE = 46.00      # cm
REFL_AX_TH = 8.80   # cm
REFL_RAD_TH = 9.30  # cm
R_FUEL = 0.585      # cm (5.85 mm)
T_CLAD = 0.048      # cm (0.48 mm)
R_CLAD = R_FUEL + T_CLAD # 0.633 cm
T_GAP = 0.084       # cm (0.84 mm)
R_CHAN = R_CLAD + T_GAP  # 0.717 cm
T_AL2O3 = 0.010     # cm (100 um)
R_AL2O3 = R_CHAN + T_AL2O3 # 0.727 cm
R_ADD = 0.165       # cm (1.65 mm, 6 aux channels/cell)
R_ADD_POS = 0.910 * (PITCH / 2.270) # 0.914 cm
R_CORE = 19.54      # cm
R_REFL = R_CORE + REFL_RAD_TH # 28.84 cm
DRUM_R = 3.95       # cm
ABS_TH = 1.80       # cm (B4C thickness)
ABS_ANGLE = 155.0   # deg (B4C arc)
DRUM_CLEARANCE = 0.180 # cm
DRUM_COUNT = 10
DRUM_RAD_POS = R_CORE + DRUM_CLEARANCE + DRUM_R # 23.67 cm


# =============================================================================
# 1. FULL RADIAL CROSS SECTION (XY MIDPLANE)
# =============================================================================
def plot_radial_xy(out_path):
    fig, ax = plt.subplots(figsize=(11, 11), dpi=300)

    # 1. BeO Radial Reflector Outer Circle
    c_outer = plt.Circle((0, 0), R_REFL, facecolor="#e0e0e0", edgecolor="#333333",
                         linewidth=1.8, label="BeO Radial Reflector (9.30 cm)")
    ax.add_patch(c_outer)

    # 2. Active Core Barrel (YH2 Matrix Boundary)
    c_core = plt.Circle((0, 0), R_CORE, facecolor="#f5efe6", edgecolor="#222222",
                        linewidth=2.0, label="Active Core (YH2.0 Moderator Matrix)")
    ax.add_patch(c_core)

    # 3. Hexagonal Lattice Pin Generation (8 rings = 217 pins)
    v1 = np.array([PITCH, 0.0])
    v2 = np.array([PITCH / 2.0, PITCH * np.sqrt(3.0) / 2.0])

    pins_pos = []
    for i in range(-8, 9):
        for j in range(-8, 9):
            ring = max(abs(i), abs(j), abs(i + j))
            if ring <= 8:
                pos = i * v1 + j * v2
                pins_pos.append((pos[0], pos[1], ring))

    # Plot pin cells
    for x, y, ring in pins_pos:
        is_inner = (ring <= 2)
        # Al2O3 Liner
        c_al = plt.Circle((x, y), R_AL2O3, facecolor="#ffe0b2", edgecolor="#fb8c00", linewidth=0.5)
        ax.add_patch(c_al)
        # Annular Coolant Gap (He)
        c_cool = plt.Circle((x, y), R_CHAN, facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.4)
        ax.add_patch(c_cool)
        # Nb-1Zr Cladding
        c_cl = plt.Circle((x, y), R_CLAD, facecolor="#cfd8dc", edgecolor="#546e7a", linewidth=0.4)
        ax.add_patch(c_cl)
        # UN Fuel Pellet (Inner vs Outer enrichment)
        fuel_col = "#e53935" if is_inner else "#8e24aa"
        c_f = plt.Circle((x, y), R_FUEL, facecolor=fuel_col, edgecolor="#212121", linewidth=0.5)
        ax.add_patch(c_f)

        # 6 Auxiliary Coolant Channels in YH2 matrix
        for k in range(6):
            ang = math.radians(k * 60.0)
            ax_x = x + R_ADD_POS * math.cos(ang)
            ax_y = y + R_ADD_POS * math.sin(ang)
            c_ax = plt.Circle((ax_x, ax_y), R_ADD, facecolor="#00bcd4", edgecolor="#006064", linewidth=0.4)
            ax.add_patch(c_ax)

    # 4. Control Drums in Radial Reflector (Drums OUT position: theta = 0°, B4C outward)
    for k in range(DRUM_COUNT):
        phi_k = k * 2.0 * math.pi / DRUM_COUNT
        xk = DRUM_RAD_POS * math.cos(phi_k)
        yk = DRUM_RAD_POS * math.sin(phi_k)

        # Drum BeO cylinder
        c_drum = plt.Circle((xk, yk), DRUM_R, facecolor="#b0bec5", edgecolor="#37474f", linewidth=1.4)
        ax.add_patch(c_drum)

        # B4C Absorber sector facing outward (0 deg)
        rot_deg = math.degrees(phi_k)
        half_ang = ABS_ANGLE / 2.0
        theta1 = rot_deg - half_ang
        theta2 = rot_deg + half_ang

        wedge = patches.Wedge((xk, yk), DRUM_R, theta1, theta2, width=ABS_TH,
                              facecolor="#212121", edgecolor="black", linewidth=1.2)
        ax.add_patch(wedge)

        # Drum center marker
        ax.plot(xk, yk, 'k+', ms=6)

    # Core Dimensions and Annotations
    ax.annotate("", xy=(-R_REFL, -R_REFL * 1.05), xytext=(R_REFL, -R_REFL * 1.05),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.5))
    ax.text(0.0, -R_REFL * 1.10, f"Total Reflector Diameter 2*R_ref = {2*R_REFL:.2f} cm",
            ha="center", fontsize=10.5, weight="bold")

    ax.annotate("", xy=(-R_CORE, R_REFL * 1.05), xytext=(R_CORE, R_REFL * 1.05),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.5))
    ax.text(0.0, R_REFL * 1.09, f"Active Core Diameter 2*R_core = {2*R_CORE:.2f} cm",
            ha="center", fontsize=10.5, weight="bold")

    ax.set_xlim(-R_REFL * 1.18, R_REFL * 1.48)
    ax.set_ylim(-R_REFL * 1.18, R_REFL * 1.18)
    ax.set_aspect("equal")
    ax.set_xlabel("X Position [cm]", fontsize=12, weight="bold")
    ax.set_ylabel("Y Position [cm]", fontsize=12, weight="bold")
    ax.set_title("C01 NEP REACTOR: FULL RADIAL CROSS-SECTION (XY Midplane, z=0)\n"
                 "217 UN Pins (8 Rings), YH2.0 Matrix + 6 Aux Channels/Cell (r=1.65 mm), 10 BeO/B4C Drums (OUT 0°)",
                 fontsize=12, weight="bold")

    # Technical Legend
    legend_elements = [
        patches.Patch(facecolor="#e53935", edgecolor="black", label="UN Fuel Inner (19.60 wt% U-235, 19 pins)"),
        patches.Patch(facecolor="#8e24aa", edgecolor="black", label="UN Fuel Outer (19.75 wt% U-235, 198 pins)"),
        patches.Patch(facecolor="#cfd8dc", edgecolor="#546e7a", label="Nb-1Zr Cladding (t = 0.48 mm)"),
        patches.Patch(facecolor="#b3e5fc", edgecolor="#0288d1", label="He Annular Hot Channel (t = 0.84 mm)"),
        patches.Patch(facecolor="#ffe0b2", edgecolor="#fb8c00", label="Al2O3 Hydrogen Barrier (t = 100 um)"),
        patches.Patch(facecolor="#00bcd4", edgecolor="#006064", label="6x Aux Coolant Channels (r = 1.65 mm)"),
        patches.Patch(facecolor="#f5efe6", edgecolor="#222222", label="YH2.0 Solid Moderator Matrix"),
        patches.Patch(facecolor="#e0e0e0", edgecolor="#333333", label="BeO Radial Reflector (t = 9.30 cm)"),
        patches.Patch(facecolor="#b0bec5", edgecolor="#37474f", label="BeO Control Drum Cylinder (r = 3.95 cm)"),
        patches.Patch(facecolor="#212121", edgecolor="black", label="B4C Absorber Sector (155°, t = 1.80 cm)")
    ]
    ax.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1.45, 1.0), fontsize=9.5, frameon=True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_path}")


# =============================================================================
# 2. FULL AXIAL CROSS SECTION (XZ PLANE)
# =============================================================================
def plot_axial_xz(out_path):
    fig, ax = plt.subplots(figsize=(10.5, 8.5), dpi=300)

    z_core_bot = -H_CORE / 2.0
    z_core_top = +H_CORE / 2.0
    z_tot_bot = z_core_bot - REFL_AX_TH
    z_tot_top = z_core_top + REFL_AX_TH

    # 1. BeO Radial Reflector Blocks (Left & Right)
    rect_refl_l = patches.Rectangle((-R_REFL, z_tot_bot), R_REFL - R_CORE, z_tot_top - z_tot_bot,
                                    facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5, label="BeO Radial Reflector (9.30 cm)")
    rect_refl_r = patches.Rectangle((+R_CORE, z_tot_bot), R_REFL - R_CORE, z_tot_top - z_tot_bot,
                                    facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5)
    ax.add_patch(rect_refl_l)
    ax.add_patch(rect_refl_r)

    # 2. Control Drums (spanning core height in reflector)
    for x_d in [-DRUM_RAD_POS, DRUM_RAD_POS]:
        rect_drum = patches.Rectangle((x_d - DRUM_R, z_core_bot), 2 * DRUM_R, H_CORE,
                                      facecolor="#b0bec5", edgecolor="#37474f", linewidth=1.2)
        ax.add_patch(rect_drum)
        sign = 1 if x_d > 0 else -1
        rect_abs = patches.Rectangle((x_d + sign * (DRUM_R - ABS_TH), z_core_bot), sign * (-ABS_TH), H_CORE,
                                     facecolor="#212121", edgecolor="black", linewidth=1.0)
        ax.add_patch(rect_abs)

    # 3. Top and Bottom BeO Axial Reflectors
    rect_ax_bot = patches.Rectangle((-R_CORE, z_tot_bot), 2 * R_CORE, REFL_AX_TH,
                                    facecolor="#d7ccc8", edgecolor="#5d4037", linewidth=1.2,
                                    label="Bottom BeO Axial Reflector (8.80 cm)")
    rect_ax_top = patches.Rectangle((-R_CORE, z_core_top), 2 * R_CORE, REFL_AX_TH,
                                    facecolor="#d7ccc8", edgecolor="#5d4037", linewidth=1.2,
                                    label="Top BeO Axial Reflector (8.80 cm)")
    ax.add_patch(rect_ax_bot)
    ax.add_patch(rect_ax_top)

    # 4. Active Core YH2.0 Matrix
    rect_core = patches.Rectangle((-R_CORE, z_core_bot), 2 * R_CORE, H_CORE,
                                  facecolor="#f5efe6", edgecolor="#222222", linewidth=1.6,
                                  label="Active Core (YH2.0 Solid Matrix)")
    ax.add_patch(rect_core)

    # 5. Representative Fuel Pins along X axis
    pin_xs = np.linspace(-R_CORE * 0.88, R_CORE * 0.88, 17)
    for px in pin_xs:
        # Cladding
        rect_clad = patches.Rectangle((px - R_CLAD, z_core_bot), 2 * R_CLAD, H_CORE,
                                      facecolor="#cfd8dc", edgecolor="#546e7a", linewidth=0.4)
        ax.add_patch(rect_clad)
        # Fuel Pellet
        rect_fuel = patches.Rectangle((px - R_FUEL, z_core_bot), 2 * R_FUEL, H_CORE,
                                      facecolor="#8e24aa" if abs(px) > 6.0 else "#e53935",
                                      edgecolor="#212121", linewidth=0.5)
        ax.add_patch(rect_fuel)
        # Axial coolant ducts through reflectors
        rect_duct_b = patches.Rectangle((px - 0.45, z_tot_bot), 0.90, REFL_AX_TH,
                                        facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.3)
        rect_duct_t = patches.Rectangle((px - 0.45, z_core_top), 0.90, REFL_AX_TH,
                                        facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.3)
        ax.add_patch(rect_duct_b)
        ax.add_patch(rect_duct_t)

    # Dimension Callouts with Arrows
    ax.annotate("", xy=(-R_REFL * 1.04, z_core_bot), xytext=(-R_REFL * 1.04, z_core_top),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.6))
    ax.text(-R_REFL * 1.06, 0.0, f"Active Core Height\nH = {H_CORE:.1f} cm",
            va="center", ha="right", weight="bold", fontsize=10)

    ax.annotate("", xy=(-R_REFL * 1.04, z_tot_bot), xytext=(-R_REFL * 1.04, z_core_bot),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(-R_REFL * 1.06, (z_tot_bot + z_core_bot)/2.0, f"Refl\n{REFL_AX_TH:.1f} cm",
            va="center", ha="right", fontsize=9)

    ax.annotate("", xy=(-R_REFL * 1.04, z_tot_top), xytext=(-R_REFL * 1.04, z_core_top),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(-R_REFL * 1.06, (z_tot_top + z_core_top)/2.0, f"Refl\n{REFL_AX_TH:.1f} cm",
            va="center", ha="right", fontsize=9)

    ax.annotate("", xy=(0, z_tot_top * 1.05), xytext=(R_CORE, z_tot_top * 1.05),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.3))
    ax.text(R_CORE / 2.0, z_tot_top * 1.08, f"R_core = {R_CORE:.2f} cm", ha="center", fontsize=9.5, weight="bold")

    ax.annotate("", xy=(0, z_tot_top * 1.15), xytext=(R_REFL, z_tot_top * 1.15),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.3))
    ax.text(R_REFL / 2.0, z_tot_top * 1.18, f"R_ref = {R_REFL:.2f} cm", ha="center", fontsize=9.5, weight="bold")

    ax.set_xlim(-R_REFL * 1.30, R_REFL * 1.30)
    ax.set_ylim(z_tot_bot * 1.25, z_tot_top * 1.26)
    ax.set_aspect("equal")
    ax.set_xlabel("X Position [cm]", fontsize=12, weight="bold")
    ax.set_ylabel("Z Axial Position [cm]", fontsize=12, weight="bold")
    ax.set_title(f"C01 NEP REACTOR: FULL CORE AXIAL CROSS-SECTION (XZ Plane, y=0)\n"
                 f"Active Core H = {H_CORE:.1f} cm, Total Height H_tot = {z_tot_top - z_tot_bot:.1f} cm, BeO Axial Reflectors = {REFL_AX_TH:.1f} cm",
                 fontsize=12, weight="bold")
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_path}")


# =============================================================================
# 3. UNIT CELL 5-LAYER DETAILED SCHEMATIC WITH AUX COOLANT HOLES
# =============================================================================
def plot_unit_cell_detailed(out_path):
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

    # 1. Hexagonal boundary for unit cell
    hex_r = PITCH / math.sqrt(3.0)
    hex_verts = []
    for k in range(6):
        a = math.radians(k * 60.0 + 30.0)
        hex_verts.append((hex_r * math.cos(a), hex_r * math.sin(a)))
    hex_patch = patches.Polygon(hex_verts, closed=True, facecolor="#f5efe6", edgecolor="#212121",
                                linewidth=2.2, label="YH2.0 Moderator Matrix")
    ax.add_patch(hex_patch)

    # 2. 6 Auxiliary Coolant Channels in YH2
    for k in range(6):
        ang = math.radians(k * 60.0)
        ax_x = R_ADD_POS * math.cos(ang)
        ax_y = R_ADD_POS * math.sin(ang)
        c_ax = plt.Circle((ax_x, ax_y), R_ADD, facecolor="#00bcd4", edgecolor="#006064", linewidth=1.6)
        ax.add_patch(c_ax)
        if k == 0:
            ax.annotate(f"6x Aux Coolant Channels in YH2\nr_aux = {R_ADD*10:.2f} mm (He at 880 K)",
                        xy=(ax_x, ax_y), xytext=(ax_x + 0.35, ax_y + 0.35),
                        arrowprops=dict(arrowstyle="->", color="#006064", lw=1.8),
                        fontsize=10, weight="bold", color="#006064")

    # 3. Layer 4: Al2O3 Liner (100 um)
    c_al = plt.Circle((0, 0), R_AL2O3, facecolor="#ffe0b2", edgecolor="#fb8c00",
                      linewidth=1.8, label="Al2O3 Hydrogen Barrier (100 um)")
    ax.add_patch(c_al)

    # 4. Layer 3: He Annular Hot Coolant Gap (0.84 mm)
    c_chan = plt.Circle((0, 0), R_CHAN, facecolor="#b3e5fc", edgecolor="#0288d1",
                        linewidth=1.5, label="He Annular Hot Channel (0.84 mm)")
    ax.add_patch(c_chan)

    # 5. Layer 2: Nb-1Zr Cladding (0.48 mm)
    c_clad = plt.Circle((0, 0), R_CLAD, facecolor="#cfd8dc", edgecolor="#546e7a",
                        linewidth=1.5, label="Nb-1Zr Cladding (0.48 mm)")
    ax.add_patch(c_clad)

    # 6. Layer 1: UN Fuel Pellet (r = 5.85 mm)
    c_fuel = plt.Circle((0, 0), R_FUEL, facecolor="#8e24aa", edgecolor="#212121",
                        linewidth=2.0, label="UN Fuel Pellet (r = 5.85 mm, 19.75% U-235)")
    ax.add_patch(c_fuel)

    # Center Marker
    ax.plot(0, 0, 'k+', ms=8)

    # Detailed Dimension Callouts with Leader Arrows
    # Fuel radius
    ax.annotate(f"Fuel Pellet Radius\nr_fuel = {R_FUEL*10:.2f} mm",
                xy=(R_FUEL * math.cos(math.radians(135)), R_FUEL * math.sin(math.radians(135))),
                xytext=(-1.15, 0.75),
                arrowprops=dict(arrowstyle="->", color="#8e24aa", lw=1.6),
                fontsize=9.5, weight="bold", color="#8e24aa")

    # Cladding
    ax.annotate(f"Nb-1Zr Cladding\nt_clad = {T_CLAD*10:.2f} mm\nr_clad = {R_CLAD*10:.2f} mm",
                xy=(R_CLAD * math.cos(math.radians(160)), R_CLAD * math.sin(math.radians(160))),
                xytext=(-1.35, 0.15),
                arrowprops=dict(arrowstyle="->", color="#546e7a", lw=1.6),
                fontsize=9.5, weight="bold", color="#37474f")

    # Hot Coolant Gap
    ax.annotate(f"Hot He Annulus (Upflow)\nt_gap = {T_GAP*10:.2f} mm (880 -> 1120 K)",
                xy=(R_CHAN * math.cos(math.radians(205)), R_CHAN * math.sin(math.radians(205))),
                xytext=(-1.35, -0.55),
                arrowprops=dict(arrowstyle="->", color="#0288d1", lw=1.6),
                fontsize=9.5, weight="bold", color="#0288d1")

    # Al2O3 Liner
    ax.annotate(f"Al2O3 Ceramic Liner\nt_liner = {T_AL2O3*1e4:.0f} um\n(H2 Diffusion Barrier)",
                xy=(R_AL2O3 * math.cos(math.radians(245)), R_AL2O3 * math.sin(math.radians(245))),
                xytext=(-1.20, -1.15),
                arrowprops=dict(arrowstyle="->", color="#fb8c00", lw=1.6),
                fontsize=9.5, weight="bold", color="#e65100")

    # Pitch dimension
    ax.annotate("", xy=(-PITCH/2.0, -1.30), xytext=(PITCH/2.0, -1.30),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.8))
    ax.text(0.0, -1.40, f"Hexagonal Lattice Pitch p = {PITCH:.2f} cm", ha="center", fontsize=10, weight="bold")

    ax.set_xlim(-1.45, 1.45)
    ax.set_ylim(-1.50, 1.45)
    ax.set_aspect("equal")
    ax.set_xlabel("X Position [cm]", fontsize=11, weight="bold")
    ax.set_ylabel("Y Position [cm]", fontsize=11, weight="bold")
    ax.set_title("C01 NEP REACTOR: 5-LAYER UNIT CELL ARCHITECTURE\n"
                 "Concentric Dual-Pass Coolant Channels & Al2O3 Permeation Barrier (p = 2.28 cm)",
                 fontsize=12, weight="bold")
    ax.legend(loc="upper left", fontsize=9.2, frameon=True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_path}")


# =============================================================================
# 4. CONTROL DRUM MECHANISM AND ROTATION TRAJECTORY
# =============================================================================
def plot_control_drum_mechanism(out_path):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), dpi=300)

    for ax, theta_deg, title_state in [
        (ax1, 0.0, "DRUMS OUT (0°) — CRITICAL OPERATIONAL STATE"),
        (ax2, 180.0, "DRUMS IN (180°) — SAFETY SHUTDOWN STATE")
    ]:
        # Reflector arc background
        arc_ref = patches.Wedge((0, 0), R_REFL, -40, 40, width=REFL_RAD_TH,
                                facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5, label="BeO Reflector")
        ax.add_patch(arc_ref)

        # Core boundary arc
        arc_core = patches.Wedge((0, 0), R_CORE, -40, 40, facecolor="#f5efe6",
                                 edgecolor="#222222", linewidth=1.8, label="Active Core (YH2 Matrix)")
        ax.add_patch(arc_core)

        # Drum cylinder in reflector
        c_drum = plt.Circle((DRUM_RAD_POS, 0), DRUM_R, facecolor="#b0bec5",
                            edgecolor="#37474f", linewidth=1.6, label="BeO Drum Cylinder (r = 3.95 cm)")
        ax.add_patch(c_drum)

        # B4C Absorber sector rotated by theta_deg
        rot = theta_deg
        half_arc = ABS_ANGLE / 2.0
        th1 = rot - half_arc
        th2 = rot + half_arc

        wedge = patches.Wedge((DRUM_RAD_POS, 0), DRUM_R, th1, th2, width=ABS_TH,
                              facecolor="#212121", edgecolor="red" if theta_deg == 180.0 else "black",
                              linewidth=1.8, label=f"B4C Absorber (155°, t = {ABS_TH:.2f} cm)")
        ax.add_patch(wedge)

        # Drum center axis
        ax.plot(DRUM_RAD_POS, 0, 'ro', ms=6, label="Drum Rotation Pivot")

        # Rotation arrow
        arrow_arc = patches.Arc((DRUM_RAD_POS, 0), 2.2*DRUM_R, 2.2*DRUM_R,
                                angle=0, theta1=20, theta2=160, color="#1976d2", lw=2.0, ls="--")
        ax.add_patch(arrow_arc)

        # Annotations
        if theta_deg == 0.0:
            ax.annotate("B4C Faces Outward\n(Maximum Neutron Reflection\nk_eff = 1.05945)",
                        xy=(DRUM_RAD_POS + DRUM_R - 0.4, 0), xytext=(DRUM_RAD_POS + 1.8, 1.8),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.6),
                        fontsize=10, weight="bold")
            ax.text(DRUM_RAD_POS, -3.2, "θ = 0.0°\nExcess: +5612 pcm", ha="center",
                    fontsize=10.5, weight="bold", color="darkgreen")
        else:
            ax.annotate("B4C Faces Active Core\n(Maximum Neutron Absorption\nk_eff = 0.99551)",
                        xy=(DRUM_RAD_POS - DRUM_R + 0.4, 0), xytext=(DRUM_RAD_POS - 4.5, 2.2),
                        arrowprops=dict(arrowstyle="->", color="red", lw=1.8),
                        fontsize=10, weight="bold", color="darkred")
            ax.text(DRUM_RAD_POS, -3.2, "θ = 180.0°\nSDM: +451 pcm", ha="center",
                    fontsize=10.5, weight="bold", color="darkred")

        ax.set_xlim(12.0, 32.0)
        ax.set_ylim(-7.0, 7.0)
        ax.set_aspect("equal")
        ax.set_xlabel("X Position [cm]", fontsize=11, weight="bold")
        ax.set_ylabel("Y Position [cm]", fontsize=11, weight="bold")
        ax.set_title(title_state, fontsize=11.5, weight="bold")
        ax.grid(True, alpha=0.35)
        ax.legend(loc="lower left", fontsize=8.5, frameon=True)

    fig.suptitle("C01 NEP REACTOR: 10-DRUM CONTROL MECHANISM & ROTATION PHYSICS\n"
                 "Drum Worth W_drum = 6063 pcm | Total Shutdown Margin SDM = +451 pcm",
                 fontsize=13, weight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_path}")


# =============================================================================
# 5. OPENMC CSG HIGH-RESOLUTION PIXEL RAYTRACING (RADIAL & AXIAL)
# =============================================================================
def plot_openmc_csg_slices(spec, out_xy, out_xz):
    """Generate 3000x3000 pixel OpenMC CSG slices with custom high-contrast material palette."""
    import openmc
    sys.path.insert(0, str(REACTOR_DIR.parent.parent / "scripts"))
    from model import build_model

    # Build base model for plotting
    model = build_model(spec, theta_deg=0.0, T_fuel=800.0, T_mod=800.0, T_cool=800.0, T_refl=800.0,
                        particles=100, batches=5, inactive=2)

    # Color map
    # Material IDs:
    # 1: UN_Inner, 2: UN_Outer, 3: Nb1Zr, 4: Helium, 5: Al2O3, 6: YH_2.0, 7: BeO, 8: B4C_Absorber
    colors = {
        model.materials[0]: (229, 57, 53),    # UN Inner: Vivid Red
        model.materials[1]: (142, 36, 170),   # UN Outer: Deep Purple
        model.materials[2]: (120, 144, 156),  # Nb-1Zr Clad: Slate Gray
        model.materials[3]: (179, 229, 252),  # Helium: Light Blue
        model.materials[4]: (255, 183, 77),   # Al2O3: Bright Orange
        model.materials[5]: (245, 239, 230),  # YH2.0: Cream White
        model.materials[6]: (207, 216, 220),  # BeO Reflector: Cool Gray
        model.materials[7]: (33, 33, 33),     # B4C: Charcoal Black
    }

    R_ref = spec['r_ref_cm']
    H_core = spec['height_cm']
    H_tot = H_core + 2.0 * spec['t_refl_ax_cm']

    # 1. Radial CSG Slice (XY)
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    model.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(R_ref * 2.15, R_ref * 2.15),
        pixels=(2400, 2400),
        basis='xy',
        color_by='material',
        colors=colors,
        axes=ax
    )
    ax.set_title("C01 Reactor — OpenMC CSG Radial Slice (XY, z=0)\n"
                 "Exact Constructive Solid Geometry Geometry Raytracing", fontsize=12, weight="bold")
    ax.set_xlabel("X [cm]", fontsize=11, weight="bold")
    ax.set_ylabel("Y [cm]", fontsize=11, weight="bold")
    fig.tight_layout()
    fig.savefig(out_xy, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_xy}")

    # 2. Axial CSG Slice (XZ)
    fig, ax = plt.subplots(figsize=(9, 11), dpi=300)
    model.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(R_ref * 2.15, H_tot * 1.12),
        pixels=(2000, 2500),
        basis='xz',
        color_by='material',
        colors=colors,
        axes=ax
    )
    ax.set_title("C01 Reactor — OpenMC CSG Axial Slice (XZ, y=0)\n"
                 "Exact Constructive Solid Geometry Geometry Raytracing", fontsize=12, weight="bold")
    ax.set_xlabel("X [cm]", fontsize=11, weight="bold")
    ax.set_ylabel("Z [cm]", fontsize=11, weight="bold")
    fig.tight_layout()
    fig.savefig(out_xz, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[+] Saved high-precision plot: {out_xz}")


def main():
    import json
    with open(REACTOR_DIR / "results.json") as f:
        spec = json.load(f)["spec"]

    print("=" * 80)
    print("GENERATING AUTHORITATIVE HIGH-PRECISION TECHNICAL FIGURES FOR C01 REACTOR")
    print("=" * 80)

    # 1. Full Radial Cross Section
    p1 = FIG_DIR / "01_reactor_radial_xy.png"
    plot_radial_xy(p1)

    # 2. Full Axial Cross Section
    p2 = FIG_DIR / "02_reactor_axial_xz.png"
    plot_axial_xz(p2)

    # 3. 5-Layer Unit Cell Detailed
    p3 = FIG_DIR / "03_unit_cell_5layer_detailed.png"
    plot_unit_cell_detailed(p3)

    # 4. Control Drum Mechanism
    p4 = FIG_DIR / "04_control_drum_mechanism.png"
    plot_control_drum_mechanism(p4)

    # 5 & 6. OpenMC CSG High-Res Raytracing
    p5 = FIG_DIR / "05_openmc_csg_slice_xy.png"
    p6 = FIG_DIR / "06_openmc_csg_slice_xz.png"
    plot_openmc_csg_slices(spec, p5, p6)

    # Copy all figures to artifacts directory for user presentation
    for p in [p1, p2, p3, p4, p5, p6]:
        dest = ARTIFACT_DIR / p.name
        import shutil
        shutil.copy(p, dest)
        print(f"[+] Copied to artifact: {dest}")

    print("\nALL HIGH-PRECISION FIGURES SUCCESSFULLY GENERATED AND COPIED!")


if __name__ == "__main__":
    main()
