#!/usr/bin/env python3
"""
generate_reactor_plots_tb35.py

Generates publication-quality technical CAD/CSG figures (300 DPI) for GEOM_OPT15_TB35:
  1. geom_tb35_01_reactor_radial_xy.png: Full radial cross section (217 pins, 6 aux channels, 10 drums, BeO reflector)
  2. geom_tb35_02_reactor_axial_xz.png: Full axial cross section (active core, BeO axial reflectors, plenums)
  3. geom_tb35_03_unit_cell_5layer_detailed.png: Close-up 5-layer unit cell with 6 auxiliary coolant holes
  4. geom_tb35_04_control_drum_mechanism.png: Control drum mechanism and rotation trajectory
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/best_configurations/GEOM_OPT15_TB35"
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "figure.titlesize": 14
})

# Geometric constants
PITCH = 2.280       # cm
H_CORE = 45.50      # cm
REFL_AX_TH = 8.80   # cm
REFL_RAD_TH = 9.30  # cm
R_FUEL = 0.585      # cm
T_CLAD = 0.048      # cm
R_CLAD = R_FUEL + T_CLAD # 0.633 cm
T_GAP = 0.084       # cm
R_CHAN = R_CLAD + T_GAP  # 0.717 cm
T_AL2O3 = 0.010     # cm (100 um)
R_AL2O3 = R_CHAN + T_AL2O3 # 0.727 cm
R_ADD = 0.150       # cm (1.50 mm)
R_ADD_POS = 0.910 * (PITCH / 2.270) # 0.914 cm
R_CORE = 19.45 * (PITCH / 2.270)    # 19.536 cm
R_REFL = R_CORE + REFL_RAD_TH       # 28.836 cm
DRUM_R = 3.95       # cm
ABS_TH = 1.50       # cm
ABS_ANGLE = 148.0   # deg
DRUM_CLEARANCE = 0.180 # cm
DRUM_COUNT = 10
DRUM_RAD_POS = R_CORE + DRUM_CLEARANCE + DRUM_R # 23.666 cm


# -----------------------------------------------------------------------------
# 1. RADIAL XY SLICE
# -----------------------------------------------------------------------------
def plot_radial_xy():
    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

    # 1. Radial Reflector Outer Circle
    c_outer = plt.Circle((0, 0), R_REFL, facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5, label="BeO Radial Reflector")
    ax.add_patch(c_outer)

    # 2. Active Core Barrel Circle
    c_core = plt.Circle((0, 0), R_CORE, facecolor="#f5efe6", edgecolor="#222222", linewidth=1.8, label="Active Core Barrel (YH2 Matrix)")
    ax.add_patch(c_core)

    # 3. Hexagonal Lattice Pin Generation
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
        # Annular Coolant Gap
        c_cool = plt.Circle((x, y), R_CHAN, facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.4)
        ax.add_patch(c_cool)
        # Cladding
        c_cl = plt.Circle((x, y), R_CLAD, facecolor="#cfd8dc", edgecolor="#546e7a", linewidth=0.4)
        ax.add_patch(c_cl)
        # Fuel Pellet
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

    # 4. Control Drums in Radial Reflector
    for k in range(DRUM_COUNT):
        phi_k = k * 2.0 * math.pi / DRUM_COUNT
        xk = DRUM_RAD_POS * math.cos(phi_k)
        yk = DRUM_RAD_POS * math.sin(phi_k)

        # Drum BeO body
        c_drum = plt.Circle((xk, yk), DRUM_R, facecolor="#b0bec5", edgecolor="#37474f", linewidth=1.2)
        ax.add_patch(c_drum)

        # B4C Absorber sector (rotated at theta = 80.4 deg)
        rot_deg = math.degrees(phi_k) + 80.40
        half_ang = ABS_ANGLE / 2.0
        theta1 = rot_deg - half_ang
        theta2 = rot_deg + half_ang

        wedge = patches.Wedge((xk, yk), DRUM_R, theta1, theta2, width=ABS_TH,
                              facecolor="#263238", edgecolor="black", linewidth=1.0)
        ax.add_patch(wedge)

    # Annotations and legend elements
    ax.set_xlim(-R_REFL * 1.08, R_REFL * 1.08)
    ax.set_ylim(-R_REFL * 1.08, R_REFL * 1.08)
    ax.set_aspect("equal")
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Y Position [cm]", fontsize=12)
    ax.set_title("GEOM_OPT15_TB35: Full Reactor Radial Cross-Section (XY Midplane, z=0)\n"
                 "217 Pins (9 Rings), 6 Aux Channels/Cell (r=1.50mm), 10 Rotating B4C Drums in BeO Reflector",
                 fontsize=12, fontweight="bold")

    # Custom Legend
    legend_elements = [
        patches.Patch(facecolor="#e53935", edgecolor="black", label="UN Fuel Inner (19.60 wt% U-235, 19 pins)"),
        patches.Patch(facecolor="#8e24aa", edgecolor="black", label="UN Fuel Outer (19.75 wt% U-235, 198 pins)"),
        patches.Patch(facecolor="#cfd8dc", edgecolor="#546e7a", label="Nb-1Zr Cladding (0.48 mm)"),
        patches.Patch(facecolor="#b3e5fc", edgecolor="#0288d1", label="He Coolant Annulus (0.84 mm)"),
        patches.Patch(facecolor="#ffe0b2", edgecolor="#fb8c00", label="Al2O3 Ceramic Liner (100 um)"),
        patches.Patch(facecolor="#00bcd4", edgecolor="#006064", label="6x Aux Coolant Channels (r=1.50 mm)"),
        patches.Patch(facecolor="#f5efe6", edgecolor="#222222", label="Canned YH2 Moderator Matrix (102.6 kg)"),
        patches.Patch(facecolor="#e0e0e0", edgecolor="#333333", label="BeO Radial Reflector (9.30 cm)"),
        patches.Patch(facecolor="#263238", edgecolor="black", label="B4C Absorber Sector (148 deg, 1.5 cm)"),
        patches.Patch(facecolor="#b0bec5", edgecolor="#37474f", label="BeO Control Drum Cylinder (r=3.95 cm)")
    ]
    ax.legend(handles=legend_elements, loc="upper right", bbox_to_anchor=(1.35, 1.0), fontsize=9, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_01_reactor_radial_xy.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[+] Saved: {out_path}")


# -----------------------------------------------------------------------------
# 2. AXIAL XZ SLICE
# -----------------------------------------------------------------------------
def plot_axial_xz():
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    z_core_bot = -H_CORE / 2.0
    z_core_top = +H_CORE / 2.0
    z_tot_bot = z_core_bot - REFL_AX_TH
    z_tot_top = z_core_top + REFL_AX_TH

    # 1. BeO Radial Reflector Blocks (Left & Right)
    rect_refl_l = patches.Rectangle((-R_REFL, z_tot_bot), R_REFL - R_CORE, z_tot_top - z_tot_bot,
                                    facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5, label="BeO Radial Reflector")
    rect_refl_r = patches.Rectangle((+R_CORE, z_tot_bot), R_REFL - R_CORE, z_tot_top - z_tot_bot,
                                    facecolor="#e0e0e0", edgecolor="#333333", linewidth=1.5)
    ax.add_patch(rect_refl_l)
    ax.add_patch(rect_refl_r)

    # 2. Control Drums (spanning core height in reflector)
    for x_d in [-DRUM_RAD_POS, DRUM_RAD_POS]:
        rect_drum = patches.Rectangle((x_d - DRUM_R, z_core_bot), 2 * DRUM_R, H_CORE,
                                      facecolor="#b0bec5", edgecolor="#37474f", linewidth=1.2)
        ax.add_patch(rect_drum)
        # Absorber strip
        sign = 1 if x_d > 0 else -1
        rect_abs = patches.Rectangle((x_d + sign * (DRUM_R - ABS_TH), z_core_bot), sign * (-ABS_TH), H_CORE,
                                     facecolor="#263238", edgecolor="black", linewidth=1.0)
        ax.add_patch(rect_abs)

    # 3. Top and Bottom BeO Axial Reflectors
    rect_ax_bot = patches.Rectangle((-R_CORE, z_tot_bot), 2 * R_CORE, REFL_AX_TH,
                                    facecolor="#d7ccc8", edgecolor="#5d4037", linewidth=1.2, label="Bottom BeO Axial Reflector (8.80 cm)")
    rect_ax_top = patches.Rectangle((-R_CORE, z_core_top), 2 * R_CORE, REFL_AX_TH,
                                    facecolor="#d7ccc8", edgecolor="#5d4037", linewidth=1.2, label="Top BeO Axial Reflector (8.80 cm)")
    ax.add_patch(rect_ax_bot)
    ax.add_patch(rect_ax_top)

    # 4. Active Core YH2 Matrix
    rect_core = patches.Rectangle((-R_CORE, z_core_bot), 2 * R_CORE, H_CORE,
                                  facecolor="#f5efe6", edgecolor="#222222", linewidth=1.5, label="Active Core (YH2 Matrix)")
    ax.add_patch(rect_core)

    # 5. Representative Fuel Pins along X axis
    pin_xs = np.linspace(-R_CORE * 0.88, R_CORE * 0.88, 17)
    for px in pin_xs:
        # Cladding
        rect_clad = patches.Rectangle((px - R_CLAD, z_core_bot), 2 * R_CLAD, H_CORE,
                                      facecolor="#cfd8dc", edgecolor="#546e7a", linewidth=0.4)
        ax.add_patch(rect_clad)
        # Fuel
        rect_fuel = patches.Rectangle((px - R_FUEL, z_core_bot), 2 * R_FUEL, H_CORE,
                                      facecolor="#8e24aa" if abs(px) > 6.0 else "#e53935",
                                      edgecolor="#212121", linewidth=0.5)
        ax.add_patch(rect_fuel)
        # Axial coolant ducts through reflectors
        rect_duct_b = patches.Rectangle((px - 0.43, z_tot_bot), 0.86, REFL_AX_TH,
                                        facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.3)
        rect_duct_t = patches.Rectangle((px - 0.43, z_core_top), 0.86, REFL_AX_TH,
                                        facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=0.3)
        ax.add_patch(rect_duct_b)
        ax.add_patch(rect_duct_t)

    # Dimension Callouts
    ax.annotate("", xy=(-R_REFL * 1.04, z_core_bot), xytext=(-R_REFL * 1.04, z_core_top),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.5))
    ax.text(-R_REFL * 1.06, 0.0, "Active Core Height\nH = 45.50 cm", va="center", ha="right", fontweight="bold", fontsize=9.5)

    ax.annotate("", xy=(-R_REFL * 1.04, z_tot_bot), xytext=(-R_REFL * 1.04, z_core_bot),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(-R_REFL * 1.06, (z_tot_bot + z_core_bot)/2.0, "Refl\n8.80 cm", va="center", ha="right", fontsize=8.5)

    ax.annotate("", xy=(0, z_tot_top * 1.05), xytext=(R_CORE, z_tot_top * 1.05),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(R_CORE / 2.0, z_tot_top * 1.08, "R_core = 19.54 cm", ha="center", fontsize=9)

    ax.annotate("", xy=(0, z_tot_top * 1.15), xytext=(R_REFL, z_tot_top * 1.15),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text(R_REFL / 2.0, z_tot_top * 1.18, "R_ref = 28.84 cm", ha="center", fontsize=9)

    ax.set_xlim(-R_REFL * 1.25, R_REFL * 1.25)
    ax.set_ylim(z_tot_bot * 1.22, z_tot_top * 1.25)
    ax.set_aspect("equal")
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Z Axial Position [cm]", fontsize=12)
    ax.set_title("GEOM_OPT15_TB35: Full Core Axial Cross-Section (XZ Plane, y=0)\n"
                 "Total Height H_tot = 63.10 cm, Active Core H = 45.50 cm, Axial BeO Reflectors = 8.80 cm",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_02_reactor_axial_xz.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[+] Saved: {out_path}")


# -----------------------------------------------------------------------------
# 3. UNIT CELL 5-LAYER DETAILED SCHEMATIC
# -----------------------------------------------------------------------------
def plot_unit_cell():
    fig, ax = plt.subplots(figsize=(9, 9), dpi=300)

    # Hexagonal boundary
    hex_r = PITCH / math.sqrt(3.0)
    hex_verts = []
    for k in range(6):
        a = math.radians(k * 60.0 + 30.0)
        hex_verts.append((hex_r * math.cos(a), hex_r * math.sin(a)))
    hex_patch = patches.Polygon(hex_verts, closed=True, facecolor="#f5efe6", edgecolor="#212121", linewidth=2.0, label="Canned YH2 Moderator Matrix")
    ax.add_patch(hex_patch)

    # 6 Auxiliary Coolant Channels
    for k in range(6):
        ang = math.radians(k * 60.0)
        ax_x = R_ADD_POS * math.cos(ang)
        ax_y = R_ADD_POS * math.sin(ang)
        c_ax = plt.Circle((ax_x, ax_y), R_ADD, facecolor="#00bcd4", edgecolor="#006064", linewidth=1.5)
        ax.add_patch(c_ax)
        if k == 0:
            ax.annotate("6x Aux Coolant Channels\nr_add = 1.50 mm (Orificed Flow)",
                        xy=(ax_x, ax_y), xytext=(ax_x + 0.35, ax_y + 0.35),
                        arrowprops=dict(arrowstyle="->", color="#006064", lw=1.5),
                        fontsize=9.5, fontweight="bold", color="#006064")

    # Layer 4: Al2O3 Liner (100 um)
    c_al = plt.Circle((0, 0), R_AL2O3, facecolor="#ffe0b2", edgecolor="#fb8c00", linewidth=1.5, label="Al2O3 Liner (t = 100 um)")
    ax.add_patch(c_al)

    # Layer 3: He Annular Coolant Gap (0.84 mm)
    c_chan = plt.Circle((0, 0), R_CHAN, facecolor="#b3e5fc", edgecolor="#0288d1", linewidth=1.2, label="He Annular Coolant (t = 0.84 mm)")
    ax.add_patch(c_chan)

    # Layer 2: Nb-1Zr Cladding (0.48 mm)
    c_clad = plt.Circle((0, 0), R_CLAD, facecolor="#cfd8dc", edgecolor="#546e7a", linewidth=1.2, label="Nb-1Zr Clad (t = 0.48 mm)")
    ax.add_patch(c_clad)

    # Layer 1: UN Fuel Pellet (r = 5.85 mm)
    c_fuel = plt.Circle((0, 0), R_FUEL, facecolor="#8e24aa", edgecolor="#212121", linewidth=1.8, label="UN Fuel Pellet (r = 5.85 mm)")
    ax.add_patch(c_fuel)

    # Callouts
    ax.annotate("UN Fuel Pellet (19.75 wt%)\nr = 5.85 mm, T_center = 1255 K", xy=(0, 0), xytext=(-1.1, -1.1),
                arrowprops=dict(arrowstyle="->", color="black", lw=1.2), fontsize=9, fontweight="bold")
    ax.annotate("Nb-1Zr Clad (t = 0.48 mm)\nT_peak = 1242 K (< 1350 K)", xy=(0, R_CLAD), xytext=(-1.2, 0.8),
                arrowprops=dict(arrowstyle="->", color="#546e7a", lw=1.2), fontsize=9)
    ax.annotate("Al2O3 Liner (100 um)\nThermal & Permeation Barrier", xy=(R_AL2O3 * 0.7, R_AL2O3 * 0.7), xytext=(0.9, 0.9),
                arrowprops=dict(arrowstyle="->", color="#fb8c00", lw=1.2), fontsize=9)
    ax.annotate("Canned YH2 Matrix (Pitch = 2.28 cm)\nT_peak = 914.8 K (Buffer = +35.2 K)", xy=(0, hex_r * 0.85), xytext=(-1.0, 1.25),
                arrowprops=dict(arrowstyle="->", color="#5d4037", lw=1.2), fontsize=9.5, fontweight="bold", color="#3e2723")

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X [cm]", fontsize=11)
    ax.set_ylabel("Y [cm]", fontsize=11)
    ax.set_title("GEOM_OPT15_TB35: 5-Layer Insulated Unit Cell Architecture\n"
                 "Decoupled Fuel (1255 K) & Moderator (914.8 K) Achieving +35.2 K Thermal Buffer",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_03_unit_cell_5layer_detailed.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[+] Saved: {out_path}")


# -----------------------------------------------------------------------------
# 4. CONTROL DRUM MECHANISM
# -----------------------------------------------------------------------------
def plot_control_drum_mechanism():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6), dpi=300)

    for ax, rot, title, desc in [
        (ax1, 0.0, "Maximum Excess Reactivity (theta = 0 deg)", "Absorber Turned Outward (Facing Space)"),
        (ax2, 180.0, "Full Core Shutdown (theta = 180 deg)", "Absorber Turned Inward (Facing Core Barrel)")
    ]:
        # Core edge
        arc_core = patches.Arc((-DRUM_RAD_POS, 0), 2*R_CORE, 2*R_CORE, angle=0, theta1=-40, theta2=40,
                               color="#222222", linewidth=2.5, linestyle="--")
        ax.add_patch(arc_core)
        ax.text(-DRUM_RAD_POS + R_CORE + 0.3, 0.0, "Core\nBarrel", va="center", ha="left", fontsize=9, fontweight="bold")

        # Reflector boundary
        arc_refl = patches.Arc((-DRUM_RAD_POS, 0), 2*R_REFL, 2*R_REFL, angle=0, theta1=-40, theta2=40,
                               color="#333333", linewidth=1.5)
        ax.add_patch(arc_refl)

        # Drum BeO cylinder
        c_drum = plt.Circle((0, 0), DRUM_R, facecolor="#b0bec5", edgecolor="#37474f", linewidth=1.8, label="BeO Drum Cylinder (r = 3.95 cm)")
        ax.add_patch(c_drum)

        # B4C Absorber sector
        half_ang = ABS_ANGLE / 2.0
        th1 = rot - half_ang
        th2 = rot + half_ang
        w = patches.Wedge((0, 0), DRUM_R, th1, th2, width=ABS_TH,
                          facecolor="#263238", edgecolor="black", linewidth=1.5, label="B4C Absorber (148 deg, 1.5 cm)")
        ax.add_patch(w)

        # Rotation arrow
        rot_arrow = patches.Arc((0, 0), DRUM_R * 0.7, DRUM_R * 0.7, angle=0, theta1=0, theta2=rot if rot > 0 else 60,
                                color="crimson", linewidth=2.0, linestyle=":")
        ax.add_patch(rot_arrow)

        ax.plot(0, 0, 'ko', markersize=6)
        ax.set_xlim(-5.5, 5.5)
        ax.set_ylim(-5.5, 5.5)
        ax.set_aspect("equal")
        ax.set_xlabel("X [cm]", fontsize=11)
        ax.set_ylabel("Y [cm]", fontsize=11)
        ax.set_title(f"{title}\n{desc}", fontsize=11, fontweight="bold")
        ax.legend(loc="lower left", fontsize=8.5, frameon=True)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_04_control_drum_mechanism.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[+] Saved: {out_path}")


if __name__ == "__main__":
    print("Generating technical reactor images for GEOM_OPT15_TB35...")
    plot_radial_xy()
    plot_axial_xz()
    plot_unit_cell()
    plot_control_drum_mechanism()
    print("All reactor images generated successfully.")
