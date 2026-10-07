"""Draw C4 geometry (ZrH1.2 moderator) from the exact dimensions passed to OpenMC."""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle, Wedge, Patch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
PITCH = 2.30
CORE_H = 47.0
CORE_R = 19.45 * (PITCH / 2.270)
REFL_R = CORE_R + 8.8
AX_REF = 8.5
DRUM_R = 3.95
DRUM_POS = CORE_R + 0.180 + DRUM_R
COLORS = dict(fuel="#ad553f", clad="#6c7380", coolant="#8ac6dd",
              liner="#ede8dc", moderator="#d9c18f", reflector="#d7dce0",
              absorber="#263b54")


def pins():
    v1 = (PITCH, 0.0)
    v2 = (PITCH / 2, PITCH * math.sqrt(3) / 2)
    points = []
    for ring in range(8, -1, -1):
        if ring == 0:
            points.append((0.0, 0.0))
            continue
        corners = [(ring*v2[0], ring*v2[1]),
                   (ring*(v2[0]-v1[0]), ring*(v2[1]-v1[1])),
                   (-ring*v1[0], -ring*v1[1]),
                   (-ring*v2[0], -ring*v2[1]),
                   (-ring*(v2[0]-v1[0]), -ring*(v2[1]-v1[1])),
                   (ring*v1[0], ring*v1[1])]
        for side in range(6):
            a, b = corners[side], corners[(side+1) % 6]
            for step in range(ring):
                f = step/ring
                points.append(((1-f)*a[0]+f*b[0], (1-f)*a[1]+f*b[1]))
    assert len(points) == 217
    return points


def save(fig, name):
    fig.savefig(OUT / (name + ".png"), dpi=220, bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def radial(theta):
    fig, ax = plt.subplots(figsize=(9, 9), layout="constrained")
    ax.add_patch(Circle((0, 0), REFL_R, color=COLORS["reflector"], ec="#354151", lw=1.5))
    ax.add_patch(Circle((0, 0), CORE_R, color=COLORS["moderator"], ec="#4e564c", lw=1.2))
    for x, y in pins():
        ax.add_patch(Circle((x, y), 0.585 + 0.048 + 0.084 + 0.010, color=COLORS["liner"]))
        ax.add_patch(Circle((x, y), 0.585 + 0.048 + 0.084, color=COLORS["coolant"]))
        ax.add_patch(Circle((x, y), 0.585 + 0.048, color=COLORS["clad"]))
        ax.add_patch(Circle((x, y), 0.585, color=COLORS["fuel"]))
        r_aux_pos = 0.910 * (PITCH / 2.270)
        for k in range(6):
            phi = k * math.pi / 3.0
            ax.add_patch(Circle((x + r_aux_pos*math.cos(phi), y + r_aux_pos*math.sin(phi)),
                                0.150, color=COLORS["coolant"], ec="#354151", lw=0.3))
    for k in range(10):
        phi_deg = k * 36.0
        phi_rad = math.radians(phi_deg)
        cx = DRUM_POS * math.cos(phi_rad)
        cy = DRUM_POS * math.sin(phi_rad)
        ax.add_patch(Circle((cx, cy), DRUM_R, color=COLORS["reflector"], ec="#222b35", lw=1.0))
        center_deg = phi_deg + theta
        ax.add_patch(Wedge((cx, cy), DRUM_R, center_deg - 74.0, center_deg + 74.0,
                           width=1.50, color=COLORS["absorber"], ec="#0b1723", lw=0.8))
    ax.set(aspect="equal", xlim=(-REFL_R-2, REFL_R+2), ylim=(-REFL_R-2, REFL_R+2),
           title=f"Quarto Tentativo (C4 · ZrH1.2) · Sezione Radiale · Tamburi {int(theta)}° ({'OUT' if theta==0 else 'IN'})",
           xlabel="x [cm]", ylabel="y [cm]")
    handles = [Patch(facecolor=COLORS["fuel"], label="UN (HALEU)"),
               Patch(facecolor=COLORS["moderator"], label="ZrH1.2 Moderatore"),
               Patch(facecolor=COLORS["coolant"], label="Elio refrigerante (2.5 MPa)"),
               Patch(facecolor=COLORS["reflector"], label="Riflettore Be metallico"),
               Patch(facecolor=COLORS["absorber"], label="B4C (90% B-10)")]
    ax.legend(handles=handles, loc="upper right", framealpha=0.9, fontsize=8)
    save(fig, f"01_radial_drums_{'out' if theta == 0 else 'in'}")


def axial():
    fig, ax = plt.subplots(figsize=(7, 9), layout="constrained")
    ax.add_patch(Rectangle((-REFL_R, -CORE_H/2 - AX_REF), 2*REFL_R, CORE_H + 2*AX_REF,
                           color=COLORS["reflector"], ec="#222b35", lw=1.2))
    ax.add_patch(Rectangle((-CORE_R, -CORE_H/2), 2*CORE_R, CORE_H,
                           color=COLORS["moderator"], ec="#4e564c", lw=1.0))
    for x, _ in pins()[:21]:
        ax.add_patch(Rectangle((x - 0.585, -CORE_H/2), 2*0.585, CORE_H,
                               color=COLORS["fuel"], ec="none"))
    for sign in (-1, 1):
        cx = sign * DRUM_POS
        ax.add_patch(Rectangle((cx - DRUM_R, -CORE_H/2), 2*DRUM_R, CORE_H,
                               color="#cad2d7", ec="#222b35", lw=0.8))
        ax.add_patch(Rectangle((cx + (DRUM_R-1.5 if sign>0 else -DRUM_R), -CORE_H/2),
                               1.5, CORE_H, color=COLORS["absorber"], ec="#0b1723", lw=0.8))
    ax.set(aspect="equal", xlim=(-REFL_R-3, REFL_R+3), ylim=(-CORE_H/2 - AX_REF - 3, CORE_H/2 + AX_REF + 3),
           title="Quarto Tentativo (C4 · ZrH1.2) · Sezione Assiale (x-z)",
           xlabel="x [cm]", ylabel="z [cm]")
    save(fig, "03_axial_section")


def unit_cell():
    fig, ax = plt.subplots(figsize=(6, 6), layout="constrained")
    corners = [(PITCH/math.sqrt(3)*math.cos(math.radians(30 + 60*i)),
                PITCH/math.sqrt(3)*math.sin(math.radians(30 + 60*i))) for i in range(6)]
    ax.add_patch(Polygon(corners, color=COLORS["moderator"], ec="#354151", lw=1.5))
    ax.add_patch(Circle((0, 0), 0.585 + 0.048 + 0.084 + 0.010, color=COLORS["liner"], ec="#7c848f"))
    ax.add_patch(Circle((0, 0), 0.585 + 0.048 + 0.084, color=COLORS["coolant"], ec="#5a889b"))
    ax.add_patch(Circle((0, 0), 0.585 + 0.048, color=COLORS["clad"], ec="#39414d"))
    ax.add_patch(Circle((0, 0), 0.585, color=COLORS["fuel"], ec="#562118"))
    r_aux_pos = 0.910 * (PITCH / 2.270)
    for k in range(6):
        phi = k * math.pi / 3.0
        ax.add_patch(Circle((r_aux_pos*math.cos(phi), r_aux_pos*math.sin(phi)),
                            0.150, color=COLORS["coolant"], ec="#354151", lw=0.6))
    ax.set(aspect="equal", xlim=(-1.5, 1.5), ylim=(-1.5, 1.5),
           title="Cella Unitaria a 5 Strati + 6 Canali Ausiliari (ZrH1.2)",
           xlabel="x [cm]", ylabel="y [cm]")
    save(fig, "04_unit_cell")


def main():
    radial(0)
    radial(180)
    axial()
    unit_cell()
    print("Geometry figures generated in", OUT)


if __name__ == "__main__":
    main()

