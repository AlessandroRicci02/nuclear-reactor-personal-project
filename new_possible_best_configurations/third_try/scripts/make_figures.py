"""Draw C3 geometry from the same dimensions passed to OpenMC."""
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
        ax.add_patch(Circle((x, y), .727, color=COLORS["liner"], lw=0))
        ax.add_patch(Circle((x, y), .717, color=COLORS["coolant"], lw=0))
        ax.add_patch(Circle((x, y), .633, color=COLORS["clad"], lw=0))
        ax.add_patch(Circle((x, y), .585, color=COLORS["fuel"], lw=0))
    for j in range(10):
        angle = j*36
        rad = math.radians(angle)
        center = (DRUM_POS*math.cos(rad), DRUM_POS*math.sin(rad))
        ax.add_patch(Circle(center, DRUM_R, facecolor="#aeb8ba", edgecolor="#354151", lw=.8))
        ax.add_patch(Wedge(center, DRUM_R, angle+theta-74, angle+theta+74,
                           width=1.5, color=COLORS["absorber"], lw=0))
        ax.add_patch(Circle(center, DRUM_R, fill=False, edgecolor="#354151", lw=.8))
    ax.set(xlim=(-30, 30), ylim=(-30, 30), aspect="equal",
           xlabel="x [cm]", ylabel="y [cm]",
           title=f"C3 · sezione radiale BOL · tamburi {theta}° ({'OUT' if theta==0 else 'IN'})")
    ax.grid(alpha=.12)
    ax.legend(handles=[Patch(color=COLORS["fuel"], label="UN (217 pin)"),
                       Patch(color=COLORS["moderator"], label="YH₂"),
                       Patch(color=COLORS["reflector"], label="Be"),
                       Patch(color=COLORS["absorber"], label="B₄C")],
              loc="upper right", fontsize=9, framealpha=.95)
    save(fig, f"01_radial_drums_{'out' if theta==0 else 'in'}")


def axial():
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    half = CORE_H/2
    total = half + AX_REF
    ax.add_patch(Rectangle((-REFL_R, -total), 2*REFL_R, 2*total,
                           facecolor=COLORS["reflector"], edgecolor="#354151"))
    ax.add_patch(Rectangle((-CORE_R, -half), 2*CORE_R, CORE_H,
                           facecolor=COLORS["moderator"], edgecolor="#4e564c"))
    for x in [j*PITCH for j in range(-8, 9)]:
        ax.add_patch(Rectangle((x-.717, -half), 1.434, CORE_H,
                               facecolor=COLORS["coolant"], edgecolor="none"))
        ax.add_patch(Rectangle((x-.585, -half), 1.170, CORE_H,
                               facecolor=COLORS["fuel"], edgecolor="none"))
    for x in [-DRUM_POS, DRUM_POS]:
        ax.add_patch(Rectangle((x-DRUM_R, -half), 2*DRUM_R, CORE_H,
                               facecolor="#aeb8ba", edgecolor="#354151", alpha=.75))
    ax.annotate("47,0 cm zona attiva", (0, 0), ha="center", fontsize=11,
                bbox={"facecolor":"white", "alpha":.88, "edgecolor":"none"})
    ax.annotate("8,5 cm Be", (0, half+AX_REF/2), ha="center", va="center", fontsize=10)
    ax.annotate("8,5 cm Be", (0, -half-AX_REF/2), ha="center", va="center", fontsize=10)
    ax.set(xlim=(-31, 31), ylim=(-35, 35), aspect="equal", xlabel="x [cm]",
           ylabel="z [cm]", title="C3 · sezione assiale schematica")
    ax.grid(alpha=.12)
    save(fig, "03_axial_section")


def unit_cell():
    fig, ax = plt.subplots(figsize=(8, 8), layout="constrained")
    circum = PITCH/math.sqrt(3)
    vertices = [(circum*math.cos(math.radians(30+60*j)),
                 circum*math.sin(math.radians(30+60*j))) for j in range(6)]
    ax.add_patch(Polygon(vertices, closed=True, facecolor=COLORS["moderator"],
                         edgecolor="#354151", lw=1.6))
    for radius, name in [(0.727,"liner"),(0.717,"coolant"),
                         (0.633,"clad"),(0.585,"fuel")]:
        ax.add_patch(Circle((0, 0), radius, color=COLORS[name], lw=0))
    aux_r = .910*(PITCH/2.270)
    for j in range(6):
        angle = math.radians(60*j)
        ax.add_patch(Circle((aux_r*math.cos(angle), aux_r*math.sin(angle)),
                            .150, color=COLORS["coolant"], ec="#59899d", lw=.6))
    ax.annotate("UN", (0, 0), ha="center", va="center", color="white", fontsize=13, weight="bold")
    ax.annotate("6 canali He ausiliari", (0, -1.49), ha="center", fontsize=10)
    ax.set(xlim=(-1.65, 1.65), ylim=(-1.7, 1.65), aspect="equal",
           xlabel="x [cm]", ylabel="y [cm]", title="C3 · cella esagonale · passo 2,30 cm")
    ax.grid(alpha=.12)
    ax.legend(handles=[Patch(color=COLORS["fuel"], label="UN"),
                       Patch(color=COLORS["clad"], label="Nb-1Zr"),
                       Patch(color=COLORS["coolant"], label="He"),
                       Patch(color=COLORS["moderator"], label="YH₂")],
              loc="upper right", fontsize=9)
    save(fig, "04_unit_cell")


if __name__ == "__main__":
    radial(0)
    radial(180)
    axial()
    unit_cell()
    print(f"Wrote reactor figures to {OUT}")
