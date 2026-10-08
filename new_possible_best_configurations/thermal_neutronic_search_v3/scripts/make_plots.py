#!/usr/bin/env python3
"""
make_plots.py - High-resolution CSG and Engineering Plot Generator for NEP Reactors
"""

from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import openmc
from model import build_model

def generate_core_figures(spec, out_dir):
    """Generate high-resolution PNG cross sections and schematics."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    H_core = spec['height_cm']
    R_ref = spec['r_ref_cm']

    # 1. Radial Cross Section (Drums OUT - 0 deg)
    m_out = build_model(spec, theta_deg=0.0, T_fuel=800.0, T_mod=800.0, T_cool=800.0, T_refl=800.0, particles=100, batches=5, inactive=2)
    fig, ax = plt.subplots(figsize=(8, 8))
    m_out.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(R_ref * 2.2, R_ref * 2.2),
        basis='xy',
        color_by='material',
        axes=ax
    )
    ax.set_title(f"{spec['id']} — Radial Section (Drums OUT, 0°)", fontsize=13, weight='bold')
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    p_out = out_dir / "01_radial_drums_out.png"
    fig.savefig(p_out, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 2. Radial Cross Section (Drums IN - 180 deg)
    m_in = build_model(spec, theta_deg=180.0, T_fuel=800.0, T_mod=800.0, T_cool=800.0, T_refl=800.0, particles=100, batches=5, inactive=2)
    fig, ax = plt.subplots(figsize=(8, 8))
    m_in.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(R_ref * 2.2, R_ref * 2.2),
        basis='xy',
        color_by='material',
        axes=ax
    )
    ax.set_title(f"{spec['id']} — Radial Section (Drums IN, 180°)", fontsize=13, weight='bold')
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    p_in = out_dir / "02_radial_drums_in.png"
    fig.savefig(p_in, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 3. Axial Cross Section (XZ plane)
    fig, ax = plt.subplots(figsize=(8, 10))
    m_out.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(R_ref * 2.2, (H_core + 2.0 * spec['t_refl_ax_cm']) * 1.15),
        basis='xz',
        color_by='material',
        axes=ax
    )
    ax.set_title(f"{spec['id']} — Axial Section (X-Z)", fontsize=13, weight='bold')
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Z [cm]")
    p_ax = out_dir / "03_axial_section.png"
    fig.savefig(p_ax, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 4. Zoomed Unit Cell
    pitch = spec['pitch_cm']
    fig, ax = plt.subplots(figsize=(7, 7))
    m_out.geometry.root_universe.plot(
        origin=(0.0, 0.0, 0.0),
        width=(pitch * 2.2, pitch * 2.2),
        basis='xy',
        color_by='material',
        axes=ax
    )
    ax.set_title(f"{spec['id']} — Unit Cell (Fuel, Clad, He, Liner, YHx)", fontsize=12, weight='bold')
    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    p_cell = out_dir / "04_unit_cell.png"
    fig.savefig(p_cell, dpi=200, bbox_inches='tight')
    plt.close(fig)

    print(f"Generated 4 figures in {out_dir}")

if __name__ == "__main__":
    from catalog import get_candidates
    c = get_candidates()[0]
    generate_core_figures(c, Path(__file__).resolve().parents[1] / "figures")
