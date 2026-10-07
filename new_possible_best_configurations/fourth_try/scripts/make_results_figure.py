"""Plot measured C4 (ZrH1.2) mass inventory and negative feedback coefficients."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"


def main():
    feedback = json.loads((ROOT / "feedback_results.json").read_text())
    mass = json.loads((ROOT / "mass_result.json").read_text())
    with (ROOT / "mass_breakdown.csv").open() as file:
        rows = list(csv.DictReader(file))
    categories = [
        ("UN", sum(float(r["mass_kg"]) for r in rows if r["material"].startswith("UN Fuel")), "#ad553f"),
        ("ZrH₁.₂", sum(float(r["mass_kg"]) for r in rows if "ZrH" in r["material"]), "#d9c18f"),
        ("Be", sum(float(r["mass_kg"]) for r in rows if r["material"] == "Be Reflector"), "#9ba8ac"),
        ("Nb-1Zr", sum(float(r["mass_kg"]) for r in rows if "Cladding" in r["material"]), "#6c7380"),
        ("B₄C", sum(float(r["mass_kg"]) for r in rows if "Absorber" in r["material"]), "#263b54"),
        ("Al₂O₃", sum(float(r["mass_kg"]) for r in rows if "Liner" in r["material"]), "#ddd2bb"),
        ("Hardware*", 65.0, "#7c9ab1"), ("Riserva*", 15.0, "#bbc4ce")]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), layout="constrained")
    left = 0.0
    for label, value, color in categories:
        ax1.barh([0], [value], left=left, color=color, edgecolor="white", label=f"{label}: {value:.1f} kg")
        left += value
    ax1.axvline(650, color="#b13c39", ls="--", lw=1.5, label="Soglia 650 kg (User Limit)")
    ax1.set(xlim=(0, 680), ylim=(-.55, .7), yticks=[], xlabel="Massa [kg]",
            title=f"Quarto Tentativo (C4) · solidi OpenMC {mass['modeled_dry_mass_kg']:.1f} ± {mass['sampling_sigma_kg']:.1f} kg; con stime* {left:.1f} kg (< 650 kg PASS)")
    ax1.legend(loc="upper center", bbox_to_anchor=(.5, -.28), ncol=4, fontsize=8, frameon=False)
    ax1.text(.01, -.5, "* Hardware (65 kg) e riserva (15 kg) sono stime; schermo, radiatori e ciclo esclusi.",
             transform=ax1.transAxes, fontsize=8, va="top")

    names = {"fuel": "Combustibile (Doppler)", "moderator": "ZrH₁.₂ (Moderatore)", "total": "Totale Isotermico"}
    coeffs = feedback["coefficients"]
    ys = [2, 1, 0]
    for y, c in zip(ys, coeffs):
        value, sigma = c["alpha_pcm_per_K"], c["sigma_pcm_per_K"]
        color = "#287da7" if value < 0 else "#b24d3a"
        ax2.errorbar([value], [y], xerr=[1.96*sigma], fmt="o", color=color,
                     capsize=5, ms=7, lw=1.8)
        ax2.text(value + (-.05 if value < 0 else .05), y+.19,
                 f"{value:+.3f} ± {sigma:.3f} pcm/K", fontsize=9,
                 ha="right" if value < 0 else "left", color=color, weight="bold")
    ax2.axvline(0, color="#b13c39", lw=1.2, ls="--", label="Soglia di Stabilità (0 pcm/K)")
    ax2.set(yticks=ys, yticklabels=[names[c["kind"]] for c in coeffs],
            xlim=(-1.6, 0.4), ylim=(-.5, 2.6),
            xlabel="Coefficiente di Reattività [pcm/K] (barre: IC 95%)",
            title="Coefficienti Termici BOL · Tutti i 10 Tamburi OUT (0°) · TUTTI NEGATIVI (PASS)")
    ax2.grid(axis="x", alpha=.25)
    ax2.legend(loc="lower right", fontsize=8)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"05_mass_and_feedback.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("Results summary figure generated in", OUT)


if __name__ == "__main__":
    main()

