"""Plot measured C3 mass inventory and fixed-density feedback estimates."""
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
        ("YH₂", sum(float(r["mass_kg"]) for r in rows if "YH2" in r["material"]), "#d9c18f"),
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
    ax1.axvline(600, color="#b13c39", ls="--", lw=1.5, label="Soglia 600 kg")
    ax1.set(xlim=(0, 625), ylim=(-.55, .7), yticks=[], xlabel="Massa [kg]",
            title=f"C3 · solidi OpenMC {mass['modeled_dry_mass_kg']:.1f} ± {mass['sampling_sigma_kg']:.1f} kg; con stime* {left:.1f} kg")
    ax1.legend(loc="upper center", bbox_to_anchor=(.5, -.28), ncol=4, fontsize=8, frameon=False)
    ax1.text(.01, -.5, "* Hardware e riserva sono stime; schermatura, radiatori e conversione esclusi.",
             transform=ax1.transAxes, fontsize=8, va="top")

    names = {"fuel": "Combustibile", "moderator": "YH₂", "total": "Totale"}
    coeffs = feedback["coefficients"]
    ys = [2, 1, 0]
    for y, c in zip(ys, coeffs):
        value, sigma = c["alpha_pcm_per_K"], c["sigma_pcm_per_K"]
        color = "#b24d3a" if value > 0 else "#287da7"
        ax2.errorbar([value], [y], xerr=[1.96*sigma], fmt="o", color=color,
                     capsize=5, ms=7, lw=1.8)
        ax2.text(value + (.08 if value > 0 else -.08), y+.19,
                 f"{value:+.3f} ± {sigma:.3f}", fontsize=9,
                 ha="left" if value > 0 else "right")
    ax2.axvline(0, color="#303c47", lw=1.2)
    ax2.set(yticks=ys, yticklabels=[names[c["kind"]] for c in coeffs],
            xlim=(-1.6, 1.65), ylim=(-.5, 2.6),
            xlabel="Coefficiente [pcm/K] (barre: IC statistico 95%)",
            title="Feedback BOL · tamburi OUT · densità dei solidi fissa")
    ax2.grid(axis="x", alpha=.18)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"05_mass_and_feedback.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)

    coupled = json.loads((ROOT / "coupled_feedback_results.json").read_text())
    fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
    fixed_map = {c["kind"]: c for c in coeffs}
    coupled_map = {c["kind"]: c for c in coupled["coefficients"]}
    for index, kind in enumerate(("fuel", "moderator", "total")):
        c = fixed_map[kind]
        ax.errorbar([c["alpha_pcm_per_K"]], [2-index+.11],
                    xerr=[1.96*c["sigma_pcm_per_K"]], fmt="o", capsize=4,
                    color="#b24d3a", label="Densità dei solidi fissa" if index==0 else None)
        if kind in coupled_map:
            c = coupled_map[kind]
            ax.errorbar([c["alpha_pcm_per_K"]], [2-index-.11],
                        xerr=[1.96*c["sampling_sigma_pcm_per_K"]],
                        fmt="s", capsize=4, color="#247b75",
                        label="Scenario dilatazione YH₂" if index==1 else None)
    ax.axvline(0, color="#303c47", lw=1.1)
    ax.set(yticks=[2, 1, 0], yticklabels=["Combustibile", "YH₂", "Totale"],
           xlim=(-1.55, 1.6), ylim=(-.5, 2.55),
           xlabel="Coefficiente [pcm/K] (barre: IC statistico 95%)",
           title="C3 · temperatura a BOL: effetto dell'ipotesi di dilatazione")
    ax.grid(axis="x", alpha=.18)
    ax.legend(loc="upper right", fontsize=9)
    ax.text(.01, -.23,
            "Scenario: β volumetrico YH₂ = 3,6×10⁻⁵/K; geometria fissa. Il β reale non è qualificato per C3.",
            transform=ax.transAxes, fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"06_feedback_scenarios.{suffix}", dpi=220,
                    bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
