"""Reproduce the third_try BOL criticality check and analytic mass estimate."""
from __future__ import annotations

import json
import os
from pathlib import Path

import openmc

import build_model as b


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "campaign" / "c3_bol"
PARAMETERS = dict(pitch=2.30, core_h=47.0, refl_rad_th=8.8,
                  refl_ax_th=8.5, particles=15000, batches=170,
                  inactive=50, n_axial_slices=6, seed=1010,
                  drum_angles=[0.0] * 10)


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    b.PITCH = PARAMETERS["pitch"]
    b.H_CORE = PARAMETERS["core_h"]
    b.REFL_RAD_TH = PARAMETERS["refl_rad_th"]
    b.REFL_AX_TH = PARAMETERS["refl_ax_th"]
    mass = b.compute_mass_breakdown()
    model = b.build_openmc_model_tb35(str(RUN), **PARAMETERS)
    # Preserve the temperature treatment used for the archived BOL statepoint.
    model["settings"].temperature = {"multipole": True, "method": "nearest",
                                     "tolerance": 1500.0}
    model["settings"].export_to_xml(RUN / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    openmc.run(cwd=str(RUN), threads=8, openmc_exec=b.OPENMC_EXEC, output=False)
    with openmc.StatePoint(RUN / "statepoint.170.h5") as statepoint:
        keff = statepoint.keff
        result = {
            "case": "third_try_C3_BOL_all_drums_OUT",
            "k_eff": float(keff.nominal_value),
            "k_eff_sigma": float(keff.std_dev),
            "confidence_95_normal": [float(keff.nominal_value - 1.96 * keff.std_dev),
                                     float(keff.nominal_value + 1.96 * keff.std_dev)],
            "modeled_dry_mass_kg_analytic": mass["M_tot_kg"],
            "mass_breakdown_kg": {"UN": mass["M_un_kg"], "Nb1Zr": mass["M_clad_kg"],
                                  "Al2O3": mass["M_al_kg"], "YH2": mass["M_yh2_kg"],
                                  "Be": mass["M_beo_kg"], "B4C": mass["M_b4c_kg"]},
            "hardware_plus_reserve_allowance_kg_not_modeled": 80.0,
            "screening_mass_with_allowance_kg": mass["M_tot_kg"] + 80.0,
            "parameters": PARAMETERS,
            "cross_sections": b.CROSS_SECTIONS_PATH,
        }
    (ROOT / "bol_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
