"""BOL criticality check (drums OUT, 0 deg) and analytic mass estimate for fourth_try."""
from __future__ import annotations

import json
import os
from pathlib import Path

import openmc

import build_model as b


ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "campaign" / "c4_bol"
PARAMETERS = dict(pitch=2.30, core_h=47.0, refl_rad_th=8.8,
                  refl_ax_th=8.5, particles=15000, batches=130,
                  inactive=40, n_axial_slices=6, seed=1010,
                  drum_angles=[0.0] * 10)


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    b.PITCH = PARAMETERS["pitch"]
    b.H_CORE = PARAMETERS["core_h"]
    b.REFL_RAD_TH = PARAMETERS["refl_rad_th"]
    b.REFL_AX_TH = PARAMETERS["refl_ax_th"]
    mass = b.compute_mass_breakdown()
    model = b.build_openmc_model_tb35(str(RUN), **PARAMETERS)
    # Use continuous interpolation on ENDF/B-VIII.0 cross sections
    model["settings"].temperature = {"multipole": False, "method": "interpolation",
                                     "tolerance": 10.0}
    model["settings"].export_to_xml(RUN / "settings.xml")
    os.environ["OPENMC_CROSS_SECTIONS"] = b.CROSS_SECTIONS_PATH
    os.environ["OMP_NUM_THREADS"] = "8"
    openmc.run(cwd=str(RUN), threads=8, openmc_exec=b.OPENMC_EXEC, output=False)
    statepoint_path = max(RUN.glob("statepoint.*.h5"), key=lambda p: int(p.name.split(".")[1]))
    with openmc.StatePoint(statepoint_path) as statepoint:
        keff = statepoint.keff
        k_val = float(keff.nominal_value)
        k_sig = float(keff.std_dev)
        rho_pcm = 1e5 * (k_val - 1.0) / k_val
        result = {
            "case": "fourth_try_C4_BOL_all_drums_OUT",
            "k_eff": k_val,
            "k_eff_sigma": k_sig,
            "reactivity_excess_pcm": rho_pcm,
            "confidence_95_normal": [k_val - 1.96 * k_sig, k_val + 1.96 * k_sig],
            "modeled_dry_mass_kg_analytic": mass["M_tot_kg"],
            "mass_breakdown_kg": {
                "UN": mass["M_un_kg"],
                "Nb1Zr": mass["M_clad_kg"],
                "Al2O3": mass["M_al_kg"],
                "ZrH1.2": mass["M_zrh_kg"],
                "Be": mass["M_beo_kg"],
                "B4C": mass["M_b4c_kg"]
            },
            "margin_to_650kg_limit": 650.0 - mass["M_tot_kg"],
            "hardware_plus_reserve_allowance_kg_not_modeled": 80.0,
            "screening_mass_with_allowance_kg": mass["M_tot_kg"] + 80.0,
            "parameters": PARAMETERS,
            "cross_sections": b.CROSS_SECTIONS_PATH,
        }
    (ROOT / "bol_result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

