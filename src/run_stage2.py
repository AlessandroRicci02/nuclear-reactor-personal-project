#!/usr/bin/env python3
"""
Executes Stage 2 intermediate evaluation on the 4 top architectural finalists.
"""
import os
import sys

from screen_architectures import ArchitectureCandidate
from evaluate_stage2_finalists import run_stage2_evaluation

def main():
    finalists = []

    # Finalist A-Star: 10-Drum Advanced Pin Core (Family A)
    finalists.append(ArchitectureCandidate(
        cand_id="A*",
        family="Family A (Square Multi-Drum)",
        name="10-Drum Optimized Pin Core",
        lattice_type='square',
        core_radius=19.50,
        active_core_height=42.00,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.95,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # Finalist B-Star: 10-Drum Hexagonal Lattice Core (Family B)
    finalists.append(ArchitectureCandidate(
        cand_id="B*",
        family="Family B (Hexagonal Lattice)",
        name="10-Drum Hexagonal Optimized Core",
        lattice_type='hex',
        core_radius=19.35,
        active_core_height=42.00,
        pitch=1.810,
        number_of_drums=10,
        drum_radius=3.95,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # Finalist C-Star: BeO-Island Annular Core (Family C)
    finalists.append(ArchitectureCandidate(
        cand_id="C*",
        family="Family C (Annular Island)",
        name="BeO-Island Annular Core",
        lattice_type='annular',
        core_radius=20.20,
        active_core_height=41.50,
        pitch=1.665,
        inner_island_radius=4.50,
        island_material='BeO',
        number_of_drums=10,
        drum_radius=4.00,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.40,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.820
    ))

    # Finalist D-Star: 12-Drum Modular Enhanced-Streaming Core (Family D / Hybrid)
    finalists.append(ArchitectureCandidate(
        cand_id="D*",
        family="Family D (Clustered/Zoned)",
        name="12-Drum Enhanced Streaming Modular Core",
        lattice_type='square',
        core_radius=19.50,
        active_core_height=42.00,
        pitch=1.665,
        number_of_drums=12,
        drum_radius=3.55,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.880
    ))

    run_stage2_evaluation(finalists, output_dir="results/stage2_evaluation")

if __name__ == "__main__":
    main()
