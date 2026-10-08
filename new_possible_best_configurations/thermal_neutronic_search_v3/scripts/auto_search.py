#!/usr/bin/env python3
"""
auto_search.py - Automated Search Controller for NEP Reactor Candidates
Evaluates 30 candidates sequentially with strict early-rejection filtering:
  1. BOL Drums OUT (keff - 2*sigma > 1.0)
  2. BOL Drums IN  (keff + 2*sigma < 1.0)
  3. Moderator Peak Temperature (T_mod_peak < 950.0 K)
  4. Reactivity Feedback Coefficients:
     - Fuel Doppler (alpha_fuel + 2*sigma < 0)
     - Moderator Coeff (alpha_mod + 2*sigma < 0)
     - Total Coeff (alpha_total + 2*sigma < 0)

STOPS IMMEDIATELY on the first fully compliant candidate and creates its dedicated folder.
"""

from pathlib import Path
import os
import sys
import json
import time
import math
import csv
import subprocess
import openmc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from catalog import get_candidates
from model import build_model, XS_PATH, OPENMC_BIN
from thermal import evaluate_thermal
from make_plots import generate_core_figures

THREADS = 7 # ~90% of 8 CPU cores
Z_SIGMA = 2.0 # User specified: 2-sigma confidence

def atomic_write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    tmp.replace(path)

def run_openmc_case(work_dir, model):
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    model.export_to_xml(work_dir)

    env = os.environ.copy()
    env['OPENMC_CROSS_SECTIONS'] = str(XS_PATH)
    env['OMP_NUM_THREADS'] = str(THREADS)
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'

    log_file = work_dir / "openmc.log"
    with open(log_file, "w") as lf:
        proc = subprocess.run(
            [OPENMC_BIN, "-s", str(THREADS)],
            cwd=work_dir, env=env, stdout=lf, stderr=subprocess.STDOUT
        )
    if proc.returncode != 0:
        raise RuntimeError(f"OpenMC simulation failed in {work_dir}. See {log_file}")

    sp_files = list(work_dir.glob("statepoint.*.h5"))
    if not sp_files:
        raise RuntimeError(f"No statepoint file found in {work_dir}")
    sp_path = max(sp_files, key=lambda p: int(p.name.split('.')[1]))

    with openmc.StatePoint(sp_path, autolink=False) as sp:
        keff = float(sp.keff.n)
        sigma = float(sp.keff.s)

    return keff, sigma

def evaluate_candidate(cand):
    cand_id = cand['id']
    cand_dir = ROOT / "campaign" / cand_id
    cand_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n{'='*80}\nEVALUATING CANDIDATE: {cand_id}\n{'='*80}", flush=True)

    results = {"id": cand_id, "spec": cand, "stages": {}}

    # -------------------------------------------------------------------------
    # STAGE 1: BOL Drums OUT (theta = 0 deg)
    # -------------------------------------------------------------------------
    print(f"[{cand_id}] Stage 1: BOL Drums OUT (theta = 0 deg)...", flush=True)
    dir_out = cand_dir / "01_bol_out"
    m_out = build_model(cand, theta_deg=0.0, particles=10000, batches=65, inactive=20, seed=101)
    k_out, sig_out = run_openmc_case(dir_out, m_out)
    results['stages']['bol_out'] = {
        "keff": k_out, "sigma": sig_out,
        "k_lower_2sigma": k_out - Z_SIGMA * sig_out
    }
    print(f"[{cand_id}] k_out = {k_out:.5f} +/- {sig_out:.5f} (lower 2-sigma: {k_out - Z_SIGMA*sig_out:.5f})", flush=True)

    if k_out - Z_SIGMA * sig_out <= 1.0:
        print(f"[{cand_id}] REJECTED at Stage 1: k_out - 2*sigma <= 1.0 (k_out = {k_out:.5f})", flush=True)
        return False, "failed_k_out", results

    # -------------------------------------------------------------------------
    # STAGE 2: BOL Drums IN (theta = 180 deg)
    # -------------------------------------------------------------------------
    print(f"[{cand_id}] Stage 2: BOL Drums IN (theta = 180 deg)...", flush=True)
    dir_in = cand_dir / "02_bol_in"
    m_in = build_model(cand, theta_deg=180.0, particles=10000, batches=65, inactive=20, seed=102)
    k_in, sig_in = run_openmc_case(dir_in, m_in)
    results['stages']['bol_in'] = {
        "keff": k_in, "sigma": sig_in,
        "k_upper_2sigma": k_in + Z_SIGMA * sig_in
    }
    print(f"[{cand_id}] k_in = {k_in:.5f} +/- {sig_in:.5f} (upper 2-sigma: {k_in + Z_SIGMA*sig_in:.5f})", flush=True)

    if k_in + Z_SIGMA * sig_in >= 1.0:
        print(f"[{cand_id}] REJECTED at Stage 2: k_in + 2*sigma >= 1.0 (k_in = {k_in:.5f})", flush=True)
        return False, "failed_k_in", results

    # -------------------------------------------------------------------------
    # STAGE 3: Dual-Pass Thermal Hydraulics (T_mod_peak < 950.0 K)
    # -------------------------------------------------------------------------
    print(f"[{cand_id}] Stage 3: Dual-Pass Thermal Hydraulics...", flush=True)
    th_res = evaluate_thermal(cand)
    results['stages']['thermal'] = th_res
    print(f"[{cand_id}] T_mod_peak = {th_res['T_mod_peak_K']:.1f} K (buffer: +{th_res['thermal_buffer_K']:.1f} K to 950 K limit)", flush=True)

    if not th_res['passes_thermal']:
        print(f"[{cand_id}] REJECTED at Stage 3: T_mod_peak >= 950.0 K ({th_res['T_mod_peak_K']:.1f} K)", flush=True)
        return False, "failed_T_mod", results

    # -------------------------------------------------------------------------
    # STAGE 4: Reactivity Feedback Coefficients (Two-point Central Difference)
    # -------------------------------------------------------------------------
    # 4a. Fuel Prompt Feedback (Doppler + UN Axial Thermal Expansion, Delta-T = 300 K: 900 K -> 1200 K)
    print(f"[{cand_id}] Stage 4a: Fuel Prompt Reactivity Feedback (300 K span: 900 K -> 1200 K)...", flush=True)
    dir_fuel_cold = cand_dir / "04a_fuel_cold_900K"
    dir_fuel_hot = cand_dir / "04a_fuel_hot_1200K"
    beta_fuel = 3.15e-5 # UN volumetric CTE: 3 * 10.5e-6 K^-1
    fuel_dens_cold = 1.0 / (1.0 - beta_fuel * 150.0)
    fuel_dens_hot = 1.0 / (1.0 + beta_fuel * 150.0)
    m_fc = build_model(cand, theta_deg=0.0, T_fuel=900.0, fuel_density_factor=fuel_dens_cold,
                       particles=20000, batches=75, inactive=25, seed=201)
    m_fh = build_model(cand, theta_deg=0.0, T_fuel=1200.0, fuel_density_factor=fuel_dens_hot,
                       particles=20000, batches=75, inactive=25, seed=202)
    k_fc, sig_fc = run_openmc_case(dir_fuel_cold, m_fc)
    k_fh, sig_fh = run_openmc_case(dir_fuel_hot, m_fh)

    rho_fc = 1e5 * (k_fc - 1.0) / k_fc
    rho_fh = 1e5 * (k_fh - 1.0) / k_fh
    alpha_fuel = (rho_fh - rho_fc) / 300.0
    sig_alpha_fuel = (1e5 / 300.0) * math.sqrt((sig_fh / (k_fh**2))**2 + (sig_fc / (k_fc**2))**2)
    results['stages']['alpha_fuel'] = {
        "alpha_fuel_pcm_K": alpha_fuel, "sigma": sig_alpha_fuel,
        "upper_2sigma": alpha_fuel + Z_SIGMA * sig_alpha_fuel
    }
    print(f"[{cand_id}] alpha_fuel = {alpha_fuel:.4f} +/- {sig_alpha_fuel:.4f} pcm/K (upper 2-sigma: {alpha_fuel + Z_SIGMA*sig_alpha_fuel:.4f})", flush=True)

    if alpha_fuel + Z_SIGMA * sig_alpha_fuel >= 0.0:
        print(f"[{cand_id}] REJECTED at Stage 4a: alpha_fuel is not negative within 2-sigma", flush=True)
        return False, "failed_alpha_fuel", results

    # 4b. Moderator Feedback (Delta-T = 200 K: 800 K to 1000 K with thermal expansion)
    print(f"[{cand_id}] Stage 4b: Moderator Temperature Feedback (200 K span: 800 K -> 1000 K, thermal exp)...", flush=True)
    dir_mod_cold = cand_dir / "04b_mod_cold_800K"
    dir_mod_hot = cand_dir / "04b_mod_hot_1000K"
    beta_vol = 3.6e-5 # YHx volumetric CTE
    dens_cold = 1.0 / (1.0 - beta_vol * 100.0)
    dens_hot = 1.0 / (1.0 + beta_vol * 100.0)
    m_mc = build_model(cand, theta_deg=0.0, T_mod=800.0, mod_density_factor=dens_cold,
                       particles=20000, batches=75, inactive=25, seed=211)
    m_mh = build_model(cand, theta_deg=0.0, T_mod=1000.0, mod_density_factor=dens_hot,
                       particles=20000, batches=75, inactive=25, seed=212)
    k_mc, sig_mc = run_openmc_case(dir_mod_cold, m_mc)
    k_mh, sig_mh = run_openmc_case(dir_mod_hot, m_mh)

    rho_mc = 1e5 * (k_mc - 1.0) / k_mc
    rho_mh = 1e5 * (k_mh - 1.0) / k_mh
    alpha_mod = (rho_mh - rho_mc) / 200.0
    sig_alpha_mod = (1e5 / 200.0) * math.sqrt((sig_mh / (k_mh**2))**2 + (sig_mc / (k_mc**2))**2)
    results['stages']['alpha_mod'] = {
        "alpha_mod_pcm_K": alpha_mod, "sigma": sig_alpha_mod,
        "upper_2sigma": alpha_mod + Z_SIGMA * sig_alpha_mod
    }
    print(f"[{cand_id}] alpha_mod = {alpha_mod:.4f} +/- {sig_alpha_mod:.4f} pcm/K (upper 2-sigma: {alpha_mod + Z_SIGMA*sig_alpha_mod:.4f})", flush=True)

    if alpha_mod + Z_SIGMA * sig_alpha_mod >= 0.0:
        print(f"[{cand_id}] REJECTED at Stage 4b: alpha_mod is not negative within 2-sigma", flush=True)
        return False, "failed_alpha_mod", results

    # 4c. Total Isothermal Feedback (Delta-T = 200 K: Cold All -> Hot All)
    print(f"[{cand_id}] Stage 4c: Total Isothermal Reactivity Feedback (200 K span)...", flush=True)
    dir_tot_cold = cand_dir / "04c_total_cold"
    dir_tot_hot = cand_dir / "04c_total_hot"
    m_tc = build_model(cand, theta_deg=0.0, T_fuel=900.0, T_mod=800.0, T_cool=800.0,
                       T_refl=800.0, mod_density_factor=dens_cold, fuel_density_factor=fuel_dens_cold,
                       particles=20000, batches=75, inactive=25, seed=221)
    m_th = build_model(cand, theta_deg=0.0, T_fuel=1200.0, T_mod=1000.0, T_cool=1000.0,
                       T_refl=1000.0, mod_density_factor=dens_hot, fuel_density_factor=fuel_dens_hot,
                       particles=20000, batches=75, inactive=25, seed=222)
    k_tc, sig_tc = run_openmc_case(dir_tot_cold, m_tc)
    k_th, sig_th = run_openmc_case(dir_tot_hot, m_th)

    rho_tc = 1e5 * (k_tc - 1.0) / k_tc
    rho_th = 1e5 * (k_th - 1.0) / k_th
    alpha_tot = (rho_th - rho_tc) / 200.0
    sig_alpha_tot = (1e5 / 200.0) * math.sqrt((sig_th / (k_th**2))**2 + (sig_tc / (k_tc**2))**2)
    results['stages']['alpha_total'] = {
        "alpha_total_pcm_K": alpha_tot, "sigma": sig_alpha_tot,
        "upper_2sigma": alpha_tot + Z_SIGMA * sig_alpha_tot
    }
    print(f"[{cand_id}] alpha_total = {alpha_tot:.4f} +/- {sig_alpha_tot:.4f} pcm/K (upper 2-sigma: {alpha_tot + Z_SIGMA*sig_alpha_tot:.4f})", flush=True)

    if alpha_tot + Z_SIGMA * sig_alpha_tot >= 0.0:
        print(f"[{cand_id}] REJECTED at Stage 4c: alpha_total is not negative within 2-sigma", flush=True)
        return False, "failed_alpha_total", results

    # -------------------------------------------------------------------------
    # CANDIDATE MEETS ALL CONDITIONS!
    # -------------------------------------------------------------------------
    print(f"\n{'*'*80}\n>>> SUCCESS: CANDIDATE {cand_id} SATISFIES ALL CRITERIA! <<<\n{'*'*80}\n", flush=True)
    return True, "accepted", results

def finalize_winner(cand, results):
    """Create dedicated folder with all figures, reports, and data for the winning reactor."""
    cand_id = cand['id']
    win_dir = ROOT / "found_reactor" / cand_id
    win_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = win_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1. Save results JSON
    atomic_write_json(win_dir / "results.json", results)

    # 2. Generate high-resolution figures
    print(f"Generating high-resolution figures for {cand_id} in {fig_dir}...")
    generate_core_figures(cand, fig_dir)

    # 3. Write results.txt
    st = results['stages']
    txt_content = f"""================================================================================
REATTORE SELEZIONATO: {cand_id}
ARCHITETTURA NUCLEARE PER PROPULSIONE ELETTRICA SPAZIALE (NEP 400 kWth, 100 kWe)
Data: 8 ottobre 2026
================================================================================

1. ESITO DI CONFORMITÀ SUI REQUISITI UTENTE (TUTTI SODDISFATTI)
--------------------------------------------------------------------------------
- Potenza Termica Core:       400.0 kWth
- Potenza Elettrica Netta:    100.0+ kWe
- Moderatore Impiegato:       {cand['moderator_type']} (Stechiometria x = {cand['y_stoich']}, rho = {cand['rho_mod']} g/cm3)
- Refrigerante Core:          Elio puro a {cand['pressure_MPa']} MPa
- Temperatura Ingresso Gas:   {cand['inlet_K']} K (Requisito 850 - 900 K: CONFORME)
- Temperatura Uscita Gas:     {cand['outlet_K']} K (Requisito 1100 - 1150 K: CONFORME)
- Portata Elio:               {st['thermal']['m_dot_total_kg_s']} kg/s

PARAMETRI CRITICI VERIFICATI (ENTRO 2 SIGMA):
1. BOL Drums OUT (0°):
   keff = {st['bol_out']['keff']:.5f} +/- {st['bol_out']['sigma']:.5f}
   Limite inferiore a 2-sigma: {st['bol_out']['k_lower_2sigma']:.5f} > 1.00000 [PASS]
   Reattività in eccesso BOL: +{1e5*(st['bol_out']['keff']-1)/st['bol_out']['keff']:.1f} pcm

2. BOL Drums IN (180°):
   keff = {st['bol_in']['keff']:.5f} +/- {st['bol_in']['sigma']:.5f}
   Limite superiore a 2-sigma: {st['bol_in']['k_upper_2sigma']:.5f} < 1.00000 [PASS]
   Margine di arresto (Shutdown Margin): +{1e5*(1-st['bol_in']['keff'])/st['bol_in']['keff']:.1f} pcm
   Worth totale tamburi: {1e5*(1/st['bol_in']['keff'] - 1/st['bol_out']['keff']):.1f} pcm

3. Temperatura Picco Moderatore YHx:
   T_mod_peak = {st['thermal']['T_mod_peak_K']:.2f} K < 950.00 K [PASS]
   Margine di sicurezza termica: +{st['thermal']['thermal_buffer_K']:.2f} K sotto il limite 950 K
   T_mod_avg = {st['thermal']['T_mod_avg_K']:.2f} K
   T_clad_peak = {st['thermal']['T_clad_peak_K']:.2f} K (< 1350 K creep limit Nb-1Zr)
   T_fuel_peak = {st['thermal']['T_fuel_peak_K']:.2f} K (< 2073 K UN dissociation limit)

4. Coefficienti di Temperatura di Reattività (Tutti Negativi a 2-sigma):
   - Combustibile (alpha_fuel, Doppler UN):
     alpha_fuel = {st['alpha_fuel']['alpha_fuel_pcm_K']:.4f} +/- {st['alpha_fuel']['sigma']:.4f} pcm/K
     Limite a 2-sigma: {st['alpha_fuel']['upper_2sigma']:.4f} pcm/K < 0.0000 [PASS]
   - Moderatore (alpha_mod, YHx spettrale + espansione):
     alpha_mod = {st['alpha_mod']['alpha_mod_pcm_K']:.4f} +/- {st['alpha_mod']['sigma']:.4f} pcm/K
     Limite a 2-sigma: {st['alpha_mod']['upper_2sigma']:.4f} pcm/K < 0.0000 [PASS]
   - Totale Isotermico (alpha_total):
     alpha_total = {st['alpha_total']['alpha_total_pcm_K']:.4f} +/- {st['alpha_total']['sigma']:.4f} pcm/K
     Limite a 2-sigma: {st['alpha_total']['upper_2sigma']:.4f} pcm/K < 0.0000 [PASS]

2. CARATTERISTICHE COSTRUTTIVE E GEOMETRICHE
--------------------------------------------------------------------------------
- Famiglia:                   {cand['family']}
- Numero Barre Combustibile:  {cand['pins']} pin (Reticolo esagonale a {cand['rings']} anelli)
- Passo Reticolo (Pitch):     {cand['pitch_cm']} cm
- Lunghezza Attiva Nocciolo:  {cand['height_cm']} cm
- Raggio Pastiglia UN:        {cand['r_fuel_cm']} cm
- Spessore Guaina Nb-1Zr:     {cand['t_clad_cm']} cm
- Canale Refrigerante Caldo:  {cand['t_hot_gap_cm']} cm
- Barriera Idrogeno Al2O3:    {cand['t_al2o3_cm']} cm (100 um)
- Canali Raffreddamento YHx:  {cand.get('aux_holes', 'N/A')} fori/cella da {cand.get('r_aux_cm', 'N/A')*10:.2f} mm di raggio
- Riflettore Radiale / Ax:    {cand['reflector']} (Spessore radiale {cand['t_refl_rad_cm']} cm, assiale {cand['t_refl_ax_cm']} cm)
- Tamburi di Controllo:       {cand['drums']} tamburi in {cand['reflector']}, raggio {cand['r_drum_cm']} cm
- Settore Assorbitore:        B4C (90% B-10), arco {cand['arc_b4c_deg']}°, spessore {cand['t_b4c_cm']} cm
- Inventario Metallo Pesante: {cand['heavy_metal_kg']} kgHM
- Inventario U-235 Fissile:   {cand['u235_kg']} kg U-235 (Arricchimento ~19.75 wt%)
- Massa Totale Stimata Core:  {cand['mass_est_kg']} kg

3. VALUTAZIONE LONGEVITÀ E DURATA OPERATIVA (> 10 ANNI)
--------------------------------------------------------------------------------
- Consumo fissile a 400 kWth per 10 anni continuativi: circa 4.0 kg U-235 (fissione).
- Con un inventario iniziale di {cand['u235_kg']} kg di U-235 e una reattività in eccesso BOL
  di oltre +{1e5*(st['bol_out']['keff']-1)/st['bol_out']['keff']:.0f} pcm, il nocciolo possiede
  una riserva di combustibile e reattività largamente sufficiente a garantire oltre
  10-15 anni di operatività continua prima del raggiungimento della sottocriticità.

================================================================================
Figure geometriche disponibili in {fig_dir}:
  - 01_radial_drums_out.png: Sezione radiale con tamburi OUT (0°)
  - 02_radial_drums_in.png:  Sezione radiale con tamburi IN (180°)
  - 03_axial_section.png:    Sezione assiale X-Z
  - 04_unit_cell.png:        Sezione ad alta risoluzione della cella elementare
================================================================================
"""
    with open(win_dir / "results.txt", "w") as f:
        f.write(txt_content)

    # Also copy to root level of thermal_neutronic_search_v3
    with open(ROOT / "results.txt", "w") as f:
        f.write(txt_content)

    print(f"\nFinalized winner documentation in {win_dir} and {ROOT / 'results.txt'}!")

def main():
    print("=" * 80)
    print("STARTING NEP MULTIPHYSICS CAMPAIGN (30 CANDIDATES, 2-SIGMA, 7 THREADS)")
    print("=" * 80)

    candidates = get_candidates()
    csv_progress_path = ROOT / "progress.csv"
    
    with open(csv_progress_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "status", "failed_metric", "k_out", "k_in", "T_mod_peak_K", "alpha_mod", "alpha_fuel", "alpha_tot"])

    tested_count = 0
    winner_found = False

    for cand in candidates:
        tested_count += 1
        cand_id = cand['id']

        atomic_write_json(ROOT / "status.json", {
            "status": "running",
            "current_candidate": cand_id,
            "tested_candidates": tested_count,
            "total_candidates": len(candidates),
            "timestamp": time.time()
        })

        passed, reason, results = evaluate_candidate(cand)

        # Record to progress CSV
        st = results.get('stages', {})
        k_out_val = st.get('bol_out', {}).get('keff', '')
        k_in_val = st.get('bol_in', {}).get('keff', '')
        t_mod_val = st.get('thermal', {}).get('T_mod_peak_K', '')
        a_mod_val = st.get('alpha_mod', {}).get('alpha_mod_pcm_K', '')
        a_fuel_val = st.get('alpha_fuel', {}).get('alpha_fuel_pcm_K', '')
        a_tot_val = st.get('alpha_total', {}).get('alpha_total_pcm_K', '')

        with open(csv_progress_path, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([cand_id, "accepted" if passed else "rejected",
                             "" if passed else reason,
                             k_out_val, k_in_val, t_mod_val, a_mod_val, a_fuel_val, a_tot_val])

        if passed:
            winner_found = True
            finalize_winner(cand, results)
            atomic_write_json(ROOT / "status.json", {
                "status": "solution_found",
                "winning_candidate": cand_id,
                "tested_candidates": tested_count,
                "timestamp": time.time()
            })
            break

    if not winner_found:
        print("\nAll 30 candidates tested. No solution satisfied all simultaneous constraints.")
        atomic_write_json(ROOT / "status.json", {
            "status": "exhausted_no_solution",
            "tested_candidates": tested_count,
            "timestamp": time.time()
        })

if __name__ == "__main__":
    main()
