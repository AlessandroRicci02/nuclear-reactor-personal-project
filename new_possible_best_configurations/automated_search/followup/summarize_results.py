"""Keep a readable, honestly qualified record of the follow-up calculations."""
from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPORT = HERE / "results.txt"


def read_json(name):
    p = HERE / name
    return json.loads(p.read_text()) if p.exists() else None


def report():
    mass = read_json("mass_result.json")
    heating = read_json("heating_result.json")
    thermal = read_json("thermal_results.json")
    life = read_json("lifespan_result.json")
    status = read_json("depletion/status.json") or {}
    csv_path = HERE / "depletion_2p5yr.csv"
    rows = list(csv.DictReader(csv_path.open())) if csv_path.exists() else []
    out = ["FOLLOW-UP DEL CANDIDATO matrix_short_strong_be",
           "Materiale moderatore: YH2.0 (idruro di ittrio, Y:H = 1:2); riflettore: Be.",
           "Potenza termica assunta: 400 kW. Elio: 800 K ingresso, 1150 K uscita miscelata.",
           ""]
    if mass:
        out += ["MASSA",
                f"Solidi esplicitamente modellati: {mass['modeled_dry_mass_kg']:.2f} +/- {mass['sampling_sigma_kg']:.2f} kg (1 sigma campionamento geometrico).",
                f"Con inventario di elio: {mass['core_plus_helium_kg']:.2f} kg.",
                f"Stima preliminare con margini hardware (+65 kg) e contingenza (+15 kg): {mass['screening_mass_with_allowances_kg']:.2f} kg.",
                "Non e' la massa completa del sistema: schermatura, attuatori, conversione elettrica, radiatori, tubazioni esterne e struttura non sono progettati/contabilizzati.",
                ""]
    if heating:
        out += ["DEPOSIZIONE DI CALORE OPENMC",
                f"Angolo BOL di mappa: {heating['theta_deg']:.3f} gradi; k_eff = {heating['k_eff']:.6f} +/- {heating['sigma_k']:.6f}.",
                f"Deposizione normalizzata: fuel {heating['fuel_power_W']:.0f} W, YH2 {heating['moderator_power_W']:.0f} W, altri materiali {heating['other_power_W']:.0f} W.",
                ""]
    if thermal and thermal.get("cases"):
        cases = thermal["cases"]
        base = next((c for c in cases if c["label"] == "a0.35_o0.00_f1.00_k10"), cases[0])
        best = min(cases, key=lambda c: c["max_moderator_K"])
        worst = max(cases, key=lambda c: c["max_moderator_K"])
        out += ["TEMPERATURA DEL MODERATORE (OpenMC + Cardinal, steady state)",
                f"Scenario di riferimento: 35% portata He ai 6 canali ausiliari, h e conducibilita' nominali.",
                f"T massima YH2 = {base['max_moderator_K']:.1f} K; T He miscelato = {base['mixed_outlet_K']:.1f} K.",
                f"Nei {len(cases)} scenari espliciti: {best['max_moderator_K']:.1f}--{worst['max_moderator_K']:.1f} K.",
                "Il limite richiesto <950 K NON e' soddisfatto negli scenari analizzati.",
                "Il valore e' condizionato alla ripartizione della portata, a coefficienti di scambio e conduzione assunti, e al modello 2D adiabatico della cella; mancano collettori e CFD di sistema.",
                "Il modello neutronico BOL assumeva 875 K per YH2. Poiche' la stima termica e' assai piu' alta, la criticita' e la durata neutronica devono essere ricalcolate con temperature coerenti prima di usarle per un giudizio di progetto.",
                ""]
    if rows:
        out += ["DEPLETION NEUTRONICA OPENMC, PASSO 2.5 ANNI",
                "0 gradi = assorbitore verso spazio; 180 gradi = verso nocciolo.",
                "Anni   k(0 gradi) +/-sigma  k(180 gradi) +/-sigma  angolo critico  k(angle) +/-sigma  stato"]
        for r in rows:
            def v(key):
                return f"{float(r[key]):.6f}" if r.get(key) else "--"
            out.append(f"{float(r['years']):4.1f}   {v('k_out')} +/-{v('sigma_k_out')}   {v('k_in')} +/-{v('sigma_k_in')}   {v('critical_angle_deg')}   {v('k_critical')} +/-{v('sigma_k_critical')}   {r['status']}")
        out += ["Ogni intervallo usa 400 kW e l'angolo critico trovato all'inizio; il modello termico neutronico e' tenuto al profilo BOL.",
                ""]
    if life:
        if life["status"] == "bracketed":
            lo, hi = life["interval_years"]
            out += [f"Durata neutronica condizionata: perdita di criticita' tra {lo:.1f} e {hi:.1f} anni.",
                    f"Interpolazione lineare di k(0 gradi): {life['linear_interpolation_years']:.2f} anni (non risoluzione temporale aggiuntiva)."]
        else:
            out += [f"Durata NON delimitata entro l'ultima simulazione: {life.get('last_sample_years','?')} anni."]
    else:
        out += [f"Durata: calcolo in corso; ultimo stato = {status.get('status','non disponibile')}, fase = {status.get('stage','--')}, passo = {status.get('step','--')}."]
        if status.get("status") == "failed":
            out += [f"Errore: {status.get('error')}"]
    out += ["", "Fonti numeriche: mass_result.json, heating_result.json, thermal_results.json, depletion_2p5yr.csv, lifespan_result.json."]
    tmp = REPORT.with_suffix(".txt.tmp")
    tmp.write_text("\n".join(out) + "\n")
    tmp.replace(REPORT)
    return bool(life) or status.get("status") == "failed"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    if args.watch:
        while True:
            if report():
                break
            time.sleep(60)
    else:
        report()
