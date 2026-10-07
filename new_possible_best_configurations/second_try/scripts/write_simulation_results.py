"""Write the requested observations only when every campaign has completed."""
from pathlib import Path
import json
from s2_model import ROOT

CAM=ROOT/'campaign'
def read(path):return json.loads(Path(path).read_text())

def main():
    mass=read(CAM/'mass/result.json')
    coeff=read(CAM/'feedback_coefficients.json')
    status=read(CAM/'status.json')
    heat=read(CAM/'thermal/deposition_summary.json')
    th=read(CAM/'thermal/thermal_results.json')
    if status['stage']!='complete' or len(coeff)!=3 or th['status']!='complete':
        raise RuntimeError('Campaign incomplete; preserve in-progress results.txt')
    pilot=max((CAM/'pilot').glob('statepoint.*.h5'),key=lambda p:int(p.name.split('.')[1]))
    cases={r['label']:r for r in th['cases']}
    b=cases['baseline']
    lines=[
        'S2-R1 — RISULTATI DELLE SIMULAZIONI RICHIESTE',
        'Configurazione: configuration.json; raffreddamento nominale 800 -> 1150 K, 400 kWth imposti.',
        'ATTENZIONE: la geometria è già nettamente sottocritica con tutti i tamburi a 0°.',
        f"OpenMC BOL con heating neutroni+fotoni: k_eff={heat['keff']:.6f} +/- {heat['sigma_k']:.6f} (1 sigma MC).",
        'Non può sostenere autonomamente 400 kW nella configurazione attuale.',
        '',
        'MASSA DEL REATTORE',
        f"Solidi modellati, OpenMC: {mass['solid_mass_kg']:.3f} +/- {mass['sampling_sigma_kg']:.3f} kg (1 sigma stocastica).",
        f"Solidi dalla geometria analitica: {mass['analytic_mass_kg']:.3f} kg.",
        f"Elio nei volumi modellati: {mass['helium_mass_kg']:.5f} kg.",
        f"Solidi + 65 kg hardware NON dimensionato + 15 kg riserva: {mass['conditional_total_kg']:.3f} kg circa.",
        'La sola incertezza stocastica NON copre massa di hardware, tolleranze, schermatura e interfacce.',
        'Esclusi: scudo, conversione, radiatori, tubi esterni e struttura del veicolo.',
        '',
        'COEFFICIENTI TERMICI (derivate matematiche del k_eff della geometria sottocritica)',
        'Definizione: differenze finite centrali della reattività rho=(k-1)/k; tamburi tutti OUT; geometria e densità dei solidi fisse.',
        'Nel totale cambiano tutte le temperature dei materiali; densità He alla stessa pressione.',
    ]
    for c in coeff:
        sign='negativo al 95%' if c['high_95_pcm_K']<0 else 'positivo al 95%' if c['low_95_pcm_K']>0 else 'segno non risolto al 95%'
        lines.append(f"alpha_{c['kind']}: {c['alpha_pcm_K']:+.4f} +/- {c['sampling_sigma_pcm_K']:.4f} pcm/K (1 sigma MC); IC95 [{c['low_95_pcm_K']:+.4f}, {c['high_95_pcm_K']:+.4f}]; {sign}.")
    lines.extend(['Le incertezze comprendono il campionamento Monte Carlo degli endpoint, non i dati nucleari o le proprietà termo-meccaniche.',
                  'Le temperature di riferimento sono in campaign/feedback/*/manifest.json. La perturbazione +100 K di YH2 supera 950 K: è un punto matematico di derivata, non una condizione operativa ammessa.',
                  '',
                  'TEMPERATURA YH2 — SCENARIO CONDIZIONALE A 400 kW IMPOSTI',
                  f"Deposizione OpenMC neutroni+fotoni normalizzata a 400 kW: {heat['statepoint']}.",
                  f"Modello seriale del flusso: uscita miscelata imposta {b['hot_mixed_outlet_K']:.3f} K; massimo canale caldo {b['max_hot_outlet_K']:.3f} K.",
                  f"Caso base: massimo YH2 analitico per tutti i pin {b['max_YH2_analytic_K']:.3f} K, pin {b['max_YH2_pin_id']} slice {b['max_YH2_slice']}.",
                  f"Cardinal sul pin selezionato: {b['cardinal_peak_pin']['T_YH2_max_K']:.3f} K; raffinato: {b['cardinal_peak_pin_refined']['T_YH2_max_K']:.3f} K.",
                  f"Uscita media passaggio freddo {b['cold_mixed_outlet_K']:.3f} K; massimo locale bulk freddo {b['max_cold_bulk_K']:.3f} K.",
                  f"Calore totale primo passaggio {b['cold_total_W']/1000:.3f} kW, di cui perdite attraverso isolamento {b['cold_leak_W']/1000:.3f} kW.",
                  'Lettura della temperatura: il valore analitico è uno screening globale che assegna conservativamente tutto il riscaldamento diretto del ramo freddo al manicotto YH2; Cardinal risolve il pin selezionato.',
                  'La temperatura non è una previsione validata: emissività, vuoto, ponti termici, Nusselt, conducibilità, mappa Be e flusso tra pin sono ipotesi.',
                  'Sensibilità, una modifica alla volta:',
                 ])
    for name in ('baseline','oxidized','hotter_wall','high_bridge','low_film','low_k_YH2','all_Be_to_cold'):
        c=cases[name]
        lines.append(f"  {name}: YH2 max screening {c['max_YH2_analytic_K']:.3f} K; uscita fredda {c['cold_mixed_outlet_K']:.3f} K; massimo He caldo {c['max_hot_outlet_K']:.3f} K.")
    lines.extend(['',
        'LIMITI CHE IMPEDISCONO UNA CONCLUSIONE DI CONFORMITÀ',
        '1. k_eff ~0.82 a tamburi OUT: la potenza imposta è uno scenario ipotetico. La candidata non supera il filtro di criticità.',
        '2. La camera evacuata a 2.5 MPa differenziale, la sua tenuta per >15 anni e le proprietà ottiche non sono validate.',
        '3. Il modello termico ridotto non risolve 3D del nocciolo, collettori, temperatura periferica Be o transitori. Neutron/gamma heating è un input reale, il campo termico resta condizionale.',
        '4. La massa <600 kg include una sola assegnazione di hardware, non un assieme CAD/strutturale dimensionato.',
        '',
        'DATI RIPRODUCIBILI',
        'campaign/mass/: volume OpenMC, tabella materiali e risultato.',
        'campaign/feedback/: sei simulazioni Monte Carlo e statepoint.',
        'campaign/heating/: trasporto accoppiato neutroni+fotoni e tally di deposizione.',
        'campaign/thermal/: campi del modello seriale, Cardinal, bilanci e sensitività.',
        'scripts/s2_model.py, run_neutronics.py, run_thermal.py: input e controller.',
    ])
    out=ROOT/'results.txt';tmp=ROOT/'results.txt.tmp';tmp.write_text('\n'.join(lines)+'\n');tmp.replace(out)
    print('Wrote',out)

if __name__=='__main__':main()
