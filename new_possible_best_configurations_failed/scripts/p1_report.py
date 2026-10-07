"""Regenerate comparison.txt from real result files, including unfinished stages."""
import csv
from datetime import datetime
import json
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT=Path(__file__).resolve().parents[1]
CAMPAIGN=ROOT/'campaign'
BASE=ROOT.parent/'best_configurations/GEOM_OPT15_TB35'


def read(path):
    p=Path(path)
    try:return json.loads(p.read_text()) if p.exists() else None
    except json.JSONDecodeError:return None # Writer is still completing this file; retry next refresh.


def rows(path):
    with Path(path).open() as f:return list(csv.DictReader(f))


def estimate(r):
    return f"{r['keff']:.7f} +/- {r['sigma']:.7f} (1 sigma, {r['sigma']*1e5:.1f} pcm in k)" if r else 'PENDING: no completed precision-qualified result'


def report():
    state=read(CAMPAIGN/'status.json') or {'status':'BOL_RUNNING','stage':'initial OpenMC production'}
    t=datetime.now(ZoneInfo('Europe/Rome')).isoformat(timespec='seconds')
    lines=['P1 ELLIPTIC BALANCED versus GEOM_OPT15_TB35',f'Updated: {t}',f"Campaign status: {state['status']}; stage: {state.get('stage','unknown')}",
           'Only completed solver outputs are reported as results. PENDING is not a result.',
           'Uncertainties on keff and feedback are Monte Carlo standard errors only; no nuclear-data covariance or thermal-property uncertainty has been propagated.','',
           'MODEL AND MASS']
    if state.get('active_transport'):lines.insert(3,'Active OpenMC calculation: '+state['active_transport'])
    mass=read(CAMPAIGN/'data/mass_summary.json')
    oldcfg=read(ROOT/'configurations/TB35_REFERENCE/geometry.json')
    newcfg=read(ROOT/'configurations/P1_ELLIPTIC_BALANCED/geometry.json')
    lines.extend([f"P1 modeled solids: {mass['analytic_solids_kg']:.3f} kg; stochastic check {mass['stochastic_solids_kg']:.3f} +/- {mass['stochastic_sigma_kg']:.3f} kg.",
                  f"TB35, same Boolean-union inventory convention: {oldcfg['metrics']['modeled_solids_kg']:.3f} kg.",
                  f"P1 minus TB35: {mass['analytic_solids_kg']-oldcfg['metrics']['modeled_solids_kg']:.3f} kg.",
                  f"Both retain {newcfg['metrics']['hm_kg']:.3f} kg initial uranium and {newcfg['metrics']['u235_kg']:.3f} kg U235.",
                  'These are modeled solid inventories, excluding pressure vessel, drives, manifolds, supports, extra hydrogen encapsulation, shielding and power conversion.','',
                  'BEGINNING-OF-LIFE NEUTRONICS AT PRESCRIBED NOMINAL TEMPERATURES'])
    historic=rows(BASE/'data/tb35_bol_neutronics_summary.csv')[0]
    for name,label in [('p1_bol_0deg','P1, drums outward (0 deg)'),('p1_bol_180deg','P1, drums inward (180 deg)'),('tb35_corrected_bol_0deg','TB35, corrected common data, outward'),('tb35_corrected_bol_180deg','TB35, corrected common data, inward')]:
        lines.append(label+': '+estimate(read(CAMPAIGN/'neutronics'/name/'result.json')))
    lines.extend([f"Historical TB35 outward: {float(historic['k_eff_0deg']):.7f} +/- {float(historic['sigma_0deg']):.7f}.",
                  f"Historical TB35 inward: {float(historic['k_eff_180deg']):.7f} +/- {float(historic['sigma_180deg']):.7f}.",
                  'Historical TB35 used the previous room-temperature/WMP library and ZrH scattering substitute. New common-data runs use ENDF/B-VIII.0, temperature interpolation and YH2 scattering. Historical-to-new differences cannot be attributed solely to geometry.'])
    for name in ('p1','tb35'):
        c=read(CAMPAIGN/f'neutronics/{name}_control/control.json')
        if c:
            lines.append(f"{name.upper()} criticality control: {c['status']}; angle: {c['theta_deg']} deg.")
            if c['critical']:lines.append('Independent verification: '+estimate(c['verification']))
        else:lines.append(name.upper()+' critical drum angle: PENDING.')
    lines.extend(['','THERMAL AND DENSITY REACTIVITY COEFFICIENTS'])
    feedback=read(CAMPAIGN/'data/feedback_coefficients.json')
    if feedback:
        for r in feedback:lines.append(f"P1 {r['name']}: {r['value']:+.5g} +/- {r['sigma']:.4g} {r['units']} (1 sigma); sign resolved at 95%: {r['negative_at_95pct'] or r['positive_at_95pct']}.")
    else:lines.append('P1 feedback coefficients: PENDING.')
    lines.append('Historical TB35 coefficients (different nuclear-data treatment):')
    for r in rows(BASE/'data/tb35_feedback_coefficients.csv'):
        v,s=float(r['coefficient']),float(r['sigma_coefficient'])
        sign='negative' if v+1.96*s<0 else 'positive' if v-1.96*s>0 else 'unresolved'
        lines.append(f"  {r['component']}: {v:+.5g} +/- {s:.4g} {r['coeff_unit']}; 95% sign: {sign}.")
    lines.extend(['','DEPLETION AND LIFETIME'])
    dep=read(CAMPAIGN/'data/depletion_summary.json')
    if dep:
        lines.append('P1 depletion status: '+dep['status'])
        if dep.get('rows'):
            lines.append(f"Last completed interval: {dep['rows'][-1]['time_years']:.5f} equivalent full-power years.")
            lines.append('2.5-year reporting points: years, drum deg, keff +/- sigma')
            for r in dep['rows']:
                if r['report_2_5yr']:lines.append(f"  {r['time_years']:.1f}, {r['theta_deg']:.4f}, {r['keff']:.7f} +/- {r['sigma']:.7f}")
        if dep.get('lifetime_bracket_years'):lines.append(f"P1 control-limit lifetime bracket: {dep['lifetime_bracket_years']} equivalent full-power years. This is a time-step bracket, not a confidence interval.")
        if dep['status']=='NO_CRITICAL_BOL_STATE':lines.append('No operable critical beginning-of-life state exists with these drums; fuel-exhaustion lifetime is undefined. No fictitious critical depletion trajectory is reported.')
        if dep.get('lifetime_lower_bound_years'):lines.append(f"No end of life observed through {dep['lifetime_lower_bound_years']} equivalent full-power years; this is a lower bound.")
    else:lines.append('P1 lifespan and 2.5-year drum schedule: PENDING.')
    lines.append('P1 method: full thermal-branching ENDF/B-VIII.0 chain, 30 fuel compositions, CE/CM, 400 kW recoverable-fission power; startup steps 0.25/0.75/2/7/20/61.3125 days, then <=0.25 yr. Each transport stage checks criticality within a 75 pcm 95% statistical band and respects the 0–180 deg limits. This is quasi-static numerical control, not reactor dynamics.')
    histdep=rows(BASE/'data/tb35_criticality_controlled_depletion.csv')
    for a,b in zip(histdep,histdep[1:]):
        ka,kb=float(a['k_eff_drums_out_0deg']),float(b['k_eff_drums_out_0deg'])
        if ka>=1>kb:
            ta,tb=float(a['time_years']),float(b['time_years']);cross=ta+(1-ka)/(kb-ka)*(tb-ta)
            lines.append(f'Historical TB35 outward-keff crossing: approximately {cross:.2f} years by linear interpolation, bracketed by {ta:g} and {tb:g} years. Not an exact or uncertainty-qualified lifespan.')
            break
    deviations=[r for r in histdep if r['status']=='CRITICAL (Operational)' and abs(float(r['k_eff_operational_crit'])-1)>.001]
    lines.append(f'Historical TB35 depletion uses chain_nep.xml (12 nuclides) and single 2.5-year predictor steps. {len(deviations)} rows labeled critical exceed its own 100 pcm mean-keff tolerance. It is not a like-for-like full-chain benchmark.')
    lines.extend(['','MOOSE/CARDINAL THERMAL-HYDRAULICS'])
    for name in ('p1','tb35'):
        th=read(CAMPAIGN/f'thermal/{name}/summary.json')
        if not th:lines.append(name.upper()+': PENDING actual heating-driven Cardinal results.');continue
        assert not th.get('validation_only')
        operating=read(CAMPAIGN/f'data/{name}_heating_operating_status.json')
        if operating and not operating['critical_operating_state']:lines.append('HYPOTHETICAL PRESCRIBED POWER: this geometry cannot reach criticality at the chosen nominal state.')
        lines.extend([name.upper()+': '+th['model'],f"  Active coolant mixed outlet: {th['mixed_active_coolant_outlet_K']:.3f} K; energy residual: {th['energy_residual_W']:.4g} W.",
                      f"  Peak pin-cell moderator/liner solid temperature: {th['peak_pin_cell_solid_K']:.3f} K (not the full-core moderator peak).",
                      f"  Largest channel friction losses: main {th['max_main_friction_Pa']/1000:.4f} kPa; auxiliary {th['max_aux_friction_Pa']/1000:.4f} kPa. Not total reactor pressure drop."])
        for r in th['pin_fem']:lines.append(f"  RZ FEM pin {r['pin_id']}: fuel max {r['T_fuel_max']:.3f} K; clad max {r['T_clad_max']:.3f} K; heat balance relative residual {r['relative_energy_residual']:.3g}.")
        lines.append(f"  Refined peak-pin mesh temperature change: {th['pin_mesh_change_peak_K']:.4f} K.")
        lines.append(f"  Refined moderator-cell mesh peak change: {th.get('cell_mesh_change_peak_K',float('nan')):.4f} K.")
        for sensitivity in th.get('property_sensitivities',[]):lines.append(f"  Sensitivity {sensitivity['case']}: peak pin-cell solid {sensitivity['peak_pin_cell_solid_K']:.3f} K; 1D peak fuel {sensitivity['peak_fuel_1D_K']:.3f} K. Scenario change, not a statistical uncertainty.")
        lines.extend('  Model limit: '+s for s in th['limitations'])
    hth=rows(BASE/'data/tb35_moose_cardinal_pin_benchmarks.csv')[0]
    lines.extend([f"Historical TB35 peak-pin Cardinal calculation: fuel {float(hth['cardinal_T_fuel_max_K']):.3f} K; clad {float(hth['cardinal_T_clad_max_K']):.3f} K.",
                  'Its pin powers and flow multipliers were prescribed. Moderator peak 914.83 K and core pressure loss 13.5 kPa are assigned constants in its script, not independently converged coupled predictions.',
                  'The new and historical thermal results differ in heat deposition and flow accounting. The new common-model TB35 calculation is the appropriate geometry comparison.',
                  'The reported depletion lifetime is a neutronic control-limit lifetime at prescribed temperatures. Thermal/material limits may restrict service life earlier.',
                  'Net electric power and launch-system performance are not newly established by these neutronic and reduced thermal calculations.','',
                  'FILES AND REPRODUCTION',
                  'campaign/status.json: stage status and error details',
                  'campaign/logs/campaign.log: queue log; per-stage logs alongside it',
                  'campaign/neutronics/*/result.json: precision-qualified keff and statepoint provenance',
                  'campaign/data/feedback_coefficients.json: propagated independent-endpoint uncertainties',
                  'campaign/depletion/interval_*/: restartable materials, full depletion HDF5 and control history',
                  'campaign/thermal/{p1,tb35}/: deposited heating, meshes, MOOSE inputs/logs, fluid balances and FEM outputs',
                  'Resume with: /home/aleri/miniforge3/envs/openmc-env/bin/python -u scripts/p1_campaign.py',
                  'Existing OpenMC, standalone MOOSE, and Cardinal installations are reused. Native Cardinal/OpenMC interface validation is recorded separately under campaign/validation/cardinal_openmc_statepoint.',
                  'No historical TB35 results are modified.'])
    if state.get('errors'):lines.extend(['','CAMPAIGN ERRORS']+[f"{k}: {v}" for k,v in state['errors'].items()])
    path=ROOT/'comparison.txt';tmp=path.with_suffix('.txt.tmp');tmp.write_text('\n'.join(lines)+'\n');tmp.replace(path)


if __name__=='__main__':report()
