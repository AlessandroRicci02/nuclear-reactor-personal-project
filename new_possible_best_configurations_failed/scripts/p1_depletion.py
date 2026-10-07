"""Checkpointed CE/CM depletion with a verified drum search at every transport stage.

30 independent fuel compositions; full chain; 400 kW recoverable-fission power.
Reports at 2.5 years; shorter integration steps resolve startup and burnup.
The controller is numerical quasi-static criticality, not a dynamics simulation.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import time
import numpy as np
import openmc
import openmc.deplete
import openmc.lib
from p1_model import build, CAMPAIGN, DEFAULT_XS

CHAIN=CAMPAIGN/'data/chain_endfb80_thermal.xml'


def atomic_json(path,value):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


class ControlUnavailable(RuntimeError):
    pass


class ControlledOperator(openmc.deplete.CoupledOperator):
    def __init__(self, model, theta, directory):
        super().__init__(model, str(CHAIN), normalization_mode='fission-q', fission_yield_mode='constant')
        self.theta=theta
        self.directory=directory
        self.drums=[(c.id,int(c.name.split('_')[1])*36.) for c in model.geometry.get_all_cells().values() if c.name.startswith('drum_')]
        assert len(self.drums)==10
        self.stage=0
        self.records=[]

    def evaluate(self, vec, power, theta):
        for cid,az in self.drums:openmc.lib.cells[cid].rotation=(0.,0.,az+theta)
        openmc.lib.settings.seed=173+7919*(len(self.records)+1)+100003*self.directory_index
        answer=super().__call__(vec,power)
        k,sd=answer.k.n,answer.k.s
        if sd>2e-4:raise RuntimeError(f'Depletion precision target missed: {sd}')
        r=dict(stage=self.stage,theta_deg=theta,keff=k,sigma=sd,reactivity_pcm=1e5*(1-1/k),reactivity_sigma_pcm=1e5*sd/k**2)
        self.records.append(r)
        (self.directory/'control_history.json').write_text(json.dumps(self.records,indent=2)+'\n')
        return answer,r

    def __call__(self,vec,power):
        if power==0:return super().__call__(vec,power)
        self.stage+=1
        first,r=self.evaluate(vec,power,self.theta)
        def accepted(r):return abs(r['reactivity_pcm'])+1.96*r['reactivity_sigma_pcm']<=75.
        if accepted(r):return first
        # Establish physical control limits at the CURRENT isotopic composition.
        _,out=self.evaluate(vec,power,0.)
        _,inside=self.evaluate(vec,power,180.)
        if out['keff']+1.96*out['sigma']<1 or inside['keff']-1.96*inside['sigma']>1:
            (self.directory/'control_unavailable.json').write_text(json.dumps({'stage':self.stage,'out':out,'in':inside},indent=2)+'\n')
            raise ControlUnavailable('Current composition cannot be made critical within 0–180 degrees')
        lo,hi=(0.,out),(180.,inside)
        if r['keff']>1:lo=(self.theta,r)
        else:hi=(self.theta,r)
        for _ in range(24):
            den=lo[1]['reactivity_pcm']-hi[1]['reactivity_pcm']
            f=lo[1]['reactivity_pcm']/den if den>0 else .5
            theta=lo[0]+np.clip(f,.15,.85)*(hi[0]-lo[0])
            answer,r=self.evaluate(vec,power,float(theta))
            if accepted(r):
                answer,r=self.evaluate(vec,power,float(theta))
                if accepted(r):self.theta=float(theta);return answer
            if r['keff']>1:lo=(float(theta),r)
            else:hi=(float(theta),r)
        raise RuntimeError('Depletion drum controller failed its verified tolerance')


def interval(index, days, fuel_file, theta):
    d=(CAMPAIGN/'depletion'/f'interval_{index:04d}').resolve();d.mkdir(parents=True,exist_ok=True)
    done=d/'complete.json'
    request=dict(index=index,days=days,theta_deg=theta,power_W=400000.,fuel_sha256=hashlib.sha256(Path(fuel_file).read_bytes()).hexdigest() if fuel_file else 'fresh',chain_sha256=hashlib.sha256(CHAIN.read_bytes()).hexdigest())
    if done.exists():
        cached=json.loads(done.read_text())
        if cached['request']!=request:raise RuntimeError(f'Depletion checkpoint input mismatch: {d}')
        return cached
    (d/'request.json').write_text(json.dumps(request,indent=2)+'\n')
    mats={m.id:m for m in openmc.Materials.from_xml(fuel_file) if m.depletable} if fuel_file else None
    model,_=build(d,theta=theta,fuel_materials=mats,particles=20000,batches=200,inactive=80,target_sigma=2e-4,seed=32001+index*101)
    model.settings.trigger_max_batches=20000
    model.settings.export_to_xml(d/'settings.xml')
    os.environ['OPENMC_CROSS_SECTIONS']=str(DEFAULT_XS)
    os.environ['OMP_NUM_THREADS']='4'
    openmc.deplete.pool.USE_MULTIPROCESSING=False
    previous=Path.cwd();os.chdir(d)
    try:
        op=ControlledOperator(model,theta,d);op.directory_index=index
        integ=openmc.deplete.CECMIntegrator(op,[days],power=400000.,timestep_units='d')
        integ.integrate(final_step=True)
    finally:os.chdir(previous)
    results=openmc.deplete.Results(d/'depletion_results.h5')
    materials=results.export_to_materials(-1,path=d/'materials.xml')
    checkpoint=d/'materials_final.xml';materials.export_to_xml(checkpoint)
    inventory={n:sum(m.get_mass(n) for m in materials if m.depletable and n in m.get_nuclides())/1000 for n in ['U235','U238','Pu239','Pu240','Xe135','Sm149']}
    last=op.records[-1]
    result=dict(index=index,days=days,fuel_file=str(checkpoint),theta_deg=op.theta,keff=last['keff'],sigma=last['sigma'],inventory_kg=inventory,control_stage_count=op.stage,request=request)
    atomic_json(done,result);return result


def run():
    control=json.loads((CAMPAIGN/'neutronics/p1_control/control.json').read_text())
    dest=CAMPAIGN/'data/depletion_summary.json'
    if not control['critical']:
        result={'status':'NO_CRITICAL_BOL_STATE','reason':control['status'],'operable_lifetime_years':0.,'fuel_exhaustion_lifetime_years':None}
        atomic_json(dest,result);return result
    rows=[];elapsed=0.;fuel_file=None;theta=control['theta_deg'];index=0
    result={'status':'RUNNING','power_W':400000,'report_interval_years':2.5,'max_internal_step_years':.25,'rows':rows,'lifetime_bracket_years':None}
    startup=[.25,.75,2.,7.,20.,61.3125] # total 91.3125 d = .25 yr
    while True:
        days=startup[index] if index<len(startup) else 91.3125
        try:rec=interval(index,days,fuel_file,theta)
        except ControlUnavailable:
            result.update(status='CONTROL_LIMIT_REACHED',lifetime_bracket_years=[elapsed/365.25,(elapsed+days)/365.25],note='Bracket reflects the last completed interval and first failed CE/CM stage; not an exact lifetime or a propagated confidence interval.')
            atomic_json(dest,result);return result
        elapsed+=days;fuel_file=rec['fuel_file'];theta=rec['theta_deg']
        rec['time_years']=elapsed/365.25;rec['report_2_5yr']=abs(rec['time_years']/2.5-round(rec['time_years']/2.5))<1e-8
        rows.append(rec);atomic_json(dest,result);index+=1


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');a=ap.parse_args()
    if a.smoke:
        raise SystemExit('Use the campaign validation input for a short depletion API test; production smoke must not use synthetic compositions.')
    run()
