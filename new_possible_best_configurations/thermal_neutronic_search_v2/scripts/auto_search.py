"""Restartable finite campaign: reject immediately; stop on first full pass."""
from pathlib import Path
import csv
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import time
import traceback
import numpy as np
import openmc
from catalog import ROOT,catalog,layers,pins
from model import build,initial_field,OPENMC,XS
from thermal import extract,evaluate

Z=4.0
THREADS=int(os.environ.get('NEP_SEARCH_THREADS','6'))
VERSION='2026-10-08-cold-first-v2'
SCRIPT_HASH=hashlib.sha256(b''.join((ROOT/'scripts'/x).read_bytes() for x in
    ('model.py','thermal.py','catalog.py','auto_search.py'))).hexdigest()
XS_INDEX_HASH=hashlib.sha256(XS.read_bytes()).hexdigest()


def atomic(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)


class Reject(Exception):
    def __init__(self,metric,detail):self.metric=metric;self.detail=detail


def transport(candidate,directory,field,theta,*,target=.00035,heating=False,
              particles=16000,perturb_kind=None,delta_K=0.,beta_mod=0.,beta_fuel=0.):
    d=Path(directory);d.mkdir(parents=True,exist_ok=True)
    request=dict(candidate=candidate['id'],case_id=str(d.relative_to(ROOT)),spec=candidate['spec'],field=field,theta=theta,
        target=target,heating=heating,particles=particles,perturb_kind=perturb_kind,
        delta_K=delta_K,beta_mod=beta_mod,beta_fuel=beta_fuel,script_hash=SCRIPT_HASH)
    request['cross_sections_index_sha256']=XS_INDEX_HASH
    sig=hashlib.sha256(json.dumps(request,sort_keys=True).encode()).hexdigest()
    if (d/'result.json').exists():
        r=json.loads((d/'result.json').read_text())
        if r['request_sha256']!=sig:raise RuntimeError(f'Cached request mismatch: {d}')
        return r,None
    seed=1+int(sig[:8],16)%1000000000
    model,info=build(d,candidate['spec'],field=field,theta=theta,particles=particles,
        batches=200 if heating else 160,inactive=50 if heating else 40,
        seed=seed,heating=heating,target=target,
        max_batches=12000 if target and target<=.000075 else 4000,
        perturb_kind=perturb_kind,delta_K=delta_K,beta_mod=beta_mod,beta_fuel=beta_fuel)
    atomic(d/'request.json',request|dict(request_sha256=sig,seed=seed))
    atomic(ROOT/'status.json',dict(status='running',candidate=candidate['id'],case=str(d.relative_to(ROOT)),
        stage='neutron_photon_heating' if heating else 'eigenvalue',started_unix=time.time(),pid=os.getpid()))
    env=os.environ.copy();env.update(OPENMC_CROSS_SECTIONS=str(XS),OMP_NUM_THREADS=str(THREADS),
        OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with (d/'run.log').open('w') as log:
        cp=subprocess.run([OPENMC,'-s',str(THREADS)],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
    if cp.returncode:raise RuntimeError(f'OpenMC returned {cp.returncode}: {d}/run.log')
    state=max(d.glob('statepoint.*.h5'),key=lambda p:int(p.name.split('.')[1]))
    with openmc.StatePoint(state,autolink=False) as sp:
        r=dict(request_sha256=sig,keff=float(sp.keff.n),sigma_k=float(sp.keff.s),
            batches=int(sp.current_batch),particles=int(sp.n_particles),statepoint=str(state),theta=theta)
    r['lower_bound']=r['keff']-Z*r['sigma_k'];r['upper_bound']=r['keff']+Z*r['sigma_k']
    atomic(d/'result.json',r)
    print(candidate['id'],d.name,r['keff'],'+/-',r['sigma_k'],flush=True)
    return r,info if heating else None


def endpoint(candidate,directory,field,theta,pilot=False):
    result=None
    targets=[.0010,.0004,.0002] if pilot else [.0003,.00015,.000075]
    for target in targets:
        d=directory/f'sigma_{target:.6f}'
        result,_=transport(candidate,d,field,theta,target=target,particles=8000 if pilot else 20000)
        ok=result['lower_bound']>1 if theta==0 else result['upper_bound']<1
        bad=result['upper_bound']<1 if theta==0 else result['lower_bound']>1
        if ok:return result
        if bad:raise Reject('k_out' if theta==0 else 'k_in',result)
    raise Reject('unresolved_k_out' if theta==0 else 'unresolved_k_in',result)


def critical_angle(candidate,directory,field,k0,k180,pilot=False):
    low=(0.,k0['keff']);high=(180.,k180['keff']);points=[]
    target=.0006 if pilot else .0003;tolerance=.0012 if pilot else .0005
    for j in range(10):
        theta=low[0]+(1-low[1])*(high[0]-low[0])/(high[1]-low[1])
        if not low[0]+.01<theta<high[0]-.01:theta=(low[0]+high[0])/2
        theta=round(theta,3)
        r,_=transport(candidate,directory/f'p{j:02d}_{theta:.3f}',field,theta,target=target,particles=12000)
        points.append(r)
        if abs(r['keff']-1)<tolerance:return r
        if r['keff']>1:low=(theta,r['keff'])
        else:high=(theta,r['keff'])
    raise Reject('unresolved_critical_angle',dict(points=points))


def feedback(candidate,directory,field,theta,kind,bm,bf,delta_override=None):
    delta=delta_override or (100. if kind=='fuel' else 50.)
    for target in (.00025,.00012,.00006):
        rr=[]
        for sign in (-1,1):
            r,_=transport(candidate,directory/f'sigma_{target:.6f}'/('minus' if sign<0 else 'plus'),
                field,theta,target=target,particles=20000,perturb_kind=kind,
                delta_K=sign*delta,beta_mod=bm,beta_fuel=bf)
            rr.append(r)
        lo,hi=rr
        rho=lambda r:1e5*(1-1/r['keff'])
        srho=lambda r:1e5*r['sigma_k']/r['keff']**2
        alpha=(rho(hi)-rho(lo))/(2*delta)
        sigma=math.hypot(srho(lo),srho(hi))/(2*delta)
        out=dict(kind=kind,alpha_pcm_K=alpha,sigma_pcm_K=sigma,
            upper_bound_pcm_K=alpha+Z*sigma,lower_bound_pcm_K=alpha-Z*sigma,
            beta_mod_per_K=bm,beta_fuel_per_K=bf,delta_K=delta,
            confidence_Z=Z,minus=lo,plus=hi)
        atomic(directory/'coefficient.json',out)
        if out['upper_bound_pcm_K']<0:return out
        if out['lower_bound_pcm_K']>=0:raise Reject('alpha_'+kind,out)
    raise Reject('unresolved_alpha_'+kind,out)


def save_progress():
    rows=[]
    for c in catalog():
        p=ROOT/'campaign'/c['id']/'assessment.json'
        a=json.loads(p.read_text()) if p.exists() else {}
        rows.append(dict(id=c['id'],status=a.get('status','pending'),
            failed_metric=a.get('failed_metric',''),detail=json.dumps(a.get('detail',{}),separators=(',',':'))))
    tmp=ROOT/'progress.csv.tmp'
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    tmp.replace(ROOT/'progress.csv')


def study(candidate):
    d=ROOT/'campaign'/candidate['id'];d.mkdir(exist_ok=True)
    atomic(d/'candidate.json',candidate)
    field=initial_field(candidate['spec'])
    pilotout=endpoint(candidate,d/'pilot/out',field,0.,pilot=True)
    pilotin=endpoint(candidate,d/'pilot/in',field,180.,pilot=True)
    critical=critical_angle(candidate,d/'pilot/critical',field,pilotout,pilotin,pilot=True)
    thermal=None;finalout=None;finalin=None
    for iteration in range(3):
        it=d/f'operating_state_{iteration}'
        theta=critical['theta']
        heat,info=transport(candidate,it/'heating',field,theta,target=None,heating=True,particles=15000)
        hd=it/'heating'
        if (hd/'heating_W.npz').exists():
            arrays=np.load(hd/'heating_W.npz');maps={k:arrays[k] for k in arrays.files if not k.startswith('sigma_')}
            errors={k:arrays['sigma_'+k] for k in maps};dep=json.loads((hd/'deposition.json').read_text())
        else:
            if info is None:
                _,info=build(hd,candidate['spec'],field=field,theta=theta,heating=True)
            maps,errors,dep=extract(hd,info,candidate['spec'])
        atomic(ROOT/'status.json',dict(status='running',candidate=candidate['id'],stage='thermal_Cardinal',iteration=iteration,pid=os.getpid()))
        thermal=evaluate(it/'thermal',candidate['spec'],maps,errors,dep)
        atomic(it/'thermal_result.json',thermal)
        if not thermal['pass_thermal']:raise Reject('T_moderator_upper',dict(upper_K=thermal['acceptance_upper_K'],result=str(it/'thermal_result.json')))
        newfield=thermal['nominal']['field']
        finalout=endpoint(candidate,it/'verified_out',newfield,0.)
        finalin=endpoint(candidate,it/'verified_in',newfield,180.)
        newcritical=critical_angle(candidate,it/'verified_critical',newfield,finalout,finalin)
        error=max(float(np.max(np.abs(np.asarray(newfield[k])-np.asarray(field[k])))) for k in ('moderator','hot_he','cold_he'))
        fuel_error=float(np.max(np.abs(np.asarray(newfield['fuel'])-np.asarray(field['fuel']))/np.asarray(field['fuel'])))
        angle_error=abs(newcritical['theta']-theta)
        atomic(it/'convergence.json',dict(max_coolant_moderator_change_K=error,max_fuel_relative_change=fuel_error,angle_change_deg=angle_error))
        field=newfield;critical=newcritical
        if iteration>=1 and error<8. and fuel_error<.05 and angle_error<3.:break
    else:raise Reject('unresolved_operating_state',dict(max_field_change_K=error,angle_change_deg=angle_error))
    coeffs=[]
    # Mass-conserving sleeve radial expansion into explicit clearance; fixed
    # axial/inner boundaries. Bracket material expansion rather than choosing
    # the expansion law that makes the coefficient look most favorable.
    bmlo=2.4e-5;bmhi=4.95e-5;bflo=2.1e-5;bfhi=3.9e-5
    combinations=[('moderator',bmlo,0.),('moderator',(bmlo+bmhi)/2,0.),('moderator',bmhi,0.),
        ('fuel',0.,bflo),('fuel',0.,(bflo+bfhi)/2),('fuel',0.,bfhi)]
    combinations += [('total',bm,bf) for bm in (bmlo,bmhi) for bf in (bflo,bfhi)]
    for j,(kind,bm,bf) in enumerate(combinations):
        out=feedback(candidate,d/f'feedback/{j:02d}_{kind}',field,critical['theta'],kind,bm,bf)
        coeffs.append(out);atomic(d/'feedback_results.json',coeffs)
    # Independent short-step confirmation prevents a favorable wide finite
    # difference from silently qualifying a non-local moderator/total slope.
    confirmation=[]
    for kind,bm,bf in [('moderator',bmlo,0.),('fuel',0.,bflo),('total',bmlo,bflo)]:
        confirmation.append(feedback(candidate,d/f'confirmation/{kind}',field,critical['theta'],kind,bm,bf,
            delta_override=50. if kind=='fuel' else 25.))
    return dict(status='qualified_for_screened_BOL_conditions',candidate=candidate,
        out=finalout,inward=finalin,critical=critical,thermal=thermal,
        feedback=coeffs,confirmation=confirmation,field=field,
        acceptance_scope='4 sigma MC eigenvalue/coefficient bounds; bounded thermo/material expansion screen. Not reactor licensing or lifespan validation.')


def publish(result):
    dest=ROOT/'found_reactor'/result['candidate']['id'];dest.mkdir(parents=True,exist_ok=True)
    source=ROOT/'campaign'/result['candidate']['id']
    atomic(dest/'results.json',result);atomic(dest/'configuration.json',result['candidate'])
    (dest/'simulation_directory.txt').write_text(str(source)+'\n')
    # No further transport after qualifying. Retain complete evidence in the
    # campaign and only the compact results/figures in this presentation folder.
    from report import write_found
    write_found(dest,result)
    atomic(ROOT/'solution.json',dict(candidate=result['candidate']['id'],folder=str(dest),result=str(dest/'results.json')))


def main():
    lock=(ROOT/'controller.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('Another controller is running')
    if (ROOT/'solution.json').exists():print('Already found; no simulations launched');return
    atomic(ROOT/'job.json',dict(pid=os.getpid(),started_unix=time.time(),version=VERSION,threads=THREADS,
        catalog_count=30,script_hash=SCRIPT_HASH,cross_sections_index_sha256=XS_INDEX_HASH,
        openmc_version=openmc.__version__,no_depletion=True))
    save_progress()
    for c in catalog():
        d=ROOT/'campaign'/c['id'];d.mkdir(exist_ok=True)
        ap=d/'assessment.json'
        if ap.exists():
            previous=json.loads(ap.read_text())
            if previous.get('status')=='rejected':continue
        try:
            result=study(c)
        except Reject as exc:
            atomic(ap,dict(status='rejected',failed_metric=exc.metric,detail=exc.detail,completed_unix=time.time()))
            print('REJECTED',c['id'],exc.metric,flush=True);save_progress();continue
        atomic(ap,result);publish(result);save_progress()
        atomic(ROOT/'status.json',dict(status='found_and_stopped',candidate=c['id'],folder=str(ROOT/'found_reactor'/c['id']),completed_unix=time.time()))
        print('FOUND; ALL SIMULATIONS STOPPED',c['id'],flush=True);return
    save_progress();atomic(ROOT/'status.json',dict(status='exhausted_no_solution',tested_candidates=30,completed_unix=time.time()))
    print('All 30 exhausted; no qualifying reactor found',flush=True)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        atomic(ROOT/'status.json',dict(status='numerical_failure',error=str(exc),failed_unix=time.time()))
        traceback.print_exc();raise
