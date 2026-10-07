"""Durable, single-transport queue. Run in openmc-env; safe to restart."""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback
from p1_model import CAMPAIGN, ROOT, DEFAULT_XS
from p1_report import report

PYTHON='/home/aleri/miniforge3/envs/openmc-env/bin/python'
SCRIPTS=ROOT/'scripts'
STATE=CAMPAIGN/'status.json'


def save(state):
    state['updated_unix']=time.time()
    mem={line.split(':')[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemAvailable:','SwapFree:'))}
    state['available_memory_MiB']=round(mem['MemAvailable']/1024,1)
    state['free_swap_MiB']=round(mem['SwapFree']/1024,1)
    transport_file=CAMPAIGN/'current_transport.json'
    if transport_file.exists():
        try:
            transport=json.loads(transport_file.read_text())
            state['active_transport']=transport['name'] if Path(f"/proc/{transport.get('pid',-1)}/cmdline").exists() else None
        except json.JSONDecodeError:pass
    tmp=STATE.with_suffix('.json.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(STATE)
    report()


def stage(name):
    # Imports occur in the stage child; no OpenMC C-library state crosses stages.
    from p1_transport import run_point,feedback
    from p1_control import critical
    def result(n):return json.loads((CAMPAIGN/'neutronics'/n/'result.json').read_text())
    def ctrl(n):return json.loads((CAMPAIGN/f'neutronics/{n}_control/control.json').read_text())
    if name=='bol':
        for angle in (0.,180.):run_point(f'p1_bol_{angle:.0f}deg',theta=angle,seed=1201+int(angle))
        run_point('tb35_corrected_bol_0deg',reference=True,theta=0.,seed=1601)
    elif name=='source_check':
        r=run_point('p1_source_check',theta=0.,seed=83047,inactive=120,particles=20000)
        a=result('p1_bol_0deg');z=abs(a['keff']-r['keff'])/math.hypot(a['sigma'],r['sigma'])
        (CAMPAIGN/'data/source_check.json').write_text(json.dumps({'independent_keff_zscore':z,'passed_3sigma':z<=3.,'note':'Independent seed and longer inactive period; entropy files retained for source diagnostics.'},indent=2)+'\n')
        if z>3:raise RuntimeError(f'Independent source convergence comparison differs by {z:.2f} sigma')
    elif name=='p1_control':critical('p1_control',endpoints=[result('p1_bol_0deg'),result('p1_bol_180deg')])
    elif name=='tb35_control':
        r=run_point('tb35_corrected_bol_180deg',reference=True,theta=180.,seed=1781)
        critical('tb35_control',reference=True,endpoints=[result('tb35_corrected_bol_0deg'),r])
    elif name in ('p1_heating','tb35_heating'):
        which=name.split('_')[0];c=ctrl(which)
        # A noncritical normalized heating case is clearly marked in its record.
        theta=c['theta_deg'] if c['critical'] else (0. if c['status'].startswith('SUB') else 180.)
        r=run_point(name,theta=theta,reference=which=='tb35',photons=True,tallies=True,inactive=80,seed=90001 if which=='p1' else 90101)
        (CAMPAIGN/f'data/{name}_operating_status.json').write_text(json.dumps({'critical_operating_state':c['critical'],'theta_deg':theta,'description':'Nominal-temperature normalized heating; noncritical cases are hypothetical prescribed-power distributions.'},indent=2)+'\n')
    elif name in ('p1_thermal','tb35_thermal'):
        from p1_thermal import run
        which=name.split('_')[0];run(result(which+'_heating')['statepoint'],reference=which=='tb35')
    elif name=='feedback':
        c=ctrl('p1');theta=c['theta_deg'] if c['critical'] else 0.;feedback(theta)
    elif name=='depletion':
        from p1_depletion import run
        run()
    else:raise ValueError(name)


def main():
    CAMPAIGN.mkdir(exist_ok=True);lock=(CAMPAIGN/'campaign.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('Another campaign controller already owns the lock.')
    state=json.loads(STATE.read_text()) if STATE.exists() else {'completed':[],'errors':{}}
    state.update(status='RUNNING',controller_pid=os.getpid(),stage='waiting for initial BOL job')
    env=os.environ.copy();env.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENMC_CROSS_SECTIONS=str(DEFAULT_XS))
    # This first job was launched before the controller was written. Do not duplicate it.
    jobfile=CAMPAIGN/'bol_job.json'
    if jobfile.exists():
        old=json.loads(jobfile.read_text());pid=old['pid']
        while Path(f'/proc/{pid}/cmdline').exists():
            command=Path(f'/proc/{pid}/cmdline').read_bytes()
            if b'p1_transport.py' not in command:break
            save(state);time.sleep(15)
    for relative in ['validation/depletion_api/validation.json','validation/photon_heating/validation.json']:
        check=CAMPAIGN/relative
        if not check.exists() or not json.loads(check.read_text()).get('passed'):
            state.update(status='PREFLIGHT_FAILED',stage='validation')
            state['errors']['validation']='Missing successful numerical validation: '+str(check)
            save(state);return
    for which in ['p1','tb35']:
        check=CAMPAIGN/f'validation/thermal_synthetic/{which}/summary.json'
        validation=json.loads(check.read_text()) if check.exists() else {}
        if not validation.get('validation_only') or abs(validation.get('pin_mesh_change_peak_K',999.))>2. or abs(validation.get('cell_mesh_change_peak_K',999.))>2.:
            state.update(status='PREFLIGHT_FAILED',stage='validation')
            state['errors']['validation']='Missing successful thermal mesh/energy validation: '+str(check)
            save(state);return
    stages=['bol','source_check','p1_control','tb35_control','p1_heating','p1_thermal','tb35_heating','tb35_thermal','feedback','depletion']
    dependencies={'source_check':['bol'],'p1_control':['bol','source_check'],'tb35_control':['bol'],'p1_heating':['p1_control'],'tb35_heating':['tb35_control'],'p1_thermal':['p1_heating'],'tb35_thermal':['tb35_heating'],'feedback':['p1_control'],'depletion':['p1_control','source_check']}
    for name in stages:
        if name in state['completed']:continue
        missing=[n for n in dependencies.get(name,[]) if n not in state['completed']]
        if missing:state['errors'][name]='Skipped because dependencies failed: '+','.join(missing);save(state);continue
        state.update(stage=name,status='RUNNING');state['errors'].pop(name,None);save(state)
        logfile=CAMPAIGN/'logs'/f'{name}.log'
        command=[PYTHON,'-u',str(Path(__file__).resolve()),'--stage',name]
        print('START',name,flush=True)
        with logfile.open('a') as log:
            log.write(f'\nRESUME {time.time()}\n');log.flush()
            p=subprocess.Popen(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
            state.update(child_pid=p.pid,stage_log=str(logfile));save(state)
            while p.poll() is None:time.sleep(15);save(state)
        if p.returncode:
            state['errors'][name]=f'Exit {p.returncode}; see {logfile}';print('FAILED',name,flush=True)
        else:state['completed'].append(name);print('DONE',name,flush=True)
        save(state)
    state.update(status='FINISHED_WITH_ERRORS' if state['errors'] else 'COMPLETED',stage='finished',child_pid=None)
    save(state)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage');args=parser.parse_args()
    if args.stage:stage(args.stage)
    else:main()
