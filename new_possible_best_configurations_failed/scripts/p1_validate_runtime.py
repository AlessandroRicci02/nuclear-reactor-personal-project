"""Reproduce the small numerical API checks; never used as reactor results.

env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 openmc-env/bin/python
    scripts/p1_validate_runtime.py {depletion,photons}
The full-chain check needs about 6 GiB, despite its short particle history.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import openmc
from p1_model import build,CAMPAIGN,DEFAULT_XS,OPENMC


def depletion():
    import openmc.deplete
    import openmc.lib
    d=(CAMPAIGN/'validation/depletion_api').resolve();d.mkdir(parents=True,exist_ok=True)
    model,_=build(d,particles=500,batches=16,inactive=8,seed=9453)
    openmc.deplete.pool.USE_MULTIPROCESSING=False
    class AuditOperator(openmc.deplete.CoupledOperator):
        def __call__(self,vec,power):
            for cell in self.model.geometry.get_all_cells().values():
                if cell.name.startswith('drum_'):
                    az=int(cell.name.split('_')[1])*36.
                    openmc.lib.cells[cell.id].rotation=(0.,0.,az+90.)
            return super().__call__(vec,power)
    old=Path.cwd();os.chdir(d)
    try:
        op=AuditOperator(model,str(CAMPAIGN/'data/chain_endfb80_thermal.xml'),normalization_mode='fission-q')
        integrator=openmc.deplete.CECMIntegrator(op,[.001],power=400000.,timestep_units='d')
        integrator.integrate()
        results=openmc.deplete.Results('depletion_results.h5')
        mats=results.export_to_materials(-1,path='materials.xml');mats.export_to_xml('materials_final.xml')
        assert len([m for m in mats if m.depletable])==30
        result={'validation_only':True,'purpose':'CE/CM full-chain API, 30 compositions, C API drum rotation, final-material export. Deliberately low statistics; not a critical reactor result.','timesteps_days':[.001],'depletable_materials':30,'chain_nuclides':len(op.chain.nuclides),'passed':True}
        (d/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    finally:os.chdir(old)


def photons():
    from p1_thermal import heating
    d=(CAMPAIGN/'validation/photon_heating').resolve();d.mkdir(parents=True,exist_ok=True)
    build(d,theta=110.,photons=True,tallies=True,particles=1000,batches=20,inactive=10,seed=96353)
    with (d/'run.log').open('w') as f:subprocess.run([OPENMC,'-s','1'],cwd=d,stdout=f,stderr=subprocess.STDOUT,check=True)
    arrays,_=heating(d/'statepoint.20.h5',d)
    assert all(a.shape==(217,10) and a.sum()>0 for a in arrays.values())
    result={'validation_only':True,'passed':True,'purpose':'Low-statistics neutron/photon transport and material/cell heating handoff; not a precision-qualified reactor result.','array_shapes':{k:list(v.shape) for k,v in arrays.items()}}
    (d/'validation.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('task',choices=['depletion','photons']);args=parser.parse_args()
    os.environ['OPENMC_CROSS_SECTIONS']=str(DEFAULT_XS);os.environ['OMP_NUM_THREADS']='1'
    (depletion if args.task=='depletion' else photons)()
