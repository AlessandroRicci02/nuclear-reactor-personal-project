#!/usr/bin/env python3
"""Generate genuine OpenMC geometry plots and independently sampled inventories."""
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path
import openmc
from p1_model import ROOT, CAMPAIGN, OPENMC, build


def main():
    d=CAMPAIGN/'validation/geometry_inventory'
    model,info=build(d,xs='/home/aleri/wmp_data/cross_sections_wmp.xml',legacy_data=True)
    colors={}
    palette={'UN':'#c85836','Nb1Zr':'#666b75','Al2O3':'#c6c8cc','YH2':'#efe2aa','BeO':'#91bac4','B4C':'#272938','Helium':'#76caef'}
    for m in model.materials:
        h=palette['UN' if m.name.startswith('UN') else m.name].lstrip('#')
        colors[m]=tuple(int(h[i:i+2],16) for i in (0,2,4))
    plots=[]
    specs=[('p1_radial_xy','xy',(0,0,0),(60,60),(2000,2000)),
           ('p1_axial_xz','xz',(0,0,0),(60,66),(1800,1980)),
           ('p1_inner_unit_cell','xy',(0,0,0),(2.8,2.8),(1400,1400))]
    for name,basis,origin,width,pixels in specs:
        plot=openmc.Plot();plot.filename=name;plot.basis=basis;plot.origin=origin
        plot.width=width;plot.pixels=pixels;plot.color_by='material';plot.colors=colors
        plots.append(plot)
    model.plots=openmc.Plots(plots);model.plots.export_to_xml(d/'plots.xml')
    with (d/'plot.log').open('w') as log:subprocess.run([OPENMC,'-p','-s','2'],cwd=d,stdout=log,stderr=subprocess.STDOUT,check=True)
    for name,*_ in specs:shutil.copy2(d/(name+'.png'),ROOT/'figure'/(name+'.png'))
    cfg=info['config'];R=cfg['metrics']['outer_radius_cm'];H=cfg['metrics']['total_height_cm']
    model.settings.run_mode='volume';model.settings.volume_calculations=[openmc.VolumeCalculation(list(model.materials),5_000_000,[-R,-R,-H/2],[R,R,H/2])]
    model.settings.export_to_xml(d/'settings.xml')
    with (d/'volume.log').open('w') as log:subprocess.run([OPENMC,'-s','4'],cwd=d,stdout=log,stderr=subprocess.STDOUT,check=True)
    vol=openmc.VolumeCalculation.from_hdf5(d/'volume_1.h5')
    groups={}
    for m in model.materials:
        name='UN fuel' if m.name.startswith('UN') else m.name
        if name not in groups:groups[name]={'mass':0.,'volume':0.,'density':m.get_mass_density()}
        groups[name]['volume']+=vol.volumes[m.id].n
    records=[]
    # Shared multinomial samples: compute group and total mass errors, accounting
    # for negative covariance between mutually exclusive material categories.
    domain=(2*R)**2*H;N=5_000_000
    for name,v in groups.items():
        prob=v['volume']/domain;mean=v['volume']*v['density']/1000
        sigma=domain/1000*v['density']*math.sqrt(prob*(1-prob)/N)
        records.append(dict(material=name,volume_cm3=v['volume'],density_g_cm3=v['density'],mass_kg=mean,sigma_kg=sigma))
    solids=[r for r in records if r['material']!='Helium']
    mean=sum(r['mass_kg'] for r in solids)
    erho=sum(r['volume_cm3']/domain*r['density_g_cm3'] for r in solids)
    erho2=sum(r['volume_cm3']/domain*r['density_g_cm3']**2 for r in solids)
    sigma=domain/1000*math.sqrt((erho2-erho**2)/N)
    with (CAMPAIGN/'data/mass_stochastic.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    result={'analytic_solids_kg':cfg['metrics']['modeled_solids_kg'],'stochastic_solids_kg':mean,'stochastic_sigma_kg':sigma,'samples':N,'helium_kg':groups['Helium']['volume']*groups['Helium']['density']/1000,'excludes':'Pressure vessel, drum drives/shrouds, end manifolds, supports, additional hydrogen encapsulation, shielding, conversion system'}
    result['analytic_stochastic_zscore']=(mean-result['analytic_solids_kg'])/sigma
    (CAMPAIGN/'data/mass_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    assert abs(result['analytic_stochastic_zscore'])<4,result
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
