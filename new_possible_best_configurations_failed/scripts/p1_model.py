#!/usr/bin/env python3
"""P1 transport model. Geometry from the reviewed JSON; units cm, g, kelvin.

Adapted from the historical TB35 builder. Explicit changes: three geometric zones,
elliptical auxiliary ducts, corrected disjoint endcap coolant union, actual YH2
thermal scattering, temperature interpolation and 30 fuel depletion regions.
Historical inputs are never imported for execution or changed.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
import openmc

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / 'campaign'
OPENMC = '/home/aleri/miniforge3/envs/openmc-env/bin/openmc'
CARDINAL = '/home/aleri/projects/cardinal/cardinal-opt'
DEFAULT_XS = CAMPAIGN / 'data/endfb-viii.0-hdf5/cross_sections.xml'
NOMINAL = dict(fuel=1195., clad=1170., moderator=904.5, reflector=870., coolant=1001.)


def load_config(reference=False):
    name = 'TB35_REFERENCE' if reference else 'P1_ELLIPTIC_BALANCED'
    config = json.loads((ROOT / 'configurations' / name / 'geometry.json').read_text())
    with (ROOT / 'configurations' / name / 'pins.csv').open() as f:
        pins = list(csv.DictReader(f))
    for p in pins:
        for k,v in p.items():
            p[k] = int(v) if k in ('pin_id','ring','lattice_a','lattice_b','zone') else float(v)
    return config, pins


def ellipse(a, b, radius, angle):
    c,s = math.cos(angle),math.sin(angle)
    x,y = radius*c,radius*s
    A,B,D = c*c/a**2+s*s/b**2, s*s/a**2+c*c/b**2, 2*c*s*(1/a**2-1/b**2)
    return openmc.Quadric(a=A,b=B,d=D,g=-2*A*x-D*y,h=-2*B*y-D*x,
                          k=A*x*x+B*y*y+D*x*y-1)


def build(directory, theta=0., reference=False, xs=None, temperatures=None,
          moderator_density=1., coolant_density=1., particles=20000, batches=130,
          inactive=40, seed=173, target_sigma=None, tallies=False, photons=False,
          fuel_materials=None, slices=10, legacy_data=False, grid=None):
    directory = Path(directory); directory.mkdir(parents=True,exist_ok=True)
    cfg,pins = load_config(reference)
    g,a,p = cfg['geometry'],cfg['auxiliary_channels'],cfg['parameters']
    rho = cfg['materials']['densities_g_cm3']
    temp = NOMINAL | (temperatures or {})
    xs = str(xs or DEFAULT_XS)
    openmc.reset_auto_ids(); openmc.config['cross_sections'] = xs
    H,R,Rout = p['height_cm'],p['core_radius_cm'],p['core_radius_cm']+p['radial_reflector_cm']
    totalH = H+2*p['axial_reflector_cm']
    dz = H/slices
    mats=[]
    def mat(mid,name,density,composition,T):
        m=openmc.Material(material_id=mid,name=name)
        m.set_density('g/cm3',density)
        for n,q,unit in composition:
            if n.isalpha():m.add_element(n,q,unit)
            else:m.add_nuclide(n,q,unit)
        m.temperature=float(T);mats.append(m);return m
    fuels={}
    # Reference geometry has the same spatial depletion discretization as P1.
    for z in cfg['zones']:
        zi=z['zone'];enr=z['uranium_u235_weight_fraction']
        for s in range(slices):
            mid=100+zi*slices+s
            if fuel_materials and mid in fuel_materials:
                m=fuel_materials[mid].clone();m.id=mid;mats.append(m)
                m.temperature=temp['fuel']
            else:
                m=mat(mid,f'UN_zone{zi}_slice{s}',rho['fuel'],
                      [('U235',enr*238/252,'wo'),('U238',(1-enr)*238/252,'wo'),('N14',14/252,'wo')],temp['fuel'])
            m.depletable=True;m.volume=z['pins']*math.pi*z['fuel_radius_cm']**2*dz
            fuels[(zi,s)]=m
    clad=mat(1,'Nb1Zr',rho['clad'],[('Nb',.99,'wo'),('Zr',.01,'wo')],temp['clad'])
    liner=mat(2,'Al2O3',rho['liner'],[('Al',2,'ao'),('O',3,'ao')],temp['moderator'])
    mod=mat(3,'YH2',rho['moderator']*moderator_density,[('Y89',1,'ao'),('H1',2,'ao')],temp['moderator'])
    mod.add_s_alpha_beta('c_H_in_ZrH' if legacy_data else 'c_H_in_YH2')
    if not legacy_data:mod.add_s_alpha_beta('c_Y_in_YH2')
    refl=mat(4,'BeO',rho['reflector'],[('Be9',1,'ao'),('O16',1,'ao')],temp['reflector'])
    refl.add_s_alpha_beta('c_Be_in_BeO');refl.add_s_alpha_beta('c_O_in_BeO')
    absorb=mat(5,'B4C',rho['absorber'],[('B10',3.6,'ao'),('B11',.4,'ao'),('C',1,'ao')],temp['reflector'])
    he=mat(6,'Helium',2.5e6*4.002602/(8314.462618*temp['coolant'])/1000*coolant_density,[('He4',1,'ao')],temp['coolant'])
    zplanes=[openmc.ZPlane(z0=float(v)) for v in np.linspace(-H/2,H/2,slices+1)]
    zb,zt=zplanes[0],zplanes[-1]
    zbot=openmc.ZPlane(z0=-totalH/2,boundary_type='vacuum')
    ztop=openmc.ZPlane(z0=totalH/2,boundary_type='vacuum')
    active=+zb & -zt;caps=(+zbot & -zb)|(+zt & -ztop)
    aux=[ellipse(a['a_cm'],a['b_cm'],a['centre_radius_cm'],math.radians(a['phase_deg']+j*360/a['count'])) for j in range(a['count'])]
    auxregion=openmc.Union([-e for e in aux])
    axialcool=(-openmc.ZCylinder(r=g['axial_main_duct_radius_cm']))|auxregion
    zone_surfs={z['zone']:[openmc.ZCylinder(r=z[key]) for key in ['fuel_radius_cm','clad_outer_radius_cm','annulus_outer_radius_cm','liner_outer_radius_cm']] for z in cfg['zones']}
    fuelcells=[];modcells=[];cladcells=[];linercells=[];metadata=[]
    universes={}
    for pin in pins:
        pi,zi=pin['pin_id'],pin['zone']
        f,c,ann,ln=zone_surfs[zi];cells=[]
        for s in range(slices):
            slab=+zplanes[s]&-zplanes[s+1]
            regions=[-f,+f&-c,+c&-ann,+ann&-ln,+ln&~auxregion,auxregion]
            fills=[fuels[(zi,s)],clad,he,liner,mod,he]
            names=['fuel','clad','main_coolant','liner','moderator','aux_coolant']
            cs=[openmc.Cell(name=f'{name}_p{pi}_s{s}',fill=m,region=reg&slab) for name,m,reg in zip(names,fills,regions)]
            if grid:
                for cell,name in zip(cs,names):
                    key='coolant' if 'coolant' in name else name
                    if key in grid:cell.temperature=float(grid[key][pi][s])
            cells.extend(cs);fuelcells.append(cs[0]);cladcells.append(cs[1]);linercells.append(cs[3]);modcells.append(cs[4])
            metadata.append(dict(pin_id=pi,zone=zi,ring=pin['ring'],slice=s,
                                 x_cm=pin['x_cm'],y_cm=pin['y_cm'],z_cm=-H/2+(s+.5)*dz,
                                 fuel_cell=cs[0].id,clad_cell=cs[1].id,liner_cell=cs[3].id,moderator_cell=cs[4].id))
        cells.extend([openmc.Cell(name=f'axial_duct_p{pi}',fill=he,region=axialcool&caps),
                      openmc.Cell(name=f'axial_reflector_p{pi}',fill=refl,region=~axialcool&caps)])
        universes[(pin['lattice_a'],pin['lattice_b'])]=openmc.Universe(cells=cells)
    outer=openmc.Universe(cells=[openmc.Cell(name='peripheral_moderator',fill=mod,region=active),
                                 openmc.Cell(fill=refl,region=caps)])
    lat=openmc.HexLattice(name='217_pin_y_orientation');lat.center=(0.,0.);lat.pitch=[g['pitch_cm']];lat.orientation='y'
    rings=[[outer]*(6*r) if r else [outer] for r in range(8,-1,-1)]
    lat.universes=rings
    for idx,u in universes.items():
        ir,iw=lat.get_universe_index(idx);rings[ir][iw]=u
    lat.universes=rings;lat.outer=outer
    corecyl=openmc.ZCylinder(r=R);outercyl=openmc.ZCylinder(r=Rout,boundary_type='vacuum')
    drumr=g['drum_radius_cm'];drumi=drumr-g['absorber_thickness_cm'];half=math.radians(g['absorber_sector_degrees']/2)
    dc=openmc.ZCylinder(r=drumr);di=openmc.ZCylinder(r=drumi)
    pl1=openmc.Plane(a=math.sin(half),b=-math.cos(half));pl2=openmc.Plane(a=math.sin(half),b=math.cos(half))
    absreg=-dc&+di&+pl1&+pl2
    drum_u=openmc.Universe(cells=[openmc.Cell(fill=absorb,region=absreg),openmc.Cell(fill=refl,region=-dc&~absreg)])
    drumcells=[];drumregions=[]
    angles=[float(theta)]*10 if np.isscalar(theta) else list(theta)
    for j in range(10):
        az=j*36;phi=math.radians(az);x=g['drum_centres_radius_cm']*math.cos(phi);y=g['drum_centres_radius_cm']*math.sin(phi)
        cr=openmc.ZCylinder(x0=x,y0=y,r=drumr)
        reg=-cr&active;cell=openmc.Cell(name=f'drum_{j}',fill=drum_u,region=reg)
        cell.translation=(x,y,0.);cell.rotation=(0.,0.,az+angles[j])
        drumcells.append(cell);drumregions.append(reg)
    main=[openmc.Cell(name='core_lattice',fill=lat,region=-corecyl&+zbot&-ztop),
          openmc.Cell(name='radial_reflector',fill=refl,region=+corecyl&-outercyl&active&~openmc.Union(drumregions)),
          openmc.Cell(name='outer_axial_reflector',fill=refl,region=+corecyl&-outercyl&caps)]
    geometry=openmc.Geometry(main+drumcells)
    settings=openmc.Settings();settings.particles=particles;settings.batches=batches;settings.inactive=inactive;settings.seed=seed
    settings.temperature=({'method':'nearest','multipole':True,'tolerance':1500.} if legacy_data else {'method':'interpolation','tolerance':10.})
    settings.photon_transport=photons
    settings.source=openmc.IndependentSource(space=openmc.stats.Box((-R,-R,-H/2),(R,R,H/2)),constraints={'fissionable':True})
    mesh=openmc.RegularMesh();mesh.lower_left=(-R,-R,-H/2);mesh.upper_right=(R,R,H/2);mesh.dimension=(10,10,10);settings.entropy_mesh=mesh
    settings.output={'tallies':False}
    if target_sigma:
        settings.keff_trigger={'type':'std_dev','threshold':target_sigma}
        settings.trigger_active=True;settings.trigger_max_batches=max(2000,batches);settings.trigger_batch_interval=50
    tlist=[]
    if tallies:
        for name,cells in [('fuel',fuelcells),('clad',cladcells),('liner',linercells),('moderator',modcells)]:
            t=openmc.Tally(name=f'{name}_heating');t.filters=[openmc.CellFilter(cells)];t.scores=['heating' if photons else 'heating-local'];tlist.append(t)
        tm=openmc.Tally(name='material_heating');tm.filters=[openmc.MaterialFilter(mats)];tm.scores=['heating' if photons else 'heating-local'];tlist.append(tm)
        tq=openmc.Tally(name='recoverable_power');tq.scores=['fission-q-recoverable'];tlist.append(tq)
        ts=openmc.Tally(name='fuel_spectrum');ts.filters=[openmc.MaterialFilter(list(fuels.values())),openmc.EnergyFilter([0.,.625,1e5,20e6])];ts.scores=['flux'];tlist.append(ts)
    materials=openmc.Materials(mats);materials.cross_sections=xs
    model=openmc.Model(geometry,materials,settings,openmc.Tallies(tlist))
    model.export_to_xml(directory)
    manifest=dict(theta_deg=angles,reference=reference,cross_sections=xs,temperatures=temp,moderator_density=moderator_density,coolant_density=coolant_density,particles=particles,batches=batches,inactive=inactive,seed=seed,target_sigma=target_sigma,photons=photons,legacy_data=legacy_data,slices=slices,geometry_source=str(ROOT/'configurations'/('TB35_REFERENCE' if reference else 'P1_ELLIPTIC_BALANCED')/'geometry.json'))
    (directory/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    with (directory/'cell_metadata.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=metadata[0]);w.writeheader();w.writerows(metadata)
    return model,dict(config=cfg,pins=pins,fuel_cells=fuelcells,moderator_cells=modcells,fuel_materials=fuels,metadata=metadata,materials={m.name:m for m in mats})


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('directory');ap.add_argument('--theta',type=float,default=0);ap.add_argument('--reference',action='store_true');ap.add_argument('--xs');ap.add_argument('--legacy-data',action='store_true');ap.add_argument('--particles',type=int,default=20000);ap.add_argument('--batches',type=int,default=130);ap.add_argument('--inactive',type=int,default=40);ap.add_argument('--tallies',action='store_true');ap.add_argument('--photons',action='store_true');args=vars(ap.parse_args())
    build(**args)
