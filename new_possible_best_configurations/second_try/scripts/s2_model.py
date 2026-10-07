"""OpenMC CSG of S2-R1; centimetres, Kelvin, grams. No implicit hardware."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import openmc

ROOT = Path(__file__).resolve().parents[1]
XS = Path('/home/aleri/OpenMC_projects/first_test/new_possible_best_configurations_failed/campaign/data/endfb-viii.0-hdf5/cross_sections.xml')
TEMP = dict(fuel=1300., clad=1200., hot_he=1050., hot_wall=1120., foil=1000.,
            cold_wall=875., cold_he=825., inner_can=875., barrier=875.,
            moderator=875., outer_can=880., be=850., absorber=850.)


def build(directory, *, temperatures=None, theta=0., particles=2500,
          batches=110, inactive=40, seed=601, heating=False, slices=6,
          volume=False, photons=False):
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((ROOT/'configuration.json').read_text());p=cfg['parameters']
    with (ROOT/'pins.csv').open() as f:pins=list(csv.DictReader(f))
    openmc.reset_auto_ids();openmc.config['cross_sections']=str(XS)
    T=TEMP | (temperatures or {})
    mats=[]
    def mat(name,rho,items,key,sabs=()):
        m=openmc.Material(name=name);m.set_density('g/cm3',rho)
        for nuclide,q,unit in items:
            if nuclide in ('Al','O','Mo','Zr','C'):
                m.add_element(nuclide,q,unit)
            else:m.add_nuclide(nuclide,q,unit)
        for sab in sabs:m.add_s_alpha_beta(sab)
        m.temperature=T[key];mats.append(m);return m
    fuel={}
    for zone,enr in [('inner',.196),('outer',.1975)]:
        fuel[zone]=mat('UN_'+zone,13.5,[('U235',enr*238/252,'wo'),
            ('U238',(1-enr)*238/252,'wo'),('N14',14/252,'wo')],'fuel')
    clad=mat('Nb1Zr_clad',8.57,[('Nb93',.99,'wo'),('Zr',.01,'wo')],'clad')
    hw=mat('Nb1Zr_hot_wall',8.57,[('Nb93',.99,'wo'),('Zr',.01,'wo')],'hot_wall')
    cw=mat('Nb1Zr_cold_wall',8.57,[('Nb93',.99,'wo'),('Zr',.01,'wo')],'cold_wall')
    ic=mat('Nb1Zr_inner_can',8.57,[('Nb93',.99,'wo'),('Zr',.01,'wo')],'inner_can')
    oc=mat('Nb1Zr_outer_can',8.57,[('Nb93',.99,'wo'),('Zr',.01,'wo')],'outer_can')
    foil=mat('Mo_shield',10.28,[('Mo',1.,'ao')],'foil')
    barrier=mat('Al2O3_barrier',3.95,[('Al',2.,'ao'),('O',3.,'ao')],'barrier')
    mod=mat('YH2',4.3,[('Y89',1.,'ao'),('H1',2.,'ao')],'moderator',
            ('c_H_in_YH2','c_Y_in_YH2'))
    be=mat('Be_reflector',1.85,[('Be9',1.,'ao')],'be',('c_Be',))
    b4c=mat('B4C',2.52,[('B10',3.6,'ao'),('B11',.4,'ao'),('C',1.,'ao')],'absorber')
    def he(name,key):
        rho=p['pressure_Pa']*.004002602/(8.314462618*T[key])/1000
        return mat(name,rho,[('He4',1.,'ao')],key)
    hot=he('He_hot', 'hot_he');cold=he('He_cold','cold_he')
    layers=cfg['radial_layers']; radii=[v['r_outer_cm'] for v in layers]
    assert radii==sorted(radii)
    zbot=openmc.ZPlane(z0=-(p['active_height_cm']/2+p['axial_reflector_cm']),boundary_type='vacuum')
    ztop=openmc.ZPlane(z0=+(p['active_height_cm']/2+p['axial_reflector_cm']),boundary_type='vacuum')
    zb=openmc.ZPlane(z0=-p['active_height_cm']/2);zt=openmc.ZPlane(z0=p['active_height_cm']/2)
    active=+zb & -zt;end=(+zbot & -zb)|(+zt & -ztop)
    r=[openmc.ZCylinder(r=x) for x in radii]
    shell=lambda i: (-r[i] if i==0 else +r[i-1]&-r[i])
    av=[fuel,clad,hot,hw,None,foil,None,cw,cold,ic,barrier,mod,barrier,oc]
    # Axial continuation: hot bore, jacket, cold annulus, inner can, then Be.
    ev=[hot,hw,None,foil,None,cw,cold,ic]
    axial_end_indices=[2,3,4,5,6,7,8,9]
    cells=[];pin_regions=[];fuelcells=[];modcells=[];hotcells=[];coldcells=[]
    layercells={i:[] for i in (3,5,7,9,10,12,13)}
    splits=[openmc.ZPlane(z0=-p['active_height_cm']/2+j*p['active_height_cm']/slices) for j in range(slices+1)]
    for pin in pins:
        idx=int(pin['id']);zone='inner' if int(pin['ring'])<=2 else 'outer'
        local=[]
        for j in range(slices):
            slab=+splits[j]&-splits[j+1]
            fc=openmc.Cell(name=f'fuel_p{idx}_z{j}',fill=fuel[zone],region=shell(0)&slab)
            mc=openmc.Cell(name=f'YH2_p{idx}_z{j}',fill=mod,region=shell(11)&slab)
            local.extend((fc,mc));fuelcells.append(fc);modcells.append(mc)
        for i,material in enumerate(av):
            if i in (0,11):continue
            c=openmc.Cell(name=f'{layers[i]["name"]}_p{idx}',fill=material,region=shell(i)&active)
            local.append(c)
            if i==2:hotcells.append(c)
            if i==8:coldcells.append(c)
            if i in layercells:layercells[i].append(c)
        local.append(openmc.Cell(name=f'core_Be_remainder_p{idx}',fill=be,region=+r[-1]&active))
        local.append(openmc.Cell(name=f'end_hot_bore_p{idx}',fill=hot,region=-r[2]&end))
        for i,material in zip(axial_end_indices[1:],ev[1:]):
            local.append(openmc.Cell(name=f'end_{layers[i]["name"]}_p{idx}',fill=material,region=shell(i)&end))
        local.append(openmc.Cell(name=f'end_Be_p{idx}',fill=be,region=+r[9]&end))
        u=openmc.Universe(cells=local)
        x=float(pin['x_cm']);y=float(pin['y_cm'])
        bounding=openmc.ZCylinder(x0=x,y0=y,r=radii[-1])
        reg=-bounding&+zbot&-ztop
        c=openmc.Cell(name=f'pin_{idx}',fill=u,region=reg)
        c.translation=(x,y,0.)
        cells.append(c);pin_regions.append(-bounding)
    core=openmc.ZCylinder(r=p['core_radius_cm'])
    outer=openmc.ZCylinder(r=p['core_radius_cm']+p['radial_reflector_cm'],boundary_type='vacuum')
    core_be_cell=openmc.Cell(name='core_Be_between_pins',fill=be,
                            region=-core&~openmc.Union(pin_regions)&+zbot&-ztop)
    cells.append(core_be_cell)
    dr=p['drum_radius_cm'];di=dr-p['absorber_thickness_cm'];half=math.radians(p['absorber_sector_deg']/2)
    dc=openmc.ZCylinder(r=dr);din=openmc.ZCylinder(r=di)
    pl1=openmc.Plane(a=math.sin(half),b=-math.cos(half));pl2=openmc.Plane(a=math.sin(half),b=math.cos(half))
    absreg=-dc&+din&+pl1&+pl2
    du=openmc.Universe(cells=[openmc.Cell(fill=b4c,region=absreg),openmc.Cell(fill=be,region=-dc&~absreg)])
    dregs=[]
    for k in range(p['drum_count']):
        az=360*k/p['drum_count'];phi=math.radians(az)
        x=p['drum_centres_radius_cm']*math.cos(phi);y=p['drum_centres_radius_cm']*math.sin(phi)
        cyl=openmc.ZCylinder(x0=x,y0=y,r=dr);reg=-cyl&active
        c=openmc.Cell(name=f'drum_{k}',fill=du,region=reg)
        c.translation=(x,y,0.);c.rotation=(0.,0.,az+theta)
        cells.append(c);dregs.append(reg)
    cells.append(openmc.Cell(name='radial_Be',fill=be,
                             region=+core&-outer&active&~openmc.Union(dregs)))
    cells.append(openmc.Cell(name='axial_outer_Be',fill=be,region=+core&-outer&end))
    settings=openmc.Settings();settings.particles=particles;settings.batches=batches;settings.inactive=inactive;settings.seed=seed
    settings.source=openmc.IndependentSource(space=openmc.stats.Box((-p['core_radius_cm'],)*2+(-p['active_height_cm']/2,),
             (p['core_radius_cm'],)*2+(p['active_height_cm']/2,)),constraints={'fissionable':True})
    settings.temperature={'method':'interpolation','tolerance':10.}
    settings.photon_transport=photons;settings.output={'tallies':heating}
    em=openmc.RegularMesh();em.lower_left=(-p['core_radius_cm'],-p['core_radius_cm'],-p['active_height_cm']/2)
    em.upper_right=(p['core_radius_cm'],p['core_radius_cm'],p['active_height_cm']/2);em.dimension=(8,8,8)
    settings.entropy_mesh=em
    tlist=[]
    if heating:
        score='heating' if photons else 'heating-local'
        for name,subset in [('fuel',fuelcells),('moderator',modcells),('hot_He',hotcells),('cold_He',coldcells)]:
            t=openmc.Tally(name=name+'_heating');t.filters=[openmc.CellFilter(subset)];t.scores=[score];tlist.append(t)
        for i,subset in layercells.items():
            t=openmc.Tally(name=f'layer_{i}_heating');t.filters=[openmc.CellFilter(subset)];t.scores=[score];tlist.append(t)
        t=openmc.Tally(name='core_Be_heating');t.filters=[openmc.CellFilter([core_be_cell])];t.scores=[score];tlist.append(t)
        t=openmc.Tally(name='all_material_heating');t.filters=[openmc.MaterialFilter(mats)];t.scores=[score];tlist.append(t)
    materials=openmc.Materials(mats);materials.cross_sections=str(XS)
    model=openmc.Model(openmc.Geometry(cells),materials,settings,openmc.Tallies(tlist))
    model.export_to_xml(directory)
    (directory/'manifest.json').write_text(json.dumps(dict(config=str(ROOT/'configuration.json'),
        cross_sections=str(XS),theta=theta,temperatures=T,particles=particles,
        batches=batches,inactive=inactive,seed=seed,slices=slices,
        heating=heating,photons=photons,volume=volume),indent=2)+'\n')
    if volume:
        model.settings.run_mode='volume'
        model.settings.volume_calculations=[openmc.VolumeCalculation(list(materials),2_000_000,
            [-p['core_radius_cm']-p['radial_reflector_cm']]*2+[-(p['active_height_cm']/2+p['axial_reflector_cm'])],
            [p['core_radius_cm']+p['radial_reflector_cm']]*2+[(p['active_height_cm']/2+p['axial_reflector_cm'])])]
        model.settings.export_to_xml(directory/'settings.xml')
    return model,dict(fuelcells=fuelcells,modcells=modcells,mats=mats,pins=pins,config=cfg)
