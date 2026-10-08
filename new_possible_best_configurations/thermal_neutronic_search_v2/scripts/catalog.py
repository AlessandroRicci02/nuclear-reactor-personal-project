"""Thirty deterministic, thermally separated YH2/UN candidate designs."""
from pathlib import Path
import csv
import json
import math

ROOT = Path(__file__).resolve().parents[1]


def pins(spec):
    result = []
    n = spec['rings']
    for q in range(-n, n+1):
        for r in range(-n, n+1):
            ring = max(abs(q), abs(r), abs(q+r))
            if ring <= n:
                result.append(dict(id=len(result), ring=ring,q=q,r=r,
                    x_cm=spec['pitch_cm']*(q+r/2),
                    y_cm=spec['pitch_cm']*math.sqrt(3)*r/2))
    return result


def layers(spec):
    # The sleeve has explicit outer expansion clearance. Gaps/cans are scored
    # in neutron/photon transport instead of using invisible insulation.
    definitions = [('fuel',spec['fuel_radius_cm'],'UN',spec['fuel_density_g_cm3']),
        ('fuel_gap',.006,'He_hot',0.),('clad',.048,'Nb1Zr',8.57),
        ('hot_He',spec['hot_gap_cm'],'He_hot',0.),
        ('hot_wall',.025,'Nb1Zr',8.57),
        ('vacuum_inner',spec['vacuum_gap_cm']/2,'void',0.),
        ('Mo_shield',.0025,'Mo',10.28),
        ('vacuum_outer',spec['vacuum_gap_cm']/2,'void',0.),
        ('cold_wall',.025,'Nb1Zr',8.57),
        ('cold_He',spec['cold_gap_cm'],'He_cold',0.),
        ('inner_can',.015,'Nb1Zr',8.57),
        ('inner_barrier',.008,'Al2O3',3.95),
        ('moderator',spec['moderator_thickness_cm'],'YH2',spec['moderator_density_g_cm3']),
        ('expansion_clearance',.015,'He_cold',0.),
        ('outer_barrier',.008,'Al2O3',3.95),
        ('outer_can',.015,'Nb1Zr',8.57)]
    if spec['outer_cold_channel']:
        definitions += [('outer_cold_He',spec['cold_gap_cm']*.65,'He_cold',0.),
                        ('outer_jacket',.015,'Nb1Zr',8.57)]
    radius = 0.; result=[]
    for name,thickness,material,rho in definitions:
        inner = radius; radius += thickness
        result.append(dict(name=name,r_inner_cm=inner,r_outer_cm=radius,
                           material=material,density_g_cm3=rho))
    return result


def metrics(spec):
    ll=layers(spec); pp=pins(spec); n=len(pp); h=spec['height_cm']
    r=spec['core_radius_cm']; rr=r+spec['radial_reflector_cm']
    fullh=h+2*spec['axial_reflector_cm']
    by={}
    for layer in ll:
        volume=n*math.pi*(layer['r_outer_cm']**2-layer['r_inner_cm']**2)*h
        by[layer['name']]=volume*layer['density_g_cm3']/1000
    # Upper preliminary mass includes reflector in axial pin passage locations;
    # the stochastic model calculation is the reported modeled mass if found.
    rho=1.85 if spec['reflector']=='Be' else 2.90
    be_volume=math.pi*rr**2*fullh-n*math.pi*ll[-1]['r_outer_cm']**2*h
    by['matrix_and_reflector_upper']=be_volume*rho/1000
    endnames={'hot_wall','Mo_shield','cold_wall','inner_can','inner_barrier',
              'outer_barrier','outer_can','outer_jacket'}
    by['axial_jacket_density_increment_upper']=sum(n*math.pi*(x['r_outer_cm']**2-x['r_inner_cm']**2)*
        2*spec['axial_reflector_cm']*max(0.,x['density_g_cm3']-rho)/1000 for x in ll if x['name'] in endnames)
    absorber_volume=spec['drum_count']*spec['absorber_sector_deg']/360*math.pi*(
        spec['drum_radius_cm']**2-(spec['drum_radius_cm']-spec['absorber_thickness_cm'])**2)*h
    by['absorber_density_increment_upper']=max(0.,2.52-rho)*absorber_volume/1000
    un=by['fuel'];hm=0.;fissile=0.
    mn=(1-spec['nitrogen15_atom_fraction'])*14.003074+spec['nitrogen15_atom_fraction']*15.000109
    for fraction,enr in [(19/n,spec['enrichment_inner']),((n-19)/n,spec['enrichment_outer'])]:
        mu=1/(enr/235.0439299+(1-enr)/238.0507884)
        zhm=un*fraction*mu/(mu+mn);hm+=zhm;fissile+=zhm*enr
    burn10=spec['power_W']*365.25*10/1e6/hm
    # 200 MeV/fission; a fuel-inventory screen, not a reactivity lifetime.
    consumption=spec['power_W']*10*365.25*86400/(200e6*1.602176634e-19)*235.0439/6.02214076e23/1000
    return dict(pin_count=n,UN_kg=un,heavy_metal_kg=hm,initial_U235_kg=fissile,
        estimated_fission_consumption_10yr_kg=consumption,
        burnup_10yr_MWd_kgHM=burn10,
        life_statement='Fuel inventory supports investigating >10 yr; no depletion or lifespan claim.',
        preliminary_solid_mass_upper_kg=sum(by.values()),
        helium_flow_kg_s=spec['power_W']/5193.16/(spec['outlet_K']-spec['inlet_K']),
        preliminary_mass_components_kg=by)


def catalog():
    # Each family samples sleeve/fuel ratio, leakage/size, reflector and drum
    # worth together. These are specified proposals, not numerical optimization.
    families=[
        ('double_cold_large',dict(rings=9,fuel_radius_cm=.64,moderator_thickness_cm=.22,
            outer_cold_channel=True,height_cm=62.,reflector='Be')),
        ('thin_sleeve_fuel_rich',dict(rings=10,fuel_radius_cm=.68,moderator_thickness_cm=.12,
            outer_cold_channel=False,height_cm=66.,reflector='Be')),
        ('thick_sleeve_double_cold',dict(rings=8,fuel_radius_cm=.62,moderator_thickness_cm=.30,
            outer_cold_channel=True,height_cm=62.,reflector='Be')),
        ('long_single_cold',dict(rings=9,fuel_radius_cm=.65,moderator_thickness_cm=.20,
            outer_cold_channel=False,height_cm=76.,reflector='Be')),
        ('beo_double_cold',dict(rings=9,fuel_radius_cm=.66,moderator_thickness_cm=.20,
            outer_cold_channel=True,height_cm=64.,reflector='BeO'))]
    variations=[(0.,0.,850.,1150.,4.8,2.0),
                (-6.,-.03,850.,1150.,5.1,2.3),
                (6.,.03,850.,1150.,5.1,2.3),
                (0.,.06,850.,1100.,5.3,2.5),
                (4.,0.,875.,1125.,5.3,2.5),
                (8.,-.02,900.,1100.,5.5,2.6)]
    result=[]
    for family,base in families:
        for j,(dh,dt,tin,tout,dr,ab) in enumerate(variations):
            spec=dict(power_W=400000.,pressure_Pa=3.0e6,moderator_x=2.0,
                fuel_density_g_cm3=13.1,moderator_density_g_cm3=4.20,
                enrichment_inner=.196,enrichment_outer=.1975,hot_gap_cm=.11,
                cold_gap_cm=.08,vacuum_gap_cm=.08,drum_count=12,
                drum_radius_cm=dr,absorber_thickness_cm=ab,
                absorber_sector_deg=175.,b10_atom_fraction=.95,
                drum_clearance_cm=.18,axial_reflector_cm=11.,
                inlet_K=tin,outlet_K=tout,**base)
            spec['nitrogen15_atom_fraction']=[.995,.00364,.95,.995,.995,.995][j]
            if j in (2,5):spec['rings']+=1
            spec['height_cm']+=dh;spec['moderator_thickness_cm']+=dt
            router=layers(spec)[-1]['r_outer_cm']
            spec['pitch_cm']=round(2*router+.08,4)
            spec['core_radius_cm']=spec['rings']*spec['pitch_cm']+router+.20
            spec['radial_reflector_cm']=2*dr+.80
            spec['drum_centres_radius_cm']=spec['core_radius_cm']+.18+dr
            result.append(dict(id=f'C{len(result)+1:02d}_{family}_v{j+1}',
                family=family,spec=spec,metrics=metrics(spec),
                rationale='Cold-first serial return flow; insulated hot bore; sleeve ratio, core size, N15 fraction and absorber worth varied.'))
    assert len(result)==30
    return result


def main():
    cc=catalog();(ROOT/'configurations.json').write_text(json.dumps(cc,indent=2)+'\n')
    columns=['id','family','pin_count','fuel_radius_cm','moderator_thickness_cm','height_cm',
             'pitch_cm','reflector','outer_cold_channel','nitrogen15_atom_fraction','inlet_K','outlet_K','initial_U235_kg',
             'burnup_10yr_MWd_kgHM','preliminary_solid_mass_upper_kg']
    with (ROOT/'configurations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader()
        for c in cc:
            merged=c|c['spec']|c['metrics'];w.writerow({k:merged[k] for k in columns})


if __name__=='__main__':main()
