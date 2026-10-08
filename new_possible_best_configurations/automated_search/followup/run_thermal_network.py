"""Coupled He-flow energy balances and Cardinal YH2 cell solves at 400 kWth.

The six auxiliary holes and fuel annulus are modeled as parallel 800 K He
channels. Flow split is uncertain because manifolds/orifices are unspecified,
so several explicit fractions are evaluated. Outer hex edges are adiabatic.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from cardinal_yh2_cell import HERE, MESH, make_mesh, solve


POWER_W=400_000.0
TIN_K=800.0
TOUT_MIXED_K=1150.0
CP_HE=5193.16
M_TOTAL=POWER_W/(CP_HE*(TOUT_MIXED_K-TIN_K))
N=217
S=6
PITCH_CM=2.16
H_M=.47
DZ=H_M/S
MAIN_INNER_M=(.585+.048)*.01
MAIN_OUTER_M=(.585+.048+.084)*.01
LINER_M=.010*.01
AUX_R_M=.10*.01
MOD_AREA_M2=(math.sqrt(3)/2*PITCH_CM**2-
             math.pi*(.585+.048+.084+.010)**2-6*math.pi*.10**2)*1e-4


def helium_k(T_K):
    pressure=2.5e6
    return 2.682e-3*(1+1.123e-8*pressure)*T_K**(.71*(1-2e-9*pressure))


def film_coefficients(mult=1.0):
    k=helium_k(1000.0)
    dh_main=2*(MAIN_OUTER_M-MAIN_INNER_M)
    h_main=5.385*k/dh_main*mult
    h_main_eff=1/(1/h_main+LINER_M/6.0)
    h_aux=3.66*k/(2*AUX_R_M)*mult
    return h_main_eff,h_aux


def calibration(label,h_main,h_aux,kmod):
    same=solve(label+"_source",t_main_K=1000,t_aux_K=1000,
               qvol_W_m3=1e6,h_main_W_m2K=h_main,
               h_aux_W_m2K=h_aux,k_YH2_W_mK=kmod)
    diff=solve(label+"_delta",t_main_K=1100,t_aux_K=1000,
               qvol_W_m3=0,h_main_W_m2K=h_main,
               h_aux_W_m2K=h_aux,k_YH2_W_mK=kmod)
    qs=1e6*MOD_AREA_M2
    fmain=same["Q_main_W_per_m"]/qs
    g=-diff["Q_main_W_per_m"]/100
    if not (0<fmain<1 and g>0):
        raise RuntimeError(f"Invalid thermal calibration: fmain={fmain}, G={g}")
    energy_err=(same["Q_main_W_per_m"]+same["Q_aux_W_per_m"]-qs)/qs
    if abs(energy_err)>.01:
        raise RuntimeError(f"Cardinal energy balance: {energy_err}")
    return dict(fraction_moderator_heat_to_main=fmain,
                conductance_main_to_aux_W_mK=g,
                calibration_energy_error=energy_err,
                source_exodus=str(HERE/"thermal/cardinal"/(label+"_source")/"case_out.e"),
                delta_exodus=str(HERE/"thermal/cardinal"/(label+"_delta")/"case_out.e"))


def nodal_final_temperature(exodus):
    with Dataset(exodus) as data:
        # Cardinal writes the initial state and then the converged steady state.
        return np.asarray(data.variables["vals_nod_var1"][-1, :])


def flow_case(fuel,mod,*,aux_fraction,other_to_mod_frac=0.,
              film_mult=1.,kmod=10.):
    hmain,haux=film_coefficients(film_mult)
    label=f"a{aux_fraction:.2f}_o{other_to_mod_frac:.2f}_f{film_mult:.2f}_k{kmod:.0f}"
    cal=calibration(label,hmain,haux,kmod)
    other=POWER_W-float(fuel.sum())-float(mod.sum())
    if other<0:raise RuntimeError("Fuel + moderator heating exceeds power")
    weights=fuel/float(fuel.sum())
    qmod=mod+other*other_to_mod_frac*weights
    qmain_base=fuel+other*(1-other_to_mod_frac)*weights
    cm=M_TOTAL*(1-aux_fraction)/N*CP_HE
    ca=M_TOTAL*aux_fraction/N*CP_HE
    if cm<=0 or ca<=0:raise ValueError("Both He streams must have positive flow")
    main=np.full(N,TIN_K)
    aux=np.full(N,TIN_K)
    main_mid=np.zeros((N,S));aux_mid=np.zeros((N,S))
    mod_source=np.zeros((N,S));qmain_net=np.zeros((N,S))
    conductance=cal["conductance_main_to_aux_W_mK"]*DZ
    f=cal["fraction_moderator_heat_to_main"]
    for s in range(S):
        q=qmod[:,s];qb=qmain_base[:,s]
        numerator=f*q-conductance*(main-aux+qb/(2*cm)-q/(2*ca))
        denominator=1+conductance*(1/(2*cm)+1/(2*ca))
        q_to_main=numerator/denominator
        m_out=main+(qb+q_to_main)/cm
        a_out=aux+(q-q_to_main)/ca
        main_mid[:,s]=(main+m_out)/2
        aux_mid[:,s]=(aux+a_out)/2
        mod_source[:,s]=q
        qmain_net[:,s]=q_to_main
        main,aux=m_out,a_out
    mixed=float(((1-aux_fraction)*main+aux_fraction*aux).mean())
    if abs(mixed-TOUT_MIXED_K)>.05:
        raise RuntimeError(f"Flow balance failed: mixed={mixed}")
    # The FEM equation is linear for each fixed film/conductivity case. Combine
    # its two Cardinal basis fields at every node and every pin/slice state,
    # then independently re-run Cardinal at the exact hottest state.
    source=nodal_final_temperature(cal["source_exodus"])
    delta=nodal_final_temperature(cal["delta_exodus"])
    source_profile=(source-1000.)/1e6
    delta_profile=(delta-1000.)/100.
    tv=aux_mid.reshape(-1,1)+(main_mid-aux_mid).reshape(-1,1)*delta_profile
    tv+=((qmod/(MOD_AREA_M2*DZ)).reshape(-1,1))*source_profile
    all_max=tv.max(axis=1)
    order=np.argsort(all_max)[::-1][:3]
    candidates=[]
    for rank,flat in enumerate(order):
        pin,s=divmod(int(flat),S)
        r=solve(label+f"_pin{pin:03d}_s{s}",
                t_main_K=float(main_mid[pin,s]),
                t_aux_K=float(aux_mid[pin,s]),
                qvol_W_m3=float(qmod[pin,s]/(MOD_AREA_M2*DZ)),
                h_main_W_m2K=hmain,h_aux_W_m2K=haux,
                k_YH2_W_mK=kmod)
        candidates.append(dict(pin=pin,slice=s,moderator_max_K=r["T_YH2_max_K"],
                               superposition_max_K=float(all_max[flat]),
                               moderator_min_K=r["T_YH2_min_K"],
                               main_bulk_mid_K=float(main_mid[pin,s]),
                               aux_bulk_mid_K=float(aux_mid[pin,s]),
                               moderator_source_W=float(qmod[pin,s]),
                               cardinal_result=str(HERE/"thermal/cardinal"/
                                                   (label+f"_pin{pin:03d}_s{s}")/"result.json")))
    peak=max(candidates,key=lambda x:x["moderator_max_K"])
    if abs(peak["moderator_max_K"]-peak["superposition_max_K"])>.05:
        raise RuntimeError("Cardinal nodal superposition disagrees with direct solve")
    return dict(label=label,aux_flow_fraction=aux_fraction,
                other_heat_to_moderator_fraction=other_to_mod_frac,
                k_YH2_W_mK=kmod,film_multiplier=film_mult,
                h_main_effective_W_m2K=hmain,h_aux_W_m2K=haux,
                calibration=cal,mixed_outlet_K=mixed,
                mean_main_outlet_K=float(main.mean()),
                max_main_outlet_K=float(main.max()),
                mean_aux_outlet_K=float(aux.mean()),
                max_aux_outlet_K=float(aux.max()),
                max_moderator_K=peak["moderator_max_K"],
                max_moderator_location=peak,
                evaluated_pin_slices=candidates,
                model_limitations=["Prescribed 400 kW and parallel 800 K He feeds",
                    "Auxiliary flow fraction is a sensitivity parameter; manifolds absent",
                    "Heat not in fuel/YH2 allocated by fuel map; see scenario fraction",
                    "Outer hex adiabatic; Be conduction not coupled",
                    "Constant YH2 conductivity, He cp and representative film coefficients",
                    "No pressure-drop, orifice or structural-can/contact solution"])


def main():
    make_mesh()
    data=np.load(HERE/"heating_maps_W.npz")
    fuel=data["fuel_W"];mod=data["moderator_W"]
    cases=[]
    specs=[(.20,0.,1.,10.),(.35,0.,1.,10.),(.55,0.,1.,10.),
           (.35,.50,1.,10.),(.55,.50,1.,10.),
           (.35,0.,.5,10.),(.55,0.,.5,10.),
           (.35,0.,1.,5.),(.55,0.,1.,5.)]
    for aux,other,film,kmod in specs:
        result=flow_case(fuel,mod,aux_fraction=aux,
                         other_to_mod_frac=other,film_mult=film,kmod=kmod)
        cases.append(result)
        (HERE/"thermal_results.json").write_text(json.dumps(dict(
            status="running",power_W=POWER_W,inlet_K=TIN_K,
            target_mixed_outlet_K=TOUT_MIXED_K,cases=cases),indent=2)+"\n")
        print(result["label"],"max YH2",result["max_moderator_K"],
              "main outlet",result["max_main_outlet_K"],flush=True)
    (HERE/"thermal_results.json").write_text(json.dumps(dict(
        status="complete",power_W=POWER_W,inlet_K=TIN_K,
        target_mixed_outlet_K=TOUT_MIXED_K,cases=cases),indent=2)+"\n")


if __name__=="__main__":main()
