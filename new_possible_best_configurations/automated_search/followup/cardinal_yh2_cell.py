"""2-D Cardinal FEM of the actual hexagonal YH2 pin cell and seven channels."""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import subprocess

import meshio
import numpy as np
import triangle


HERE = Path(__file__).resolve().parent
CARDINAL = Path("/home/aleri/projects/cardinal/cardinal-opt")
MESH = HERE / "thermal/cell_yh2.msh"
PITCH_CM = 2.16
R_MAIN_CM = .585+.048+.084+.010
R_AUX_CM = .10
R_AUX_POS_CM = .910*(PITCH_CM/2.270)


def make_mesh(area_cm2=0.0004):
    MESH.parent.mkdir(parents=True, exist_ok=True)
    pts=[];segments=[];tags=[]
    def ring(xy,tag):
        first=len(pts)
        pts.extend(xy)
        for j in range(len(xy)):
            segments.append((first+j,first+(j+1)%len(xy)))
            tags.append(tag)
    # Hexagonal Voronoi cell with radial pitch PITCH_CM.
    circumradius=PITCH_CM/math.sqrt(3)
    corners=np.array([[circumradius*math.cos(math.radians(30+60*j)),
                       circumradius*math.sin(math.radians(30+60*j))]
                       for j in range(6)])
    outer=[]
    for j in range(6):
        a=corners[j];b=corners[(j+1)%6]
        for q in range(12):outer.append(a+(b-a)*q/12)
    ring(outer,4)
    circle=lambda x,y,r,n:[[x+r*math.cos(2*math.pi*j/n),
                            y+r*math.sin(2*math.pi*j/n)] for j in range(n)]
    ring(circle(0,0,R_MAIN_CM,96),2)
    holes=[[0.,0.]]
    for j in range(6):
        x=R_AUX_POS_CM*math.cos(j*math.pi/3)
        y=R_AUX_POS_CM*math.sin(j*math.pi/3)
        ring(circle(x,y,R_AUX_CM,28),3)
        holes.append([x,y])
    domain=dict(vertices=np.asarray(pts),segments=np.asarray(segments),
                segment_markers=np.asarray(tags).reshape(-1,1),
                holes=np.asarray(holes))
    tri=triangle.triangulate(domain,f"pq25a{area_cm2}")
    if "triangles" not in tri:raise RuntimeError("Triangle mesher failed")
    vertices=np.c_[tri["vertices"]*.01,np.zeros(len(tri["vertices"]))]
    lines=np.asarray(tri["segments"],dtype=int)
    markers=np.asarray(tri["segment_markers"],dtype=int).ravel()
    mesh=meshio.Mesh(vertices,
        cells=[("triangle",np.asarray(tri["triangles"],dtype=int)),
               ("line",lines)],
        cell_data={"gmsh:physical":[np.ones(len(tri["triangles"]),dtype=int),markers],
                   "gmsh:geometrical":[np.ones(len(tri["triangles"]),dtype=int),markers]},
        field_data={"YH2":np.array([1,2]),"main_He":np.array([2,1]),
                    "aux_He":np.array([3,1]),"outer_hex":np.array([4,1])})
    meshio.write(MESH,mesh,file_format="gmsh22",binary=False)
    va=tri["vertices"][tri["triangles"][:,1]]-tri["vertices"][tri["triangles"][:,0]]
    vb=tri["vertices"][tri["triangles"][:,2]]-tri["vertices"][tri["triangles"][:,0]]
    areas=.5*np.abs(va[:,0]*vb[:,1]-va[:,1]*vb[:,0])
    info=dict(triangles=int(len(areas)),nodes=int(len(vertices)),
              area_cm2=float(areas.sum()),
              exact_area_cm2=float(math.sqrt(3)/2*PITCH_CM**2-
                                   math.pi*R_MAIN_CM**2-6*math.pi*R_AUX_CM**2),
              boundary_segment_counts={str(i):int(sum(markers==i)) for i in (2,3,4)})
    (MESH.parent/"mesh_info.json").write_text(json.dumps(info,indent=2)+"\n")
    return info


def solve(label, *, t_main_K, t_aux_K, qvol_W_m3,
          h_main_W_m2K=700., h_aux_W_m2K=450., k_YH2_W_mK=10.):
    if not MESH.exists():make_mesh()
    d=HERE/"thermal/cardinal"/label
    d.mkdir(parents=True,exist_ok=True)
    input_text=f'''[Mesh]
  type = FileMesh
  file = {MESH}
[]
[Variables]
  [T]
    initial_condition = {max(t_main_K,t_aux_K)}
  []
[]
[Functions]
  [heat]
    type = ConstantFunction
    value = {qvol_W_m3}
  []
  [main_bulk]
    type = ConstantFunction
    value = {t_main_K}
  []
  [aux_bulk]
    type = ConstantFunction
    value = {t_aux_K}
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = heat
  []
[]
[BCs]
  [main]
    type = ConvectiveFluxFunction
    variable = T
    boundary = '2'
    coefficient = {h_main_W_m2K}
    T_infinity = main_bulk
  []
  [aux]
    type = ConvectiveFluxFunction
    variable = T
    boundary = '3'
    coefficient = {h_aux_W_m2K}
    T_infinity = aux_bulk
  []
[]
[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = {k_YH2_W_mK}
  []
  [main_bulk]
    type = GenericConstantMaterial
    prop_names = main_bulk
    prop_values = {t_main_K}
  []
  [aux_bulk]
    type = GenericConstantMaterial
    prop_names = aux_bulk
    prop_values = {t_aux_K}
  []
[]
[Postprocessors]
  [T_YH2_max_K]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
  [T_YH2_min_K]
    type = NodalExtremeValue
    variable = T
    value_type = min
  []
  [Q_main_W_per_m]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = main_bulk
    htc = {h_main_W_m2K}
    boundary = '2'
  []
  [Q_aux_W_per_m]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = aux_bulk
    htc = {h_aux_W_m2K}
    boundary = '3'
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
[]
[Outputs]
  csv = true
  exodus = true
[]
'''
    (d/"case.i").write_text(input_text)
    env=os.environ.copy()
    env.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (d/"run.log").open("w") as log:
        cp=subprocess.run([str(CARDINAL),"-i","case.i","--n-threads=1"],
                          cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
    if cp.returncode:raise RuntimeError(f"Cardinal failed: {d/'run.log'}")
    with (d/"case_out.csv").open() as file:
        result={key:float(value) for key,value in list(csv.DictReader(file))[-1].items()}
    result.update(label=label,t_main_K=t_main_K,t_aux_K=t_aux_K,
                  qvol_W_m3=qvol_W_m3,h_main_W_m2K=h_main_W_m2K,
                  h_aux_W_m2K=h_aux_W_m2K,k_YH2_W_mK=k_YH2_W_mK)
    (d/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":
    print(json.dumps(make_mesh(),indent=2))
    print(json.dumps(solve("smoke",t_main_K=1100,t_aux_K=800,
                           qvol_W_m3=1e6),indent=2))
