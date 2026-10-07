[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = /home/aleri/projects/cardinal/test/tests/neutronics/meshes/sphere.e
  []
[]

[Problem]
  type = OpenMCCellAverageProblem
  batches = 130 # this number needs to match the value in settings.xml for the purpose of the test
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Postprocessors]
  [k]
    type = KEigenvalue
  []
[]

[Outputs]
  csv = true
[]
