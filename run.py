#!/usr/bin/env python3
from pyrite.solver import Solver
from pyrite.mesh import Mesher

solver = Solver()
mesher = Mesher()

nodes, elements = mesher.mesh(
    input_csv="vertices.csv",
    input_file="boundary.json",
    characteristic_length=0.5,
    characteristic_length_variance=0.06,
)

solver.run(nodes, elements)
