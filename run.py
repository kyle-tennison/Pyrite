#!/usr/bin/env python3
from pyrite.solver import Solver
from pyrite.mesh import Mesher

solver = Solver()
mesher = Mesher()

elements, nodes = mesher.parse_nodes("nodes.csv")

for node in nodes:
    print(str(node))

solver.run(nodes, elements)