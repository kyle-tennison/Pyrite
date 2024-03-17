#!/usr/bin/env python3
from pyrite.solver import Solver
from pyrite.mesh import Mesher

solver = Solver()
mesher = Mesher()


nodes, elements = mesher.mesh(
    input_csv="vertices.csv",
    characteristic_length=0.5,
    characteristic_length_variance=0.1,
)

mesher.plot_elements(elements)

# elements, nodes = mesher.parse_nodes("nodes.csv")

# for node in nodes:
#     print(str(node))

# solver.run(nodes, elements)
