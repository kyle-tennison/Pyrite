#!/usr/bin/env python3
from pyrite.solver import Solver
from pyrite.mesh import Mesher
from pyrite.datatypes import Node
from pyrite.element import Element

solver = Solver()
mesher = Mesher()


nodes, elements = mesher.mesh(
    input_csv="vertices.csv",
    boundary_file="boundary.json",
    characteristic_length=0.5,
    characteristic_length_variance=0.06,
)


mesher.plot_elements(elements)

# nodes = [
#     Node(3, 0, None, 0, 0, None, 0),
#     Node(3, 2, None, None, 0, -1000, 1),
#     Node(0, 2, 0, 0, None, None, 2),
#     Node(0, 0, 0, 0, None, None, 3)
# ]

# elements = [
#     Element(nodes[0], nodes[1], nodes[3]),
#     Element(nodes[2], nodes[1], nodes[3]),
#     # Element(nodes[2], nodes[1], nodes[3]),
# ]

# for row in elements[1].global_stiffness_matrix / 1E7:
#     print([f"{float(i):.3}" for i in row])

# mesher.plot_elements(elements)


solver.run(nodes, elements)
