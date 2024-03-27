"""
Entry point for cli interface.

Kyle Tennison
March 26, 2024
"""

import argparse
from pyrite.solver import Solver
from pyrite.mesh import Mesher


def entry():
    """Entry point of the cli"""

    # Collect cli arguments
    parser = argparse.ArgumentParser(
        prog="pyrite",
        usage="pyrite vertices.csv boundary.json",
        description="A 2D Linear Elastic Finite Element Mechanical Simulator",
    )

    parser.add_argument(
        "vertices", help="A CSV file that contains the vertices of the 2D geometry"
    )
    parser.add_argument(
        "boundary",
        help="A JSON file that defines the boundary conditions for the study",
    )
    parser.add_argument(
        "-cl",
        "--characteristic_length",
        help="Characteristic length of the mesh",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "-cv",
        "--characteristic_length_variance",
        help="Characteristic length variance of the mesh",
        default=0.06,
        type=float,
    )

    args = parser.parse_args()

    # Run pyrite
    solver = Solver()
    mesher = Mesher()

    nodes, elements = mesher.mesh(
        vertex_csv=args.vertices,
        input_file=args.boundary,
        characteristic_length=args.characteristic_length,
        characteristic_length_variance=args.characteristic_length_variance,
    )

    solver.run(nodes, elements)


if __name__ == "__main__":
    entry()
