from pyrite.datatypes import Node, Element
from pyrite.mesh import try_float

import numpy as np
import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


class PostProcessor:

    def compute_strain(self, elements: list[Element], nodal_displacements: np.ndarray):
        """Computes the strain from the solved nodal displacements and
        loads into element."""

    def output_solved(
        self,
        nodal_forces: np.ndarray,
        nodal_displacements: np.ndarray,
        nodes: list[Node],
    ):
        """Outputs the solved system to output.csv"""

        with open("output.csv", "w") as f:
            f.write("x, y, ux, uy, fx, fy\n")

            i = 0
            for node in nodes:

                fx = nodal_forces[i]
                fy = nodal_forces[i + 1]

                ux = nodal_displacements[i]
                uy = nodal_displacements[i + 1]

                f.write(f"{node.x},{node.y},{ux},{uy},{fx},{fy}\n")

                i += 2

    def show(self, input_file: str):
        """Shows post-processed results"""

        fig, axs = plt.subplots(2)
        fig.suptitle("Simulation Results")

        solved_plot = axs[0]
        initial_plot = axs[1]

        begin_points = []
        with open(input_file, "r") as f:
            header_skipped = False

            for line in f.readlines():
                if not header_skipped:
                    header_skipped = True
                    continue

                items = line.split(",")

                x = float(items[0])
                y = float(items[1])
                ux = try_float(items[2])
                uy = try_float(items[3])

                if ux is None:
                    ux = 0
                if uy is None:
                    uy = 0

                begin_points.append((x + ux, y + uy))

        begin_points = np.array(begin_points)
        tri_begin = Delaunay(begin_points)

        initial_plot.set_title("Initial Model:")
        initial_plot.triplot(
            begin_points[:, 0], begin_points[:, 1], tri_begin.simplices
        )
        # initial_plot.plot(begin_points[:,0], begin_points[:,1], 'o')

        solved_points = []
        with open("output.csv", "r") as f:
            header_skipped = False

            for line in f.readlines():
                if not header_skipped:
                    header_skipped = True
                    continue

                items = line.split(",")

                x = float(items[0])
                y = float(items[1])
                ux = float(items[2])
                uy = float(items[3])

                solved_points.append((x + ux, y + uy))
        solved_points = np.array(solved_points)
        tri_solved = Delaunay(solved_points)
        print("solved points:\n", solved_points)

        solved_plot.set_title("Deformed Results:")
        solved_plot.triplot(
            solved_points[:, 0], solved_points[:, 1], tri_begin.simplices
        )
        # solved_plot.plot(solved_points[:,0], solved_points[:,1])

        if not (solved_plot.get_xlim() > initial_plot.get_xlim()):
            print("1")
            initial_plot.set_xlim(solved_plot.get_xlim())
        else:
            print("2")
            solved_plot.set_xlim(initial_plot.get_xlim())

        if not (solved_plot.get_ylim() > initial_plot.get_ylim()):
            print("3")
            initial_plot.set_ylim(solved_plot.get_ylim())
        else:
            print("4")
            solved_plot.set_ylim(initial_plot.get_ylim())

        fig.tight_layout(pad=2.0)
        plt.show()
