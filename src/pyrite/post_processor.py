import math

from matplotlib.colors import rgb2hex
from pyrite.datatypes import Node
from pyrite.mesh import try_float
from pyrite.element import Element

import numpy as np
import csv
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt


class PostProcessor:

    @staticmethod
    def load_displacements_into_elements(
        nodal_displacements: np.ndarray, elements: list[Element]
    ):
        """Loads the nodal displacements into each element object"""

        for element in elements:
            n1_ux = nodal_displacements[2 * element.n1.index]
            n1_uy = nodal_displacements[2 * element.n1.index + 1]

            n2_ux = nodal_displacements[2 * element.n2.index]
            n2_uy = nodal_displacements[2 * element.n2.index + 1]

            element.n1_u = (n1_ux, n1_uy)
            element.n2_u = (n2_ux, n2_uy)

    def compute_strain(self, elements: list[Element]):
        """Computes the strain from the solved nodal displacements and
        loads into element."""

        for element in elements:

            n1_ux, n1_uy = element.n1_u
            n2_ux, n2_uy = element.n2_u

            dx = n1_ux + n2_ux
            dy = n1_uy + n2_uy

            d = math.sqrt(dx**2 + dy**2)

            element.normal_strain = d / element.length

    @staticmethod
    def color_gradient(min_value: float, max_value: float, value: float) -> str:
        """Computes a color for a given value within an interval. Low will be
        blue, high will be red.

        Args:
            value: The value to convert to a color (float).
            min_value: The minimum value in the range (float).
            max_value: The maximum value in the range (float).

        Returns:
            A string representing the RGB color in hexadecimal format.
        """
        norm_value = (value - min_value) / (max_value - min_value)

        blue = 0.0
        red = 0.0

        if norm_value > 0.5:
            blue = 2 * (norm_value - 0.5)

        else:
            red = 2 * norm_value

        # Combine RGB values into a hex string
        color = "#%02x%02x%02x" % (int(255 * red), 0, int(255 * blue))

        print(
            f"max: {max_value}, min: {min_value}, value: {value}, red: {red}, blue:{blue} norm: {norm_value}"
        )

        return color

    @staticmethod
    def find_max_stresses(elements: list[Element]) -> tuple[float, float]:
        """Searches for the highest and lowest stresses in a list
        of elements. Note, the elements must have stresses computed
        before calling this function.

        Args:
            elements: A list of Element objects

        Returns:
            A tuple of (min, max) floats representing the min & max
                stresses
        """

        min = 0.0
        max = 0.0

        for element in elements:

            if element.normal_stress > max:
                max = element.normal_stress

            elif element.normal_strain < min:
                min = element.normal_stress

        return (min, max)

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

    def show(self, input_file: str, elements: list[Element]):
        """Shows post-processed results"""

        # Setup plot
        fig, axs = plt.subplots(2)
        fig.suptitle("Simulation Results")

        solved_plot = axs[0]
        initial_plot = axs[1]

        # Read input csv and plot initial results
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

        min_stress, max_stress = self.find_max_stresses(elements)

        # Draw deformed results with stress colors
        for element in elements:

            x = np.empty(2)
            y = np.empty(2)

            x[0] = element.n1.x + element.n1_u[0]
            x[1] = element.n2.x + element.n2_u[0]

            y[0] = element.n1.y + element.n1_u[1]
            y[1] = element.n2.y + element.n2_u[1]

            color = self.color_gradient(min_stress, max_stress, element.normal_stress)

            print("color is:", color)

            solved_plot.plot(x, y, color=color)

        # Adjust axes to be equal to each other, and to fit each other
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
