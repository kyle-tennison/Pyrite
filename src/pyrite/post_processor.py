"""
Post-processor to solve for planar stresses and display results

Kyle Tennison
March 24, 2024
"""

from pyrite.datatypes import Node, DOF
from pyrite.element import Element

from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
import math


class PostProcessor:

    @staticmethod
    def load_displacements_into_elements(
        nodal_displacements: np.ndarray, elements: list[Element]
    ) -> None:
        """Loads the nodal displacements into each element object

        Args:
            nodal_displacements: An array of nodal displacements
            elements: A list of Elements
        """

        for element in elements:
            n1_ux = nodal_displacements[2 * element.n1.index]
            n1_uy = nodal_displacements[2 * element.n1.index + 1]

            n2_ux = nodal_displacements[2 * element.n2.index]
            n2_uy = nodal_displacements[2 * element.n2.index + 1]

            n3_ux = nodal_displacements[2 * element.n3.index]
            n3_uy = nodal_displacements[2 * element.n3.index + 1]

            element.n1_u = (n1_ux, n1_uy)
            element.n2_u = (n2_ux, n2_uy)
            element.n3_u = (n3_ux, n3_uy)

    def compute_stress(self, elements: list[Element]) -> None:
        """Computes the stress from the solved nodal displacements and
        loads into element. This means load_displacements_into_elements
        must be called first.

        Args:
            elements: A list of Elements
        """

        for element in elements:

            nodal_displacements = np.array(
                [element.n1_u, element.n2_u, element.n3_u]
            ).reshape(3 * DOF, 1)

            stress = (
                element.stress_strain_matrix
                @ element.strain_displacement_matrix
                @ nodal_displacements
            ).flatten()

            sigma_x = stress[0]
            sigma_y = stress[1]

            element.stress = math.sqrt(sigma_x**2 + sigma_y**2)

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

        if max_value == min_value:
            return "#4C4C4C"

        norm_value = (value - min_value) / (max_value - min_value)

        blue = 0.0
        red = 0.0

        if norm_value > 0.5:
            red = 2 * (norm_value - 0.5)

        else:
            blue = 2 * norm_value

        # Combine RGB values into a hex string
        color = "#%02x%02x%02x" % (int(255 * red), 0, int(255 * blue))

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

            if element.stress > max:
                max = element.stress

            elif element.stress < min:
                min = element.stress

        return (min, max)

    def output_solved(
        self,
        nodal_forces: np.ndarray,
        nodal_displacements: np.ndarray,
        nodes: list[Node],
    ) -> None:
        """Outputs the solved system to output.csv

        Args:
            nodal_forces: The vector of solved nodal forces
            nodal_displacements: The vector of solved nodal displacements
            nodes: The list of Nodes
        """

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

    def show(self, elements: list[Element]) -> None:
        """Shows post-processed results with matplotlib

        Args:
            elements: The list of elements to plot
        """

        plt.style.use("seaborn-v0_8")

        self.compute_stress(elements)

        # Setup plot
        fig, axs = plt.subplots(2)
        fig.suptitle("Simulation Results")

        solved_plot = axs[0]
        initial_plot = axs[1]

        # Show initial plot
        triangles = np.empty((len(elements), 3, 2))

        for i, element in enumerate(elements):

            n1 = element.n1
            n2 = element.n2
            n3 = element.n3
            triangles[i, 0] = (n1.x, n1.y)
            triangles[i, 1] = (n2.x, n2.y)
            triangles[i, 2] = (n3.x, n3.y)

        for triangle in triangles:

            polygon = Polygon(
                triangle, closed=True, edgecolor="black", linewidth=0.2, alpha=0.7
            )

            polygon.set_facecolor("#4C4C4C")

            initial_plot.add_patch(polygon)

        initial_plot.set_title("Initial Model")

        # Show final plot
        triangles = np.empty((len(elements), 3, 2))
        triangle_colormap: list[str] = []

        min_stress, max_stress = self.find_max_stresses(elements)

        for i, element in enumerate(elements):

            n1 = element.n1
            n2 = element.n2
            n3 = element.n3

            n1_ux, n1_uy = element.n1_u
            n2_ux, n2_uy = element.n2_u
            n3_ux, n3_uy = element.n3_u

            triangles[i, 0] = (n1.x + n1_ux, n1.y + n1_uy)
            triangles[i, 1] = (n2.x + n2_ux, n2.y + n2_uy)
            triangles[i, 2] = (n3.x + n3_ux, n3.y + n3_uy)

            color = self.color_gradient(min_stress, max_stress, element.stress)
            triangle_colormap.append(color)

        for i, triangle in enumerate(triangles):

            polygon = Polygon(
                triangle, closed=True, edgecolor="black", linewidth=0.2, alpha=0.7
            )

            polygon.set_facecolor(triangle_colormap[i])

            solved_plot.add_patch(polygon)

        solved_plot.set_title("Solved Model")

        solved_plot.autoscale()
        initial_plot.autoscale()

        # Adjust axes to be equal to each other, and to fit each other
        if not (solved_plot.get_xlim() > initial_plot.get_xlim()):
            initial_plot.set_xlim(solved_plot.get_xlim())
        else:
            solved_plot.set_xlim(initial_plot.get_xlim())

        if not (solved_plot.get_ylim() > initial_plot.get_ylim()):
            initial_plot.set_ylim(solved_plot.get_ylim())
        else:
            solved_plot.set_ylim(initial_plot.get_ylim())

        fig.tight_layout(pad=2.0)
        plt.show()
