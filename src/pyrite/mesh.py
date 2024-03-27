"""
Meshing utility build around Gmsh. Handles mesh generation and boundary
condition application.

March 25, 2024
Kyle Tennison
"""

from pyrite.datatypes import Node, MshState, PartMetadata
from pyrite.element import Element
from pyrite.error import InputError

from matplotlib.patches import Polygon
from matplotlib import pyplot as plt
import json
import os
from typing import Callable
import numpy as np
import subprocess


def try_float(string: str):

    if string.strip() == "null":
        return None
    else:
        return float(string)


class Mesh:

    def __init__(self, vertices: np.ndarray):
        self.vertices = vertices


class BoundaryRule:

    SUPPORTED_TARGETS = ("ux", "uy", "fx", "fy")

    def __init__(self, name: str, targets: dict[str, float], region: dict[str, float]):

        self.name = name
        self.checks: list[Callable] = []
        self.targets = targets
        self._region = region
        self.part_metadata: PartMetadata

        self._parse_region(region)

    def _register_check(self, check_id: str, value: float):
        """Registers a check for the boundary rule

        Args:
            check_id: The check identifier, like 'x_target_min', for example.
            value: The corresponding value for the check_id
        """

        if check_id == "x_target_min":
            self.checks.append(lambda point: point.x >= value)
        elif check_id == "x_target_max":
            self.checks.append(lambda point: point.x <= value)
        elif check_id == "y_target_min":
            self.checks.append(lambda point: point.y >= value)
        elif check_id == "y_target_max":
            self.checks.append(lambda point: point.y <= value)
        else:
            raise RuntimeError(f"Unrecognized check identifier {check_id}")

    def _parse_region(self, region: dict[str, float]):
        """Parses a region into a list of checks

        Args:
            region: A dictionary mapping check ids to their corresponding values
        """

        for check_id, check_value in region.items():
            self._register_check(check_id, check_value)

    def maybe_apply_to_node(self, node: Node):
        """Applies targets to the node if it meets the criteria

        Args:
            node: The node to modify
        """

        if self.check(node):
            node.ux = self.targets["ux"]
            node.uy = self.targets["uy"]
            node.Fx = self.targets["fx"]
            node.Fy = self.targets["fy"]

    def check(self, point) -> bool:
        """Runs a point through all the checks for the boundary rule.

        Args:
            point: A value that gives x & y or r & θ values (depending on
                the coord system) when called with point.x, point.y, etc.
        """

        for check in self.checks:
            if not check(point):
                return False

        return True

    def __str__(self):
        return f"BoundaryRule({self.targets})"

    def __repr__(self) -> str:
        return str(self)


class Mesher:

    def __init__(self): ...

    def _generate_geo(
        self,
        outer_vertices: list[tuple],
        output_file: str,
        characteristic_length: float,
        characteristic_length_variance: float,
    ):
        """Generates a .geo file from a list of vertices.

        Args:
            outer_vertices: the ordered list of vertices to generate the
                outermost surface boundary.
            output_file: The output .geo filepath
            characteristic_length: The mean characteristic element length
            characteristic_length_variance: The ± variance allowed
                in the characteristic length.

        """

        ELEMENT_ORDER = 1
        ALGORITHM = 1  # delaunay

        with open("geom.geo", "w") as f:

            # define points
            f.write("// Define Points\n")
            for i, vertex in enumerate(outer_vertices):
                x, y, ux, uy, fx, fy = vertex

                f.write(f"Point({i}) = {{{x}, {y}, 0, 1.0}};\n")

            # connect points
            f.write("\n\n// Connect Points\n")
            for i in range(1, len(outer_vertices)):
                f.write(f"Line({i-1}) = {{{i-1}, {i}}};\n")
            f.write(
                f"Line({len(outer_vertices)-1}) = {{{len(outer_vertices)-1}, 0}};\n"
            )

            # define outer loop
            f.write("\n\n// Register outer loop\n")
            f.write("Line Loop(1) = {")
            for i in range(len(outer_vertices)):
                f.write(("," if i != 0 else "") + f"{i}")
            f.write("};\n")
            f.write("Plane Surface(1) = {1};")

            # define meshing settings
            f.write("\n\n// Define Mesh Settings\n")
            f.write(f"Mesh.ElementOrder = {ELEMENT_ORDER};\n")
            f.write(f"Mesh.Algorithm = {ALGORITHM};\n")
            f.write(
                f"Mesh.CharacteristicLengthMax = {characteristic_length + characteristic_length_variance};\n"
            )
            f.write(
                f"Mesh.CharacteristicLengthMin = {characteristic_length - characteristic_length_variance};\n"
            )
            f.write(f"Mesh 2;\n")

    def _compute_mesh(self, geo_file: str, output_file: str = "output.msh"):
        """Runs gmsh to compute mesh from .geo file

        Args:
            geo_file: The .geo file to target
        """

        try:
            subprocess.run(["gmsh", geo_file, "-2", "-o", output_file])
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to generate mesh: {str(e)}")

    def _parse_mesh(self, mesh_file: str) -> tuple[list[Node], list[Element]]:
        """Parses a mesh file, building nodes and elements.

        Args:
            mesh_file: The .msh file to target

        Returns:
            A list of node objects and a list of three-pair tuples that
            match elements to nodes.
        """

        nodes: dict[int, Node] = {}
        triangles = []
        state = MshState.LIMBO
        parsed_metadata = False
        skipped_elements = 0

        with open(mesh_file, "r") as f:

            while f.readable():

                line = f.readline()

                # detect state change
                if state is MshState.LIMBO:

                    parsed_metadata = False

                    if line == "":
                        break

                    if line.startswith("$Entities"):
                        state = MshState.ENTITIES
                    elif line.startswith("$Nodes"):
                        state = MshState.NODES
                    elif line.startswith("$Elements"):
                        state = MshState.ELEMENTS

                    continue

                # handle sections
                if line.startswith("$End"):
                    state = MshState.LIMBO
                    continue

                if state is MshState.ENTITIES:
                    continue

                if state is MshState.NODES:

                    # read metadata on first line
                    if not parsed_metadata:
                        num_nodes = int(line.strip().split(" ")[1])
                        parsed_metadata = True
                        continue

                    else:

                        metadata = [int(i) for i in line.strip().split(" ")]
                        entity_dim, entity_tag, parametric, num_nodes_local = metadata

                        # read the corresponding node tags
                        node_tags = []
                        for i in range(num_nodes_local):
                            tag = int(f.readline())
                            node_tags.append(tag)

                        assert len(node_tags) == num_nodes_local

                        # read proceeding dims and match to node tag
                        for i in range(num_nodes_local):
                            x, y, z = [
                                float(i) for i in f.readline().strip().split(" ")
                            ]
                            node = Node(
                                x=x,
                                y=y,
                                ux=None,
                                uy=None,
                                Fx=0,
                                Fy=0,
                                index=node_tags[i] - 1,
                            )
                            nodes[node_tags[i] - 1] = node

                if state is MshState.ELEMENTS:
                    if not parsed_metadata:
                        metadata = [int(i) for i in line.strip().split(" ")]
                        parsed_metadata = True
                    else:
                        metadata = [int(i) for i in line.strip().split(" ")]

                        entity_dim, entity_tag, element_type, num_elements = metadata

                        for i in range(num_elements):
                            data = [int(i) for i in f.readline().strip().split(" ")]

                            if entity_dim != 2:
                                skipped_elements += 1

                            else:
                                element_tag, n1, n2, n3 = data
                                triangles.append((n1, n2, n3))

        Element.material_elasticity = self.part_metadata.material_elasticity
        Element.poisson_ratio = self.part_metadata.poisson_ratio
        Element.part_thickness = self.part_metadata.part_thickness
        
        elements: list[Element] = []
        for n1_idx, n2_idx, n3_idx in triangles:
            n1 = nodes[n1_idx - 1]
            n2 = nodes[n2_idx - 1]
            n3 = nodes[n3_idx - 1]

            # NOTE: Triangles are given by note tags, which start at index 1
            # whereas our node indices start at index 0

            elements.append(Element(n1, n2, n3))

        print(f"info: skipped {skipped_elements} elements")
        print(f"info: registered {len(nodes)} nodes")
        print(f"info: registered {len(triangles)} triangles")

        return list(nodes.values()), elements

    def plot_elements(self, elements: list[Element]):
        """Plots the triangles formed from a list of elements using
        matplotlib.

        Args:
            elements: A list of Element objects to plot
        """

        triangles = np.empty((len(elements), 3, 2))

        for i, element in enumerate(elements):

            n1 = element.n1
            n2 = element.n2
            n3 = element.n3
            triangles[i, 0] = (n1.x, n1.y)
            triangles[i, 1] = (n2.x, n2.y)
            triangles[i, 2] = (n3.x, n3.y)

        fig, ax = plt.subplots()
        for triangle in triangles:

            polygon = Polygon(
                triangle, closed=True, edgecolor="black", linewidth=2, alpha=0.7
            )

            polygon.set_facecolor(tuple(np.random.random(3)))

            ax.add_patch(polygon)

        # Add labels and title
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_title("Triangles")

        # Show the plot
        plt.axis("equal")  # Equal aspect ratio
        plt.show()

    def _parse_csv(self, input_file: str) -> list[tuple]:
        """Parses input CSV that contains vertices

        Args:
            input_file: The filepath of the CSV file to reference

        Returns:
            A list of tuples that contain vertex coordinates and metadata
        """

        vertices = []

        with open(input_file, "r") as f:

            headers = [i.strip() for i in f.readline().strip().split(",")]

            for line in [i.strip().split(",") for i in f.readlines()]:

                x = float(line[headers.index("x")])
                y = float(line[headers.index("y")])
                ux = try_float(line[headers.index("ux")])
                uy = try_float(line[headers.index("uy")])
                fx = try_float(line[headers.index("fx")])
                fy = try_float(line[headers.index("fy")])

                vertices.append((x, y, ux, uy, fx, fy))

        return vertices
    

    def _load_metadata(self, input_file: str) -> PartMetadata:
        """Loads the metadata from the input file"""

        with open(input_file, "r") as f:
            try:
                input_file_data = json.load(f)
            except json.JSONDecodeError as e:
                raise InputError(f"Error in input file json: {str(e)}")

        try:
            metadata = input_file_data["metadata"]

            material_elasticity = metadata["material_elasticity"]
            poisson_ratio = metadata["poisson_ratio"]
            part_thickness = metadata["part_thickness"]

            return PartMetadata(
                material_elasticity, poisson_ratio, part_thickness
            )
        except KeyError as e:
            raise InputError(f"Input file missing field {str(e)}")
        

    def _apply_boundary_conditions(self, input_file: str, nodes: list[Node]):
        """Applies boundary conditions from the input file onto a list of nodes
        
        Args:
            input_file: The path to the input file
            nodes: A list of nodes from the mesh to apply the boundary conditions   
                to.
        """

        input_file_data = {}
        boundary_rules: list[BoundaryRule] = []

        with open(input_file, "r") as f:
            input_file_data = json.load(f)

        try:
            for check_name, check_data in input_file_data["boundary_conditions"].items():

                region = check_data["region"]
                targets = check_data["targets"]

                rule = BoundaryRule(check_name, targets, region)
                boundary_rules.append(rule)
        except KeyError as e:
            raise InputError(f"Input file missing field {str(e)}")

        for node in nodes:
            for rule in boundary_rules:
                rule.maybe_apply_to_node(node)

    def mesh(
        self,
        input_csv: str,
        input_file: str,
        characteristic_length: float,
        characteristic_length_variance: float,
    ) -> tuple[list[Node], list[Element]]:
        """Creates a mesh from a list of 2D vertices

        Args:
            input_csv: The filepath to the input CSV that contains a list of
                vertices
            characteristic_length: The mean characteristic element length.
            characteristic_length_variance: The ± variance allowed
                in the characteristic length.

        Returns:
            A list of Nodes and Elements
        """

        if not os.path.exists(input_file):
            raise InputError(f"Could not find input file at {os.path.abspath(input_file)}")
        
        if not os.path.exists(input_csv):
            raise InputError(f"Could not find vertices csv at {os.path.abspath(input_csv)}")

        vertices = self._parse_csv(input_csv)

        geo_filename = "geom.geo"
        mesh_filename = "geom.msh"

        self.part_metadata = self._load_metadata(input_file)

        self._generate_geo(
            vertices,
            geo_filename,
            characteristic_length,
            characteristic_length_variance,
        )
        self._compute_mesh(geo_filename, mesh_filename)
        nodes, elements = self._parse_mesh(mesh_filename)

        self._apply_boundary_conditions(input_file, nodes)

        # cleanup files
        os.remove(geo_filename)
        os.remove(mesh_filename)
        return nodes, elements
