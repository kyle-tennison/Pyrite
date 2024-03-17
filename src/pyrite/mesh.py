from pyrite.datatypes import Node, ElementLight

from scipy.spatial import Delaunay
import numpy as np
from itertools import combinations


def try_float(string: str):

    if string.strip() == "null":
        return None
    else:
        return float(string)

class Mesh:

    def __init__(self, vertices: np.ndarray):
        self.vertices = vertices

class Mesher:

    def __init__(self): ...

    def generate_geo(self, outer_vertices: list[tuple]):
        """Generates a .geo file from a list of vertices.
        
        Args:
            outer_vertices: the ordered list of vertices to generate the 
                outermost surface boundary.
        
        """

        ELEMENT_ORDER = 1
        ALGORITHM = 1 # delaunay
        CHARACTERISTIC_LENGTH_MAX = 1 # min element size
        CHARACTERISTIC_LENGTH_MIN = 5 # max element size

        with open("geom.geo", "w") as f:


            # define points
            f.write("// Define Points\n")
            for i, vertex in enumerate(outer_vertices):
                x = vertex[0]
                y = vertex[1]

                f.write(f"Point({i}) = {{{x}, {y}, 0, 1.0}};\n")

            # connect points
            f.write("\n\n// Connect Points\n")
            for i in range(1, len(outer_vertices)):
                f.write(f"Line({i-1}) = {{{i-1}, {i}}};\n")
            f.write(f"Line({len(outer_vertices)-1}) = {{{len(outer_vertices)-1}, 0}};\n")

            # define outer loop
            f.write("\n\n// Register outer loop\n")
            f.write("Line Loop(1) = {")
            for i in range(len(outer_vertices)):
                f.write(("," if i!=0 else "") + f"{i}")
            f.write("};\n")
            f.write("Plane Surface(1) = {1};")

            # define meshing settings
            f.write("\n\n// Define Mesh Settings\n")
            f.write(f"Mesh.ElementOrder = {ELEMENT_ORDER};\n")
            f.write(f"Mesh.Algorithm = {ALGORITHM};\n")
            f.write(f"Mesh.CharacteristicLengthMax = {CHARACTERISTIC_LENGTH_MAX};\n")
            f.write(f"Mesh.CharacteristicLengthMin = {CHARACTERISTIC_LENGTH_MIN};\n")
            f.write(f"Mesh 2;\n")

            

    def parse_csv(self, input_file: str):
        """Parses input CSV that contains vertices"""

        vertices = []

        with open(input_file, 'r') as f:

            headers = [i.strip() for i in f.readline().split(",")]

            for line in [i.split(",") for i in f.readlines()]:

                x = float(line[headers.index("x")])
                y = float(line[headers.index("y")])
                ux = try_float(line[headers.index("ux")])
                uy = try_float(line[headers.index("uy")])
                fx = try_float(line[headers.index("fx")])
                fy = try_float(line[headers.index("fy")])

                vertices.append(
                    (x,y,ux,uy,fx,fy)
                )

        self.generate_geo(vertices)

