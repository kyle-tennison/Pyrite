from pyrite.datatypes import Node, ElementLight

from scipy.spatial import Delaunay
import numpy as np
from itertools import combinations


def try_float(string: str):

    if string.strip() == "null":
        return None
    else:
        return float(string)


class Mesher:

    def __init__(self): ...

    @staticmethod
    def _pair_equals(p1: tuple, p2: tuple) -> bool:
        """Returns True if the paris contain the same items,
        regardless of their order; False otherwise.
        """

        if p1[0] == p2[0] and p1[1] == p2[1]:
            return True

        elif p1[0] == p2[1] and p1[1] == p2[0]:
            return True

        else:
            return False

    def parse_nodes(self, nodes_csv: str) -> tuple[list[tuple[int, int]], list[Node]]:
        """Parses the nodes in a csv_file. Creates a list of node-pairs
        that can be used to build every element in the mesh.

        Args:
            nodes_csv: The filepath to the CSV containing node definitions

        Returns:
            A list of tuples that contain the node ID to target and a list
            of parsed nodes.
        """

        nodes: list[Node] = []

        # Read CSV file, extract does
        with open(nodes_csv, "r") as f:
            header_skipped = False

            id = 0
            for line in f.readlines():
                if not header_skipped:
                    header_skipped = True
                    continue

                items = line.split(",")

                try:
                    x = float(items[0])
                    y = float(items[1])
                    ux = try_float(items[2])
                    uy = try_float(items[3])
                    fx = try_float(items[4])
                    fy = try_float(items[5])
                except Exception:
                    print("Error loading nodes from nodes.csv!")
                    raise

                nodes.append(Node(x, y, ux, uy, fx, fy, id))
                id += 1

        points = np.empty((len(nodes), 2))

        for i, node in enumerate(nodes):
            points[i] = [node.x, node.y]

        tri = Delaunay(points)

        elements = []

        for simp in tri.simplices:
            pairs = combinations(simp, 2)

            for pair in pairs:

                has_match = False
                for existing_pair in elements:
                    if self._pair_equals(pair, existing_pair):
                        print(f"Found duplicate: {pair} == {existing_pair}")
                        has_match = True

                if has_match:
                    break
                else:
                    elements.append(pair)

        return elements, nodes
