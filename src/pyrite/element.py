"""
Element class for each triangular planar element. Contains code to compute
strain-displacement and stress-strain matrices for each element, and 
provides tools for indexing the element's global stiffness matrix.

March 25, 2024
Kyle Tennison
"""

from pyrite.datatypes import Node, MatrixIndex, Axis
import numpy as np


class Element:

    poisson_ratio: float
    part_thickness: float
    material_elasticity: int

    def __init__(self, n1: Node, n2: Node, n3: Node):
        self.n1: Node = n1
        self.n2: Node = n2
        self.n3: Node = n3

        self._check_ccw()

        self.n1_u = (0.0, 0.0)
        self.n2_u = (0.0, 0.0)
        self.n3_u = (0.0, 0.0)

        self._global_stiffness_matrix = None

        self.row_indexes = []
        self.column_indexes = []

        self.stress = 0.0

    def _check_ccw(self) -> None:
        """2D elements must be arranged counter-clockwise. Reassigns
        nodes 1, 2, and 3 if needed.
        """

        if self.area < 0:
            self.n1, self.n3 = (self.n3, self.n1)

    @property
    def area(self) -> float:
        """Returns the area of the triangular element"""

        return 0.5 * (
            self.n1.x * (self.n2.y - self.n3.y)
            + self.n2.x * (self.n3.y - self.n1.y)
            + self.n3.x * (self.n1.y - self.n2.y)
        )

    @property
    def strain_displacement_matrix(self) -> np.ndarray:
        """Calculates the strain-displacement matrix of the element"""

        beta_1 = self.n2.y - self.n3.y
        beta_2 = self.n3.y - self.n1.y
        beta_3 = self.n1.y - self.n2.y

        gamma_1 = self.n3.x - self.n2.x
        gamma_2 = self.n1.x - self.n3.x
        gamma_3 = self.n2.x - self.n1.x

        B = np.array(
            [
                [beta_1, 0, beta_2, 0, beta_3, 0],
                [0, gamma_1, 0, gamma_2, 0, gamma_3],
                [gamma_1, beta_1, gamma_2, beta_2, gamma_3, beta_3],
            ],
            dtype=float,
        )

        # Prevent zero division
        if self.area == 0:
            B *= 0
        else:
            B *= 1 / (2 * self.area)

        return B

    @property
    def stress_strain_matrix(self) -> np.ndarray:
        """Calculates the stress-strain matrix for the planar element"""

        D = np.array(
            [
                [1, self.poisson_ratio, 0],
                [self.poisson_ratio, 1, 0],
                [0, 0, (1 - self.poisson_ratio) / 2],
            ]
        )

        D *= self.material_elasticity / (1 - self.poisson_ratio**2)

        return D

    def _calculate_global_stiffness_matrix(self) -> np.ndarray:
        """Computes the global stiffness matrix for the element"""

        self.row_indexes = [
            MatrixIndex(self.n1, Axis.X),
            MatrixIndex(self.n1, Axis.Y),
            MatrixIndex(self.n2, Axis.X),
            MatrixIndex(self.n2, Axis.Y),
            MatrixIndex(self.n3, Axis.X),
            MatrixIndex(self.n3, Axis.Y),
        ]

        self.column_indexes = [
            MatrixIndex(self.n1, Axis.X),
            MatrixIndex(self.n1, Axis.Y),
            MatrixIndex(self.n2, Axis.X),
            MatrixIndex(self.n2, Axis.Y),
            MatrixIndex(self.n3, Axis.X),
            MatrixIndex(self.n3, Axis.Y),
        ]

        return (
            (
                (
                    self.strain_displacement_matrix.transpose()
                    @ self.stress_strain_matrix
                )
                @ self.strain_displacement_matrix
            )
            * self.area
            * self.part_thickness
        )

    @property
    def global_stiffness_matrix(self) -> np.ndarray:
        """Returns the global stiffness matrix for the element"""

        if self._global_stiffness_matrix is None:
            self._global_stiffness_matrix = self._calculate_global_stiffness_matrix()

        return self._global_stiffness_matrix

    def index_GSM(self, row: MatrixIndex, column: MatrixIndex):
        """Index the element's global stiffness matrix"""

        self.global_stiffness_matrix

        # print(f"debug: indexing element {self} at {row.node.index}, {column.node.index}")

        if row in self.row_indexes and column in self.column_indexes:

            return self.global_stiffness_matrix[
                self.row_indexes.index(row), self.column_indexes.index(column)
            ]

        else:
            return 0

    def __repr__(self) -> str:
        return f"Element({self.n1.index}, {self.n2.index}, {self.n3.index})"
