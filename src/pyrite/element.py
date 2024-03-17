from pyrite.datatypes import Node, MATERIAL_ELASTICITY, CROSS_AREA, MatrixIndex, Axis
import numpy as np
import math


class Element:

    def __init__(self, n1: Node, n2: Node, n3: Node):
        self.n1: Node = n1
        self.n2: Node = n2
        self.n3: Node = n3

        self.n1_u = (0.0, 0.0)
        self.n2_u = (0.0, 0.0)
        self.n3_u = (0.0, 0.0)

        self._global_stiffness_matrix = None

        self.row_indexes = []
        self.column_indexes = []

        self._normal_stress = 0.0
        self._normal_strain = 0.0

    @property
    def length(self) -> float:
        """Calculates the length of the element"""

        dx = abs(self.n1.x - self.n2.x)
        dy = abs(self.n1.y - self.n2.y)

        return math.sqrt(dx**2 + dy**2)

    @property
    def local_stiffness_matrix(self) -> np.ndarray:
        """Calculates the local stiffness matrix"""

        line_element_matrix = np.array(
            [[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]]
        )

        return line_element_matrix * ((MATERIAL_ELASTICITY * CROSS_AREA) / self.length)

    @property
    def element_angle(self) -> float:
        """Calculates the angle of the element relative to the global
        coordinate system. Returns a float in the interval [0, 360)"""

        dx = self.n1.x - self.n2.x
        dy = self.n1.y - self.n2.y

        if dx < 0:
            dx *= -1
            dy *= -1

        return math.degrees(math.atan2(dy, dx)) % 360

    @property
    def transition_matrix(self) -> np.ndarray:
        """Computes the transition matrix for the global stiffness matrix"""

        t = math.radians(self.element_angle)

        C = math.cos(t)
        S = math.sin(t)

        # print(((MATERIAL_ELASTICITY * CROSS_AREA)/self.length))

        return np.array(
            [
                [C**2, C * S, -(C**2), -C * S],
                [C * S, S**2, -C * S, -(S**2)],
                [-(C**2), -C * S, C**2, C * S],
                [-C * S, -(S**2), C * S, S**2],
            ]
        )

    def _calculate_global_stiffness_matrix(self) -> np.ndarray:
        """Computes the global stiffness matrix for the element"""

        self.row_indexes = [
            MatrixIndex(self.n1, Axis.X),
            MatrixIndex(self.n1, Axis.Y),
            MatrixIndex(self.n2, Axis.X),
            MatrixIndex(self.n2, Axis.Y),
        ]

        self.column_indexes = [
            MatrixIndex(self.n1, Axis.X),
            MatrixIndex(self.n1, Axis.Y),
            MatrixIndex(self.n2, Axis.X),
            MatrixIndex(self.n2, Axis.Y),
        ]

        # return self.transition_matrix  # TODO <-- Delete
        return self.transition_matrix * (
            (MATERIAL_ELASTICITY * CROSS_AREA) / self.length
        )

    @property
    def global_stiffness_matrix(self) -> np.ndarray:
        """Returns the global stiffness matrix for the element"""

        if self._global_stiffness_matrix is None:
            self._global_stiffness_matrix = self._calculate_global_stiffness_matrix()

        return self._global_stiffness_matrix

    def index_GSM(self, row: MatrixIndex, column: MatrixIndex):

        self._calculate_global_stiffness_matrix()

        if row in self.row_indexes and column in self.column_indexes:

            return self.global_stiffness_matrix[
                self.row_indexes.index(row), self.column_indexes.index(column)
            ]

        else:
            return 0

    @property
    def normal_strain(self) -> float:
        """Returns the normal strain of the element. This is usually
        set after solving by the post-processor.
        """
        return self._normal_strain

    @normal_strain.setter
    def normal_strain(self, value: float) -> None:
        """Sets the normal strain in the element. This should only be
        set by the post-processor.
        """

        self._normal_strain = value

    @property
    def normal_stress(self) -> float:
        """Returns the average stress from the calculated strain. Note,
        if the strain attribute was not set, this will default to zero.

        Returns:
            The stress, in whichever unit MATERIAL_ELASTICITY is set.
        """

        return MATERIAL_ELASTICITY * self.normal_strain
