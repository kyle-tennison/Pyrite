from copy import deepcopy
from dataclasses import dataclass
from enum import Enum
import math
from typing import Optional
import numpy as np

CROSS_AREA = 10
MATERIAL_ELASTICITY = 200E6
v = 0.3
DOF = 2


def magnitude(vector: np.ndarray) -> float:
    """Calculates the magnitude of a np array"""

    return math.sqrt(vector[0] ** 2 + vector[1] ** 2)


@dataclass
class Node:
    x: float
    y: float
    ux: Optional[float]
    Fx: Optional[float]
    uy: Optional[float]
    Fy: Optional[float]


class Axis(Enum):
    X = 'X'
    Y = 'Y'


@dataclass
class MatrixIndex:
    node: Node
    axis: Axis


class Element:

    def __init__(self, n1: Node, n2: Node):
        self.n1: Node = n1
        self.n2: Node = n2

        self._global_stiffness_matrix = None

        self.row_indexes = []
        self.column_indexes = []

    @property
    def length(self) -> float:
        """Calculates the length of the element"""

        dx = abs(self.n1.x - self.n2.x)
        dy = abs(self.n1.y - self.n2.y)

        return math.sqrt(dx ** 2 + dy ** 2)

    @property
    def local_stiffness_matrix(self) -> np.ndarray:
        """Calculates the local stiffness matrix"""

        line_element_matrix = np.array(
            [
                [1, 0, -1, 0],
                [0, 0, 0, 0],
                [-1, 0, 1, 0],
                [0, 0, 0, 0]
            ]
        )

        return line_element_matrix * ((MATERIAL_ELASTICITY * CROSS_AREA)/self.length)

    @property
    def element_angle(self) -> float:
        """Calculates the angle of the element relative to the global 
        coordinate system. Returns a float in the interval [0, 360) """

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
            [[C ** 2, C*S, - C**2, -C*S],
             [C*S, S**2, -C*S, -S**2],
             [-C**2, -C*S, C**2, C*S],
             [-C*S, -S**2, C*S, S**2]]
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
        return self.transition_matrix * ((MATERIAL_ELASTICITY * CROSS_AREA)/self.length)

    @property
    def global_stiffness_matrix(self) -> np.ndarray:
        """Returns the global stiffness matrix for the element"""

        if self._global_stiffness_matrix is None:
            self._global_stiffness_matrix = self._calculate_global_stiffness_matrix()

        return self._global_stiffness_matrix

    def index_GSM(self, row: MatrixIndex, column: MatrixIndex):

        self._calculate_global_stiffness_matrix()

        if row in self.row_indexes and column in self.column_indexes:

            return self.global_stiffness_matrix[self.row_indexes.index(row), self.column_indexes.index(column)]

        else:
            return 0


def create_total_stiffness_matrix(nodes, elements) -> np.ndarray:

    total_stiffness_matrix = np.zeros((DOF * len(nodes), DOF * len(nodes)))

    node_stack = []

    for node in nodes.values():
        node_stack.append(MatrixIndex(node, Axis.X))
        node_stack.append(MatrixIndex(node, Axis.Y))

    row_stack = node_stack.copy()
    row = 0
    while row_stack:
        n1 = row_stack.pop(0)

        column_stack = node_stack.copy()
        column = 0
        while column_stack:
            n2 = column_stack.pop(0)

            ki = 0
            for element in elements.values():
                ki += element.index_GSM(n1, n2)

            total_stiffness_matrix[row, column] = ki

            column += 1

        row += 1

    return total_stiffness_matrix


def calculate_force_vector(nodes: dict[int, Node]) -> np.ndarray:
    """Calculates column vector of forces."""
    F = np.zeros((len(nodes) * DOF))

    i = 0
    for node in nodes.values():
        F[i] = node.Fx
        F[i + 1] = node.Fy
        i += 2

    return F


def calculate_displacement_vector(nodes: dict[int, Node]) -> np.ndarray:
    """Calculates column vector of displacements."""
    U = np.zeros((len(nodes) * DOF))

    i = 0
    for node in nodes.values():
        U[i] = node.ux
        U[i + 1] = node.uy
        i += 2

    return U


def check_unconstrained(f: np.ndarray, u: np.ndarray) -> None:
    """Raises an error if the model is over or under constrained"""

    for i in range(len(f)):
        if np.isnan(f[i]) and np.isnan(u[i]):
            raise Exception(f"Model Underconstrained at index {i}")

        if (not np.isnan(f[i])) and (not np.isnan(u[i])):
            raise Exception(f"Model Overconstrained at index {i}")


def count_unknown(column_vector: np.ndarray) -> tuple[int, int]:
    """Counts the number of NaN values in a vector.
    
    Returns:
        Tuple of (# known, # unknown)
    """

    unknowns = 0
    for i in column_vector:
        if np.isnan(i):
            unknowns += 1

    knowns = len(column_vector) - unknowns

    return (knowns, unknowns)


def display_total_matrix(nodal_forces, nodal_displacements, total_stiffness_matrix) -> None:

    for i in range(len(total_stiffness_matrix)):
        print("| {f:<6} |   {k:<10}  | {u:0.3e} | ".format(f=round(nodal_forces[i], 3), k=str(
            [f'{i:0.2e}' if i < 0 else f'{i:0.3e}' for i in total_stiffness_matrix[i]]), u=nodal_displacements[i]))


def assemble_known_and_unknown_matrices(
        num_known_displacements: int, 
        num_unknown_displacements: int, 
        nodal_forces: np.ndarray, 
        nodal_displacements: np.ndarray, 
        total_stiffness_matrix: np.ndarray
        ) -> tuple[np.ndarray, np.ndarray]:
    """Unknown/Known Matrices are used to solve for displacements in a 
    linear system. They point to partitions of the total stiffness matrix.
    This function will generate the matrix of total stiffness matrix elements
    that correspond to unknown displacements, and it will do the same for
    known displacements.
    
    Returns:
        A tuple of (known_matrix, unknown_matrix)
    """

    # init empty matrices
    unknown_matrix=np.empty((num_unknown_displacements, num_unknown_displacements))
    known_matrix=np.empty((num_known_displacements, num_known_displacements))
    local_row = 0

    # iterate for each row in TSM. Filter TSM indices that correspond to a 
    # known displacement versus the ones that correspond to an unknown displacement.
    for tsm_row in range(len(nodal_forces)):
        if np.isnan(nodal_forces[tsm_row]):
            continue

        
        uk_buffer = []
        k_buffer = []
        for column in range(len(nodal_forces)):
            if np.isnan(nodal_displacements[column]):
                uk_buffer.append(total_stiffness_matrix[tsm_row, column])
            else:
                k_buffer.append(
                    total_stiffness_matrix[tsm_row, column] * nodal_displacements[column])

        known_matrix[local_row] = k_buffer
        unknown_matrix[local_row] = uk_buffer

        local_row += 1

    return known_matrix, unknown_matrix


def solve(nodal_forces: np.ndarray, nodal_displacements: np.ndarray, total_stiffness_matrix: np.ndarray):
    """Solves for unknown forces and nodal displacements.
    
    Returns:
        Nothing. nodal_forces and nodal_displacements arrays will be updated
    """

    check_unconstrained(nodal_forces, nodal_displacements)

    num_known_displacements, num_unknown_displacements = count_unknown(nodal_displacements)

    # Assemble known and unknown matrices
    known_matrix, unknown_matrix = assemble_known_and_unknown_matrices(
        num_known_displacements,
        num_unknown_displacements,
        nodal_forces,
        nodal_displacements,
        total_stiffness_matrix
    )

    known_forces = np.array([i for i in nodal_forces if not np.isnan(i)]).reshape(3, 1)
    known_matrix_summed = np.array([sum(i) for i in known_matrix]).reshape((num_known_displacements,1))
    
    displacement_solution = np.linalg.solve(unknown_matrix, (known_forces + known_matrix_summed))

    solution_cursor = 0
    for i, u in enumerate(nodal_displacements):
        if np.isnan(u):
            nodal_displacements[i] = displacement_solution[solution_cursor][0]
            solution_cursor += 1


    for i, f in enumerate(nodal_forces):
        if np.isnan(f):

            solved_force = 0

            for c in range(len(total_stiffness_matrix)):
                solved_force += total_stiffness_matrix[i, c] * nodal_displacements[c]
                

            nodal_forces[i] = solved_force


    print("The solved matrix is:")
    display_total_matrix(nodal_forces, nodal_displacements, total_stiffness_matrix)

    

def main():

    # Load nodes and elements
    nodes = {
        1: Node(x=0,   y=0,                          ux=None, uy=0,    Fx=0, Fy=None),
        2: Node(x=0.5, y=math.sin(math.radians(60)), ux=None, uy=None, Fx=0,    Fy=-100),
        3: Node(x=1,   y=0,                          ux=0,    uy=0,    Fx=None,    Fy=None),
    }

    elements = {
        1: Element(nodes[1], nodes[2]),
        2: Element(nodes[2], nodes[3]),
        3: Element(nodes[1], nodes[3])
    }

    # Create total stiffness matrix from nodes and elements
    total_stiffness_matrix = create_total_stiffness_matrix(nodes, elements)

    # Create nodal forces/displacement vectors
    nodal_forces = calculate_force_vector(nodes)
    nodal_displacements = calculate_displacement_vector(nodes)

    # Display pre-solved matrix
    print("The total matrix is:")
    display_total_matrix(nodal_forces, nodal_displacements, total_stiffness_matrix)

    solve(nodal_forces, nodal_displacements, total_stiffness_matrix)

    


if __name__ == "__main__":
    main()
