"""
Solves the finite element problem once elements and nodes have been generated.
Assumes linear elastic behavior on an isotropic material.

March 25, 2023
Kyle Tennison
"""

from pyrite.datatypes import Node, DOF
from pyrite.element import Element
from pyrite.post_processor import PostProcessor

from typing import Optional
import math
import numpy as np


def magnitude(vector: np.ndarray) -> float:
    """Calculates the magnitude of a np array"""

    return math.sqrt(vector[0] ** 2 + vector[1] ** 2)


class Solver:
    """
    This class contains tools to solve the FEM equations.
    """

    def __init__(self):
        self.nodes: list[Node] = []
        self.elements: list[Element] = []
        self.matrix_size = DOF * len(self.nodes)

    def create_total_stiffness_matrix(
        self, nodes: list[Node], elements: list[Element]
    ) -> np.ndarray:
        """Creates a total stiffness matrix from each global stiffness
        matrix on each element.

        Returns:
            A 2D matrix containing the matrix
        """

        n = DOF * len(nodes)

        total_stiffness_matrix = np.zeros((n, n))

        print("info: building total stiffness matrix:")
        for i, element in enumerate(elements):
            print(f"{100*(i/(len(elements)-1)):.3f}%")

            element.global_stiffness_matrix

            for row_index in element.row_indexes:
                for col_index in element.column_indexes:

                    row_num = row_index.to_num()
                    col_num = col_index.to_num()

                    value = element.index_GSM(row_index, col_index)

                    total_stiffness_matrix[row_num, col_num] += value

        return total_stiffness_matrix

    def calculate_force_vector(self, nodes: list[Node]) -> np.ndarray:
        """Calculates column vector of forces."""
        F = np.zeros((len(nodes) * DOF))

        i = 0
        for node in nodes:
            F[i] = node.Fx
            F[i + 1] = node.Fy
            i += 2

        return F

    def calculate_displacement_vector(self, nodes: list[Node]) -> np.ndarray:
        """Calculates column vector of displacements."""
        U = np.zeros((len(nodes) * DOF))

        i = 0
        for node in nodes:
            U[i] = node.ux
            U[i + 1] = node.uy
            i += 2

        return U

    def check_unconstrained(self, f: np.ndarray, u: np.ndarray) -> None:
        """Raises an error if the model is over or under constrained"""

        for i in range(len(f)):
            if np.isnan(f[i]) and np.isnan(u[i]):
                raise Exception(f"Model Underconstrained at index {i}")

            if (not np.isnan(f[i])) and (not np.isnan(u[i])):
                raise Exception(f"Model Overconstrained at index {i}")

    def count_unknown(self, column_vector: np.ndarray) -> tuple[int, int]:
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

    def display_total_matrix(
        self,
        nodal_forces: np.ndarray,
        nodal_displacements: np.ndarray,
        total_stiffness_matrix: np.ndarray,
    ) -> None:
        """Displays the total matrix equation in a readable format."""

        for i in range(len(total_stiffness_matrix)):
            print(
                "| {f:<6} |   {k:<10}  | {u:0.3e} | ".format(
                    f=round(nodal_forces[i], 3),
                    k=str(
                        [
                            f"{i:0.2e}" if i < 0 else f"{i:0.3e}"
                            for i in total_stiffness_matrix[i]
                        ]
                    ),
                    u=nodal_displacements[i],
                )
            )

    def assemble_known_and_unknown_matrices(
        self,
        num_known_displacements: int,
        num_unknown_displacements: int,
        nodal_forces: np.ndarray,
        nodal_displacements: np.ndarray,
        total_stiffness_matrix: np.ndarray,
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
        unknown_matrix = np.empty(
            (num_unknown_displacements, num_unknown_displacements)
        )
        known_matrix = np.empty((num_unknown_displacements, num_known_displacements))
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
                        total_stiffness_matrix[tsm_row, column]
                        * nodal_displacements[column]
                    )

            known_matrix[local_row] = k_buffer
            unknown_matrix[local_row] = uk_buffer

            local_row += 1

        return known_matrix, unknown_matrix

    @staticmethod
    def redundant_solve(A, b) -> Optional[np.ndarray]:
        """
        Solves the linear system of equations Ax = b using the pseudoinverse.

        Args:
            A: The coefficient matrix
            b: The column vector solution

        Returns:
            A numpy array representing the solution to the system,
            or None if the system has no solution.
        """
        try:
            # Use np.linalg.solve for non-singular matrices
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            # If singular, use the pseudoinverse
            return np.linalg.pinv(A) @ b

    def solve(
        self,
        nodal_forces: np.ndarray,
        nodal_displacements: np.ndarray,
        total_stiffness_matrix: np.ndarray,
    ) -> None:
        """Solves for unknown forces and nodal displacements.

        Returns:
            Nothing. nodal_forces and nodal_displacements arrays will be updated
        """

        self.check_unconstrained(nodal_forces, nodal_displacements)

        num_known_displacements, num_unknown_displacements = self.count_unknown(
            nodal_displacements
        )

        # Assemble known and unknown matrices
        known_matrix, unknown_matrix = self.assemble_known_and_unknown_matrices(
            num_known_displacements,
            num_unknown_displacements,
            nodal_forces,
            nodal_displacements,
            total_stiffness_matrix,
        )

        known_forces = np.array([i for i in nodal_forces if not np.isnan(i)]).reshape(
            num_unknown_displacements, 1
        )

        known_matrix_summed = np.array([sum(i) for i in known_matrix]).reshape(
            (num_unknown_displacements, 1)
        )

        displacement_solution = self.redundant_solve(
            unknown_matrix, (known_forces + known_matrix_summed)
        )

        # print("known forces:\n", known_forces)
        # print("known matrix summed:\n", known_matrix_summed)

        if displacement_solution is None:
            raise Exception("No solution.")

        solution_cursor = 0
        for i, u in enumerate(nodal_displacements):
            if np.isnan(u):
                nodal_displacements[i] = displacement_solution[solution_cursor][0]
                solution_cursor += 1

        for i, f in enumerate(nodal_forces):
            if np.isnan(f):

                solved_force = 0

                for c in range(len(total_stiffness_matrix)):
                    solved_force += (
                        total_stiffness_matrix[i, c] * nodal_displacements[c]
                    )

                nodal_forces[i] = solved_force

        # print("The solved matrix is:")
        # self.display_total_matrix(
        #     nodal_forces, nodal_displacements, total_stiffness_matrix
        # )

    def run(self, nodes: list[Node], elements: list[Element]):
        """Runs the FEA simulation for a set of given nodes and elements.
        Runs post-processor when complete.

        Args:
            nodes: A list of nodes
            elements: A list of tuples that contain the node indexes to target for
                each element.
        """

        # Load nodes and elements
        self.nodes = nodes
        self.elements = elements

        print("info: loaded the following elements")
        for element in elements:
            print(f"\t- {element}")

        print(f"\ninfo: {len(elements)} total elements")

        # Create total stiffness matrix from nodes and elements
        print("info: assembling total stiffness matrix")
        total_stiffness_matrix = self.create_total_stiffness_matrix(
            self.nodes, self.elements
        )

        # Create nodal forces/displacement vectors
        print("info: assembling column vectors")
        nodal_forces = self.calculate_force_vector(self.nodes)
        nodal_displacements = self.calculate_displacement_vector(self.nodes)

        # Display pre-solved matrix
        # print("info: the total matrix is:")
        # self.display_total_matrix(
        #     nodal_forces, nodal_displacements, total_stiffness_matrix
        # )

        # Solve matrix
        try:
            print("info: solving system...")
            self.solve(nodal_forces, nodal_displacements, total_stiffness_matrix)
        except Exception as e:
            print(f"error: solve failed: {e}")
            raise

        # Load nodal displacements into each element object

        post_processor = PostProcessor()
        post_processor.load_displacements_into_elements(
            nodal_displacements, self.elements
        )
        post_processor.compute_stress(self.elements)
        post_processor.output_solved(nodal_forces, nodal_displacements, nodes)
        post_processor.show("nodes.csv", self.elements)
