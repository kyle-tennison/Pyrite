"""
Miscellaneous datatypes used in Pyrite

Kyle Tennison
March 20, 2024
"""

from dataclasses import dataclass
from typing import Optional
from enum import Enum

DOF = 2
@dataclass 
class PartMetadata:
    """Contains metadata defined in input json"""
    material_elasticity: int
    poisson_ratio: float
    part_thickness: float

@dataclass
class Node:
    """Represents a node in the mesh"""
    x: float
    y: float
    ux: Optional[float]
    uy: Optional[float]
    Fx: Optional[float]
    Fy: Optional[float]
    index: int

    def __eq__(self, other: "Node") -> bool:
        return self.index == other.index


class Axis(Enum):
    """Different axes available"""
    X = 0
    Y = 1
    Z = 2

class MshState(Enum):
    """Different states the MeshParser may be in"""
    NODES = 0
    ELEMENTS = 1
    ENTITIES = 2
    LIMBO = 3


@dataclass
class MatrixIndex:
    """Used to index a stiffness matrix given a node and axis"""
    node: Node
    axis: Axis

    def to_num(self) -> int:
        """Converts the index to the corresponding number in the total 
        stiffness matrix.
        """
        return (2 * self.node.index) + (1 if self.axis == Axis.Y else 0)
