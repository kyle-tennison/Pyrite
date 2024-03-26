from dataclasses import dataclass
from enum import Enum
import math
from typing import Optional

import numpy as np

CROSS_AREA = 10
MATERIAL_ELASTICITY = 30e6
POISSON_RATIO = 0.25  # poisson ratio
DOF = 2
PART_THICKNESS = 0.5


@dataclass
class Node:
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
    X = 0
    Y = 1
    Z = 2


class BoundaryTarget(Enum):
    UX = "ux"
    UY = "uy"
    FX = "fx"
    FY = "fy"


class MshState(Enum):
    NODES = 0
    ELEMENTS = 1
    ENTITIES = 2
    LIMBO = 3


@dataclass
class MatrixIndex:
    node: Node
    axis: Axis

    def to_num(self):
        """Converts the index to the corresponding number in the TSM"""
        return (2 * self.node.index) + (1 if self.axis == Axis.Y else 0)


@dataclass
class ElementLight:
    n1: Node
    n2: Node
