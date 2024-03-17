from dataclasses import dataclass
from enum import Enum
import math
from typing import Optional

import numpy as np

CROSS_AREA = 10
MATERIAL_ELASTICITY = 200e6
v = 0.3
DOF = 2


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
    X = "X"
    Y = "Y"


class MshState(Enum):
    NODES = 0
    ELEMENTS = 1
    ENTITIES = 2
    LIMBO = 3


@dataclass
class MatrixIndex:
    node: Node
    axis: Axis


@dataclass
class ElementLight:
    n1: Node
    n2: Node
