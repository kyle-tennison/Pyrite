from dataclasses import dataclass
from enum import Enum
from typing import Optional


@dataclass
class Node:
    x: float
    y: float
    ux: Optional[float]
    Fx: Optional[float]
    uy: Optional[float]
    Fy: Optional[float]
    id: int

    def __repr__(self) -> str:
        return f"Node({self.id})"
    
    def __eq__(self, other: 'Node') -> bool:
        return self.id == other.id


class Axis(Enum):
    X = 'X'
    Y = 'Y'


@dataclass
class MatrixIndex:
    node: Node
    axis: Axis

@dataclass
class ElementLight:
    n1: Node 
    n2: Node 