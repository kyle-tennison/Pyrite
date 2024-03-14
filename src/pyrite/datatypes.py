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

    def __eq__(self, other: 'ElementLight'):
        
        match = False 

        if self.n1 == other.n1 and self.n2 == other.n2:
            match = True 

        if self.n1 == other.n2 and self.n2 == other.n1:
            match = True 

        return True