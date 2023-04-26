from dataclasses import dataclass
from typing import Set, Tuple, Dict


@dataclass
class Instance:
    name: str
    n: int
    m: int
    Vext: Set[int]
    Eext: Set[Tuple[int, int]]
    weights: Dict[Tuple[int, int], int]
    V: Set[int] = None
    E: Set[Tuple[int, int]] = None
    A: Set[Tuple[int, int]] = None
    Aext: Set[Tuple[int, int]] = None
    k: int = None

    def __post_init__(self):
        self.V = self.Vext - {0}
        self.E = self.Eext - {(0, v) for v in self.V} - {(v, 0) for v in self.V}
        self.A = {(u, v) for (u, v) in self.E}.union({(v, u) for (u, v) in self.E})
        self.Aext = {(u, v) for (u, v) in self.Eext}.union({(v, u) for (u, v) in self.Eext})
        self.weights = {**self.weights, **{(v, u): w for (u, v), w in self.weights.items()}}
