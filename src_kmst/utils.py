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

    def find_tree(self):
        start = 1
        tree = {start}
        visited = {i: False for i in self.V}
        q = [start]
        while len(q) > 0:
            u = q.pop(0)
            tree.add(u)
            if len(tree) == self.k:
                return tree
            visited[u] = True
            for v in self.V:
                if (u, v) in self.E and not visited[v]:
                    q.append(v)
        return tree




OPTIMAL_SOLUTIONS = {
    '1_0': 46, '1_1': 477, '2_0': 373, '2_1': 1390,
    '3_0': 725, '3_1': 3074, '4_0': 909, '4_1': 3292,
    '5_0': 1235, '5_1': 4898, '6_0': 2068, '6_1': 6705,
    '7_0': 1335, '7_1': 4534, '8_0': 1620, '8_1': 5787,
    '9_0': 2289, '9_1': 7595, '10_0': 4182, '10_1': 14991,
}