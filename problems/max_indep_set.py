import numpy as np
from ortools.linear_solver import pywraplp
import networkx as nx
import matplotlib.pyplot as plt


def f():
    solver = pywraplp.Solver.CreateSolver('SCIP')
    V = set(i + 1 for i in range(6))
    E = [(1,2), (1,4), (1,3), (1,6), (2,4), (2,5), (3,4), (3,6), (4,5), (4,6), (5,6)]

    G = nx.Graph()
    G.add_nodes_from(V)
    G.add_edges_from(E)

    # Draw the graph
    cliques = nx.find_cliques(G)

    x = {}
    for i in V:
        x[i] = solver.IntVar(0, 1, f'x_{i}')

    for clique in cliques:
        solver.Add(sum(x[i] for i in clique) <= 1)

    solver.Maximize(sum(x[i] for i in V))

    solver.Solve()

    print('Objective value =', solver.Objective().Value())
    print('Vars =', [x[i].solution_value() for i in V])


if __name__ == "__main__":
    f()

