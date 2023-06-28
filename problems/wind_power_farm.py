import numpy as np
from ortools.linear_solver import pywraplp


def f():
    solver = pywraplp.Solver.CreateSolver('SCIP')

    L = {1, 2, 3, 4}
    x, z = {}, {}
    for i in L:
        x[i] = solver.BoolVar(f'x_{i}')
        z[i] = solver.BoolVar(f'z_{i}')

    solver.Add(20 * sum(x[i] for i in L) + 20 * sum(z[i] for i in L) <= 50)
    for i in L:
        solver.Add(x[i] + z[i] <= 1)
    neighbor = {1: [2, 3], 2: [1, 3], 3: [1, 2, 4], 4: [3]}
    for i in L:
        solver.Add(z[i] <= sum(x[j] for j in neighbor[i]) + sum(z[j] for j in neighbor[i]))
        solver.Add((1 - x[i]) * 100 >= sum(x[j] for j in neighbor[i]) + sum(z[j] for j in neighbor[i]))

    e = {
        1: 20,
        2: 24,
        3: 40,
        4: 20
    }
    solver.Maximize(sum(e[i] * x[i] for i in L) + 0.75 * sum(e[i] * z[i] for i in L))

    solver.Solve()
    print('Objective value =', solver.Objective().Value())
    print('Vars x =', [x[i].solution_value() for i in L])
    print('Vars z =', [z[i].solution_value() for i in L])


def g():
    solver = pywraplp.Solver.CreateSolver('SCIP')

    D = {1, 2, 3, 4, 5}
    L = {'a', 'b', 'c', 'd', 'e'}
    x, y = {}, {}
    for i in L:
        x[i] = solver.BoolVar(f'x_{i}')
    for i in D:
        y[i] = solver.BoolVar(f'y_{i}')

    costs = {
        'a': 20,
        'b': 40,
        'c': 30,
        'd': 10,
        'e': 50
    }
    near = {
        'a': [1, 5],
        'b': [1, 3, 4],
        'c': [2, 3],
        'd': [5],
        'e': [3, 4]
    }
    inverse_near = {
        1: ['a', 'b'],
        2: ['c'],
        3: ['b', 'c', 'e'],
        4: ['b', 'e'],
        5: ['a', 'd']
    }

    solver.Add(sum(y[j] for j in D) >= 4)

    for i in L:
        for j in near[i]:
            solver.Add(x[i] <= y[j])

    for j in D:
        solver.Add(y[j] <= sum(x[i] for i in L if i in inverse_near[j]))

    solver.Minimize(sum(costs[i] * x[i] for i in L))

    # Print solution
    solver.Solve()
    print('Objective value =', solver.Objective().Value())
    print('Vars x =', [x[i].solution_value() for i in L])
    print('Vars y =', [y[i].solution_value() for i in D])


if __name__ == "__main__":
    # f()
    g()
