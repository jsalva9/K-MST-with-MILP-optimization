import numpy as np
from ortools.linear_solver import pywraplp


def f(u):
    solver = pywraplp.Solver.CreateSolver('GUROBI')
    y_1 = solver.IntVar(0, 1, 'y_1')
    y_2 = solver.IntVar(0, 1, 'y_2')
    y_3 = solver.IntVar(0, 1, 'y_3')
    # u = solver.NumVar(-solver.infinity(), solver.infinity(), 'u')
    # solver.Add(3 * y_1 + y_2 + 4 * y_3 <= 4)
    solver.Maximize(10 * y_1 + 4 * y_2 + 14 * y_3 + u * (4 - 3 * y_1 - y_2 - 4 * y_3))

    solver.Solve()
    print('Objective value =', solver.Objective().Value())
    print('Vars =', y_1.solution_value(), y_2.solution_value(), y_3.solution_value())
    return round(y_1.solution_value()), round(y_2.solution_value()), round(y_3.solution_value()), \
        solver.Objective().Value()


def g():
    map = {}
    for u in np.linspace(0, 10, 101):
        print(u)
        map[round(u, 4)] = round(f(u)[3], 4)
    # visualize in a graph
    import matplotlib.pyplot as plt
    plt.plot(map.keys(), map.values())
    # Add title and axis names
    plt.title('Objective value vs. u')
    plt.xlabel('u')
    plt.ylabel('max z(u)')
    plt.show()
    print(map)


def subgradient():
    u = 1
    k = 0
    mu0 = 1
    rho = 0.9
    # Write results to a file
    with open('results.txt', 'w') as file:
        while k <= 200:
            mu = mu0 * rho ** k
            file.write(f'{k} & {round(u, 5)} & {mu:.2e} ')
            y1, y2, y3, opt = f(u)
            file.write(f'& {y1} & {y2} & {y3} & {round(opt, 5)} \\\\ \n')
            u = max(0, u - mu * (4 - 3 * y1 - y2 - 4 * y3))
            k += 1


if __name__ == '__main__':
    # f(3.5)
    # subgradient()
    g()



