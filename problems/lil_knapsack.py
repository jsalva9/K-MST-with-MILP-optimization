import numpy as np
from ortools.linear_solver import pywraplp


def f(u):
    solver = pywraplp.Solver.CreateSolver('SCIP')
    y_1 = solver.NumVar(0, 1, 'y_1')
    y_2 = solver.NumVar(0, 1, 'y_2')
    y_3 = solver.NumVar(0, 1, 'y_3')
    # u = solver.NumVar(-solver.infinity(), solver.infinity(), 'u')
    # solver.Add(2 * y_1 + 3 * y_2 + y_3 == 4)
    solver.Add(4 * y_1 + 2 * y_2 + 3 * y_3 <= 5)
    u[1] = 0
    solver.Maximize(3 * y_1 + 1 * y_2 + 1 * y_3 + u[0] * (4 - 2 * y_1 - 3 * y_2 - 1 * y_3) + u[1] * (5 - 3 * y_1 - 1 * y_2 - 2 * y_3))

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
    u = np.array([3, 1])
    k = 0
    mu0 = 1
    rho = 1/2
    # Write results to a file
    with open('results.txt', 'w') as file:
        while k <= 200:
            mu = mu0 * rho ** k
            file.write(f'{k} & {(u)} & {mu:.2e} ')
            y1, y2, y3, opt = f(u)
            file.write(f'& {y1} & {y2} & {y3} & {opt} \\\\ \n')

            u[0] = max(0, u[0] - mu * (4 - 2 * y1 - 3 * y2 - 1 * y3))
            u[1] = max(0, u[1] - mu * (5 - 3 * y1 - 1 * y2 - 2 * y3))
            k += 1


if __name__ == '__main__':
    f([0, 0])
    wld = 100
    for i in np.linspace(0, 5, 101):
        for j in np.linspace(0, 5, 101):
            y1, y2, y3, opt = f([i, j])
            wld = min(opt, wld)
    #
    # print(wld)
    # subgradient()
    # g()



