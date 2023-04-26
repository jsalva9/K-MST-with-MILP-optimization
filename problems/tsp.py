import pandas as pd
from ortools.linear_solver import pywraplp


class TSP:
    def __init__(self, n, tol=0.1):
        self.n = n
        self.tol = tol
        self.solver = pywraplp.Solver.CreateSolver('SCIP')
        self.df_dist = pd.read_csv('../data/df_dist.csv')

        self.df_dist = self.df_dist[self.df_dist['city_id1'] < n]

        self.d1 = {(row['city_id1'] + 1, row['city_id2'] + 1): row['dist'] for i, row in self.df_dist.iterrows()}
        self.d2 = {(row['city_id2'] + 1, row['city_id1'] + 1): row['dist'] for i, row in self.df_dist.iterrows()}
        self.d = {**self.d1, **self.d2}

        self.n = len(self.df_dist['city_id1'].unique())
        self.cities = [i + 1 for i in range(self.n)]

        self.x = {}
        self.u = {}

        self.M = self.n

    def define_variables(self):
        for i in self.cities:
            for j in self.cities:
                if i == j:
                    continue
                self.x[i, j] = self.solver.IntVar(0, 1, f'x_{i}_{j}')
        for i in self.cities:
            if i == 1:
                continue
            self.u[i] = self.solver.IntVar(1, self.n - 1, f'u_{i}')

    def define_constraints(self):
        for i in self.cities:
            self.solver.Add(sum(self.x[i, j] for j in self.cities if j != i) == 1)
            self.solver.Add(sum(self.x[j, i] for j in self.cities if j != i) == 1)

        for i in self.cities:
            for j in self.cities:
                if i == j or i == 1 or j == 1:
                    continue
                self.solver.Add(self.u[i] + self.x[i, j] <= self.u[j] + self.n * (1 - self.x[i, j]))

    def define_objective(self):
        self.solver.Minimize(sum(self.d[i, j] * self.x[i, j] for i in self.cities for j in self.cities if i != j))

    def define_tolerance(self):
        parameters = pywraplp.MPSolverParameters()
        parameters.SetDoubleParam(pywraplp.MPSolverParameters.RELATIVE_MIP_GAP, self.tol)
        return parameters

    def run(self):
        self.define_variables()
        self.define_constraints()
        self.define_objective()
        parameters = self.define_tolerance()
        print('Number of variables =', self.solver.NumVariables())
        print('Number of constraints =', self.solver.NumConstraints())
        print('Solving the problem ...')
        status = self.solver.Solve(parameters)
        self.display_result(status)

    def display_result(self, status):
        print(f'Solving time = {self.solver.WallTime()} milliseconds')
        if status == pywraplp.Solver.OPTIMAL:
            print('Objective value =', self.solver.Objective().Value())
            for i in self.cities:
                for j in self.cities:
                    if i == j:
                        continue
                    if self.x[i, j].solution_value() == 1:
                        print(f'x_{i}_{j} = {self.x[i, j].solution_value()}')
        else:
            print('The problem does not have an optimal solution.')


if __name__ == '__main__':
    tsp = TSP(n=2000, tol=1.5)
    tsp.run()
