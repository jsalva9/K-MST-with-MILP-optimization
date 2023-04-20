from ortools.linear_solver import pywraplp
import pandas as pd


class SportsLeagueOpti:
    def __init__(self, n, k):
        self.n = n
        self.teams = [i + 1 for i in range(n)]
        self.k = k
        self.solver = pywraplp.Solver.CreateSolver('SAT')
        self.x, self.y, self.z, self.w, self.p = {}, {}, {}, {}, {}
        self.M = 1e5

    def define_variables(self):
        # Define x, y, z
        for i in self.teams:
            for j in self.teams:
                if i == j:
                    continue
                self.x[i, j] = self.solver.IntVar(0, 1, f'x_{i}_{j}')
                self.y[i, j] = self.solver.IntVar(0, 1, f'y_{i}_{j}')
                self.z[i, j] = self.solver.IntVar(0, 1, f'z_{i}_{j}')

        # Define w
        for j in self.teams:
            if j == 1:
                continue
            self.w[j] = self.solver.IntVar(0, 1, f'w_{j}')

        # Define p
        for i in self.teams:
            self.p[i] = self.solver.IntVar(0, 3 * 2 * (self.n - 1), 'p')

    def define_constraints(self):
        # Constrain x, y, z
        for i in self.teams:
            for j in self.teams:
                if i == j:
                    continue
                self.solver.Add(self.x[i, j] + self.y[i, j] + self.z[i, j] == 1)

        # Link p and x, y, z
        for i in self.teams:
            self.solver.Add(
                self.p[i] == sum(3 * self.x[i, j] + self.y[i, j] + self.y[j, i] + 3 * self.z[j, i]
                                 for j in self.teams if i != j)
            )

        # Constrain w
        self.solver.Add(sum(self.w[j] for j in self.teams if j != 1) == self.k - 1)

        # Link p and w
        for j in self.teams:
            if j == 1:
                continue
            self.solver.Add(self.p[j] <= self.p[1] + self.M * (1 - self.w[j]))
            self.solver.Add(self.p[1] <= self.p[j] + self.M * (self.w[j]))

    def define_objective(self):
        self.solver.Maximize(self.p[1])

    def solve_and_report(self):
        status = self.solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            print('Objective value =', self.solver.Objective().Value())
            print('Problem solved in %f milliseconds' % self.solver.wall_time())
            print('Problem solved in %d iterations' % self.solver.iterations())
            print('Problem solved in %d branch-and-bound nodes' % self.solver.nodes())
        else:
            print('The problem does not have an optimal solution.')

        league = [
            (self.p[i].solution_value(), i,
             sum(self.x[i, j].solution_value() + self.z[j, i].solution_value() for j in self.teams if i != j),
             sum(self.y[i, j].solution_value() + self.y[j, i].solution_value() for j in self.teams if i != j),
             sum(self.z[i, j].solution_value() + self.x[j, i].solution_value() for j in self.teams if i != j))
            for i in self.teams
        ]
        # Build the league table as dataframe
        df = pd.DataFrame(league, columns=['points', 'team', 'wins', 'draws', 'losses'])
        df['team'] = df['team'].astype(int)
        df['wins'] = df['wins'].astype(int)
        df['draws'] = df['draws'].astype(int)
        df['losses'] = df['losses'].astype(int)
        df = df.sort_values(by=['points', 'wins', 'draws', 'losses'], ascending=False)
        df = df.reset_index(drop=True)
        print('\nLeague table:')
        print(df)

    def run(self):
        self.define_variables()
        self.define_constraints()
        self.define_objective()
        self.solve_and_report()


if __name__ == '__main__':
    sports_league = SportsLeagueOpti(18, 3)
    sports_league.run()
