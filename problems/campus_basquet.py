from ortools.linear_solver import pywraplp


class CampusOpti:
    def __init__(self, T):
        self.T = T
        self.team_letters = ['A', 'B', 'C', 'D']
        self.team_index = [1, 2]
        self.teams = [(i, j) for i in self.team_letters for j in self.team_index]
        self.times = [i for i in range(1, T + 1)]
        self.stations = [1, 2, 3, 4, 5, 6]

        self.solver = pywraplp.Solver.CreateSolver('SAT')
        self.y = {}
        self.z = None

    def define_variables(self):
        # Define y, z
        for (i1, j1) in self.teams:
            for (i2, j2) in self.teams:
                for s in self.stations:
                    for t in self.times:
                        self.y[i1, j1, i2, j2, s, t] = self.solver.IntVar(0, 1, f'y_{i1}{j1}_{i2}{j2}_{s}_{t}')

        self.z = self.solver.IntVar(0, self.T, 'z')

    def define_constraints(self):

        # Max one game per court per time
        for s in self.stations:
            for t in self.times:
                self.solver.Add(sum(self.y[i1, j1, i2, j2, s, t] for (i1, j1) in self.teams for (i2, j2) in self.teams) <= 2)

        for (i1, j1) in self.teams:
            for s in self.stations:
                # Team does not play against itself
                for t in self.times:
                    self.solver.Add(self.y[i1, j1, i1, j1, s, t] == 0)

                # Each team plays once in each station
                self.solver.Add(sum(self.y[i1, j1, i2, j2, s, t] for (i2, j2) in self.teams for t in self.times) == 1)

            # Each team plays max once in each time slot
            for t in self.times:
                self.solver.Add(sum(self.y[i1, j1, i2, j2, s, t] for (i2, j2) in self.teams for s in self.stations) <= 1)

            for (i2, j2) in self.teams:
                # Symmetry
                for s in self.stations:
                    for t in self.times:
                        self.solver.Add(self.y[i1, j1, i2, j2, s, t] == self.y[i2, j2, i1, j1, s, t])

                # No games between teams with same team letter, one game between teams with different team letter
                if i1 == i2:
                    for s in self.stations:
                        for t in self.times:
                            self.solver.Add(self.y[i1, j1, i2, j2, s, t] == 0)
                else:
                    self.solver.Add(sum(self.y[i1, j1, i2, j2, s, t] for s in self.stations for t in self.times) == 1)

                # z is the time of the last game played
                for s in self.stations:
                    for t in self.times:
                        self.solver.Add(self.z >= self.y[i1, j1, i2, j2, s, t] * t)

    def define_objective(self):
        self.solver.Minimize(self.z)

    def solve_and_report(self):
        status = self.solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            print('Objective value =', self.solver.Objective().Value())
            print('Problem solved in %f milliseconds' % self.solver.wall_time())
            print('Problem solved in %d iterations' % self.solver.iterations())
            print('Problem solved in %d branch-and-bound nodes' % self.solver.nodes())
        else:
            print('The problem does not have an optimal solution.')

        # For each time step until self.z, print the games played
        for t in range(1, int(self.solver.Objective().Value()) + 1):
            print(f'\nTime {t}:')
            for (i1, j1) in self.teams:
                for (i2, j2) in self.teams:
                    for s in self.stations:
                        if i1 < i2 and self.y[i1, j1, i2, j2, s, t].solution_value() == 1:
                            print(f'Court {s}: {i1}{j1} vs {i2}{j2}')

    def run(self):
        self.define_variables()
        self.define_constraints()
        self.define_objective()
        self.solve_and_report()


if __name__ == '__main__':
    sports_league = CampusOpti(10)
    sports_league.run()
