from ortools.linear_solver import pywraplp
import itertools, random


class BinPacking:
    def __init__(self, bin_size, items, formulation, relax=False):
        self.bin_size = bin_size
        self.items = items
        self.model = pywraplp.Solver.CreateSolver('SCIP')
        self.formulation = formulation
        self.relax = relax

        self.x = {}
        self.y = {}

        self.A = None

    def define_variables(self):
        if self.formulation == 'poly':
            for i in range(len(items)):
                if self.relax:
                    self.y[i] = self.model.NumVar(0, 1, f'y_{i}')
                else:
                    self.y[i] = self.model.IntVar(0, 1, f'y_{i}')

            for i in range(len(items)):
                for j in range(len(items)):
                    if self.relax:
                        self.x[i, j] = self.model.NumVar(0, 1, f'x_{i}_{j}')
                    else:
                        self.x[i, j] = self.model.IntVar(0, 1, f'x_{i}_{j}')
        elif self.formulation == 'exp':
            # Calling A[w][j] = 1 if item j selected in bin w
            self.A = [Aw for Aw in itertools.product([0, 1], repeat=len(items))
                      if sum(Aw[j] * self.items[j] for j in range(len(items))) <= self.bin_size]
            for w in range(len(self.A)):
                if self.relax:
                    self.x[w] = self.model.NumVar(0, 1, f'x_{w}')
                else:
                    self.x[w] = self.model.IntVar(0, 1, f'x_{w}')

        else:
            raise ValueError('Formulation not implemented')

    def define_constraints(self):
        if self.formulation == 'poly':
            for j in range(len(items)):
                self.model.Add(sum(self.x[i, j] for i in range(len(items))) == 1)
            for i in range(len(items)):
                self.model.Add(
                    sum(self.x[i, j] * self.items[j] for j in range(len(items))) <= self.bin_size * self.y[i])
        elif self.formulation == 'exp':
            for j in range(len(items)):
                self.model.Add(sum(self.x[w] * self.A[w][j] for w in range(len(self.A))) == 1)

    def define_objective(self):
        if self.formulation == 'poly':
            self.model.Minimize(sum(self.y[i] for i in range(len(items))))
        elif self.formulation == 'exp':
            self.model.Minimize(sum(self.x[w] for w in range(len(self.A))))
        else:
            raise ValueError('Formulation not implemented')

    def print_solution(self):
        print('\nNumber of bins:', self.model.Objective().Value())
        if self.formulation == 'poly':
            for i in range(len(items)):
                for j in range(len(items)):
                    if self.x[i, j].solution_value() > 0:
                        print(f'x_{i}_{j} = {self.x[i, j].solution_value()}')
                if self.y[i].solution_value() > 0:
                    print(f'y_{i} = {self.y[i].solution_value()}')

        elif self.formulation == 'exp':
            for w in range(len(self.A)):
                if self.x[w].solution_value() > 0:
                    print(f'Bin {w}:', [j for j in range(len(items)) if self.A[w][j] == 1])
        else:
            raise ValueError('Formulation not implemented')

    def run(self):
        self.define_variables()
        self.define_constraints()
        self.define_objective()
        self.model.Solve()
        self.print_solution()
        return self.model.Objective().Value()


def run_combinations(bin_size, items):
    return [BinPacking(bin_size, items, 'poly').run(), BinPacking(bin_size, items, 'poly', relax=True).run(),
            BinPacking(bin_size, items, 'exp').run(), BinPacking(bin_size, items, 'exp', relax=True).run()]


if __name__ == '__main__':
    bin_size = 4
    items = [1, 1]
    d = run_combinations(bin_size, items)
    print(d)

    # for i in range(1000):
    #     bin_size = random.randint(1, 15)
    #     items = [random.randint(1, bin_size) for _ in range(random.randint(1, 10))]
    #     d = run_combinations(bin_size, items)
    #     if round(d[0]) != round(d[2]):
    #         print('formulations do not agree on optimum')
    #     if round(d[1], 4) < round(d[3], 4):
    #         print('found!')
    #         print(bin_size, items)
    #         print(d)
