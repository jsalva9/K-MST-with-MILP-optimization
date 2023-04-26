from math import ceil
from src_kmst.config import Config
from ortools.linear_solver import pywraplp
from src_kmst.utils import Instance, OPTIMAL_SOLUTIONS
import networkx as nx
from copy import copy


class KMST:
    def __init__(self, config: Config):
        self.instance_names = config.config['instances']
        self.formulation = config.config['formulation']
        self.solver_name = config.config['solver']
        self.validate_solution = config.config['validate_solution']
        self.data_path = config.data_path

        self.instances = {}
        solvers_map = {
            'ortools': 'SCIP',
            'gurobi': 'GUROBI'
        }
        self.solver = pywraplp.Solver.CreateSolver(solvers_map[self.solver_name])

        self.x, self.u, self.f, self.z = None, None, None, None  # Variables

    def load_instances(self):
        # For each instance in self.instances, load the instance and create two new instances with k = m/5 and k = m/2
        for instance_name in self.instance_names:
            instance_0 = self.load_instance(instance_name)
            instance_0.k = ceil(instance_0.n / 5)
            instance_0.name = f'{instance_name}_0'
            self.instances[instance_0.name] = instance_0
            instance_1 = copy(instance_0)
            instance_1.k = ceil(instance_1.n / 2)
            instance_1.name = f'{instance_name}_1'
            self.instances[instance_1.name] = instance_1

    def load_instance(self, instance_name: str) -> Instance:
        instance_string = str(instance_name) if len(str(instance_name)) == 2 else f'0{instance_name}'
        f = open(f'{self.data_path}/g{instance_string}.dat', 'r')
        lines = f.readlines()
        f.close()
        n, m = int(lines[0]) - 1, int(lines[1]) - 1
        Vext = {i for i in range(0, n + 1)}
        Eext = set()
        weights = {}
        for line in lines[2:]:
            if len(line) <= 1:
                continue
            _, a, b, w = line.split()
            Eext.add((int(a), int(b)))
            weights[(int(a), int(b))] = int(w)
        return Instance(f'instance_{instance_string}', n, m, Vext, Eext, weights)

    def define_variables(self, instance: Instance):
        self.x, self.u, self.f, self.z = {}, {}, {}, {}
        if self.formulation == 'MTZ':
            self.define_variables_mtz(instance)
        elif self.formulation == 'SCF':
            self.define_variables_scf(instance)
        elif self.formulation == 'MCF':
            self.define_variables_mcf(instance)
        else:
            raise ValueError('Formulation not recognized')

    def define_constraints(self, instance: Instance):
        if self.formulation == 'MTZ':
            self.define_constraints_mtz(instance)
        elif self.formulation == 'SCF':
            self.define_constraints_scf(instance)
        elif self.formulation == 'MCF':
            self.define_constraints_mcf(instance)
        else:
            raise ValueError('Formulation not recognized')

    def define_variables_mtz(self, instance: Instance):
        for (i, j) in instance.A:
            self.x[i, j] = self.solver.IntVar(0, 1, f'x_{i}_{j}')
        for i in instance.V:
            self.u[i] = self.solver.IntVar(0, instance.k, f'u_{i}')
            self.z[i] = self.solver.IntVar(0, 1, f'z_{i}')

    def define_constraints_mtz(self, instance: Instance):
        self.solver.Add(sum(self.x[i, j] + self.x[j, i] for (i, j) in instance.E) == instance.k - 1)
        self.solver.Add(sum(self.z[i] for i in instance.V) == instance.k)
        for (i, j) in instance.A:
            # self.solver.Add(self.x[i, j] + self.x[j, i] <= 1)
            self.solver.Add(self.u[i] + self.x[i, j] - self.u[j] <= instance.k * (1 - self.x[i, j]))
            # Can be rewritten with z_i instead of x_ij
        for i in instance.V:
            # self.solver.Add(self.z[i] <= self.u[i])
            # self.solver.Add(self.u[i] <= instance.k * self.z[i])
            # self.solver.Add(self.z[i] <= sum(self.x[i, j] for j in instance.V if (i, j) in instance.A) + sum(
            #     self.x[j, i] for j in instance.V if (j, i) in instance.A))
            self.solver.Add((sum(self.x[i, j] for j in instance.V if (i, j) in instance.A) + sum(
                self.x[j, i] for j in instance.V if (j, i) in instance.A)) <= (instance.k - 1) * self.z[i])
            # self.solver.Add(sum(self.x[j, i] for j in instance.V if (j, i) in instance.A) <= self.z[i])
            self.solver.Add(sum(self.x[i, j] for j in instance.V if (i, j) in instance.A) <= self.z[i])

    def define_variables_scf(self, instance: Instance):
        for e in instance.Eext:
            self.x[e] = self.solver.IntVar(0, 1, f'x_{e}')
        for i, j in instance.Aext:
            self.f[i, j] = self.solver.IntVar(0, instance.k, f'f_{i}_{j}')
        for i in instance.Vext:
            self.z[i] = self.solver.IntVar(0, 1, f'z_{i}')

    def define_variables_mcf(self, instance: Instance):
        pass

    def define_constraints_scf(self, instance: Instance):
        self.solver.Add(sum(self.x[e] for e in instance.E) == instance.k - 1)
        self.solver.Add(sum(self.z[i] for i in instance.V) == instance.k)
        self.solver.Add(instance.k == sum(self.f[0, i] for i in instance.V))
        for i in instance.V:
            self.solver.Add(sum(self.f[j, i] for j in instance.V if (j, i) in instance.Aext) - sum(
                self.f[i, j] for j in instance.V if (i, j) in instance.Aext) <= self.z[i])
        for i, j in instance.Eext:
            self.solver.Add(self.f[i, j] + self.f[j, i] <= (instance.k - 1) * self.x[i, j]) # Might need to split this into two constraints

    def define_constraints_mcf(self, instance: Instance):
        pass

    def define_objective(self, instance: Instance):
        if self.formulation == 'MTZ':
            self.solver.Minimize(sum(self.x[i, j] * instance.weights[i, j] for (i, j) in instance.A))
        elif self.formulation == 'SCF' or self.formulation == 'MCF':
            self.solver.Minimize(sum(self.x[e] * instance.weights[e] for e in instance.E))
        else:
            raise ValueError('Formulation not recognized')

    def solve(self, instance: Instance):
        status = self.solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            print(f' - solved with objective {self.solver.Objective().Value()}')
            print(f' - solved in {self.solver.WallTime() / 1000} s')
            return True
        print(f' - not solved. Status: {status}')
        return False

    def validate(self, instance: Instance):
        if not self.validate_solution:
            return
        if self.validate_solution == 'easy':
            # Check if the optimal value is the same as the one obtained by the solver
            equal = round(self.solver.Objective().Value()) == OPTIMAL_SOLUTIONS[instance.name]
            print(f' - solved with optimal value: {equal}')
            return

        print('-' * 5, 'Validating solution', '-' * 5)
        print(f'Instance n: {instance.n}')
        print(f'Instance m: {instance.m}')
        print(f'Instance k: {instance.k}')

        # Build the solution in networkx
        G = nx.Graph()
        if self.formulation == 'MTZ':
            for (i, j) in instance.A:
                if self.x[i, j].solution_value() == 1:
                    G.add_edge(i, j)
        elif self.formulation in ['SCF', 'MCF']:
            for e in instance.E:
                if self.x[e].solution_value() == 1:
                    G.add_edge(*e)
        # Check if the solution is a tree
        print(f'Subgraph nodes: {G.nodes()}')
        print(f'Subgraph edges: {G.edges()}')
        print(f'Subgraph a tree: {nx.is_tree(G)}')
        print(f'Subgraph is connected: {nx.is_connected(G)}')
        if self.formulation == 'MTZ':
            # Check u values for the solution nodes
            u_values = {v: self.u[v].solution_value() for v in G.nodes()}
            print(f'Subgraph u values: {u_values}')
        if self.formulation == 'SCF':
            # Check f values for the solution edges
            f_values = {e: self.f[e].solution_value() for e in G.edges()}
            print(f'Subgraph f values: {f_values}')

        print('-' * 5, 'End validation', '-' * 5)

    def run(self):
        self.load_instances()
        for instance_name, instance in self.instances.items():
            print(f'\nRunning instance {instance_name}')
            self.define_variables(instance)
            self.define_constraints(instance)
            self.define_objective(instance)
            feasible = self.solve(instance)
            if feasible:
                self.validate(instance)


if __name__ == '__main__':
    config = Config()
    kmst = KMST(config)
    kmst.run()
