from src_kmst.utils import Instance, OPTIMAL_SOLUTIONS
from src_kmst.lazy_constraints import cycle_elimination_constraint_fractional, cycle_elimination_constraint, \
    cycle_elimination_constraint_both, directed_cutset_constraint_fractional, directed_cutset_constraint, \
    directed_cutset_constraint_both
from ortools.linear_solver import pywraplp
import pandas as pd
import networkx as nx
from math import ceil
from copy import copy
import time
import gurobipy as gp
from gurobipy import GRB


class KMST:
    def __init__(self, config_dict: dict):
        """
        Initialize KMST class with config_dict. Define solver, instances, variables and results.

        Args:
            config_dict:
        """
        self.solver = None
        self.instance_names = config_dict['instances']
        self.formulation = config_dict['formulation']
        self.solver_name = config_dict['solver']
        self.validate_solution = config_dict['validate_solution']
        self.time_limit = config_dict['time_limit']
        self.hint_solution = config_dict['hint_solution']
        self.tighten = config_dict['tighten']
        self.data_path = config_dict['data_path']
        self.cuts = config_dict['cuts']
        self.exp_form = ['CEC', 'DCC']

        self.start_time = time.time()

        self.instances = {}

        self.x, self.u, self.f, self.z = None, None, None, None  # Variables

        self.results = pd.DataFrame()  # Dataframe with results

    def define_model(self):
        """
        Create the correspondent solver instance according to self.solver_name.
        """
        solvers_map = {
            'ortools': 'SCIP',
            'gurobi': 'GUROBI',
            'cplex': 'CPLEX'
        }
        if self.formulation in self.exp_form:
            assert self.solver_name == 'gurobi', \
                f'Only gurobi solver is available for the exponential formulations {self.exp_form}.'
            self.solver = gp.Model('KMST')
            self.solver.setParam('TimeLimit', self.time_limit)
            self.solver.Params.LazyConstraints = 1
        else:
            self.solver = pywraplp.Solver.CreateSolver(solvers_map[self.solver_name])
            self.solver.set_time_limit(self.time_limit * 1000)

    def load_instances(self):
        """
        Load the graph instances stored in data_path and for each of them, create two problem instances with k = m/5
        and k = m/2.
        """
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
        """
        Load the graph instance stored in data_path with name instance_name. Reads the file and creates an Instance
        object.

        Args:
            instance_name: name of the instance. Example: '01', '02', '03', etc.

        Returns:
            Instance object with the graph instance.
        """
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
        """
        Define the variables according to the formulation. For each formulation, the variables are defined in a
        different method.

        Args:
            instance: Instance object with the problem instance.
        """
        self.x, self.u, self.f, self.z = {}, {}, {}, {}
        if self.formulation == 'MTZ':
            self.define_variables_mtz(instance)
        elif self.formulation == 'SCF':
            self.define_variables_scf(instance)
        elif self.formulation == 'MCF':
            self.define_variables_mcf(instance)
        elif self.formulation == 'CEC':
            self.define_variables_cec(instance)
        elif self.formulation == 'DCC':
            self.define_variables_dcc(instance)
        else:
            raise ValueError('Formulation not recognized')

    def define_constraints(self, instance: Instance):
        """
        Define the constraints according to the formulation. For each formulation, the constraints are defined in a
        different method.

        Args:
            instance: Instance object with the problem instance.
        """
        if self.formulation == 'MTZ':
            self.define_constraints_mtz(instance)
        elif self.formulation == 'SCF':
            self.define_constraints_scf(instance)
        elif self.formulation == 'MCF':
            self.define_constraints_mcf(instance)
        elif self.formulation == 'CEC':
            self.define_constraints_cec(instance)
        elif self.formulation == 'DCC':
            self.define_constraints_dcc(instance)
        else:
            raise ValueError('Formulation not recognized')

    def define_variables_cec(self, instance: Instance):
        """
        Define the variables for the CEC formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        for e in instance.A:
            self.x[e] = self.solver.addVar(lb=0, ub=1, obj=instance.weights[e], vtype=GRB.BINARY, name='x_{e}')
        for i in instance.V:
            self.z[i] = self.solver.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY, name='z_{i}')
        self.solver._x = self.x
        self.solver._z = self.z

    def define_variables_dcc(self, instance: Instance):
        """
        Define the variables for the DCC formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        for e in instance.Aext:
            self.x[e] = self.solver.addVar(lb=0, ub=1, obj=instance.weights[e], vtype=GRB.BINARY, name='x_{e}')
        for i in instance.V:
            self.z[i] = self.solver.addVar(lb=0, ub=1, obj=0, vtype=GRB.BINARY, name='z_{i}')
        self.solver._x = self.x
        self.solver._z = self.z

    def define_constraints_cec(self, instance: Instance):
        """
        Define the constraints for the CEC formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        self.solver.addConstr(gp.quicksum(self.x[e] for e in instance.E) == instance.k - 1, name='total_edges')
        self.solver.addConstr(gp.quicksum(self.z[i] for i in instance.V) == instance.k, name='total_nodes')

        for i in instance.V:
            self.solver.addConstr(gp.quicksum(self.x[e] for e in instance.E if e[0] == i or e[1] == i) <= (instance.k - 1) * self.z[i], name=f'node_1_{i}')
            self.solver.addConstr(self.z[i] <= gp.quicksum(self.x[e] for e in instance.E if e[0] == i or e[1] == i), name=f'node_2_{i}')

    def define_constraints_dcc(self, instance: Instance):
        """
        Define the constraints for the DCC formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        self.solver.addConstr(gp.quicksum(self.x[e] for e in instance.A) == instance.k - 1, name='total_edges')
        self.solver.addConstr(gp.quicksum(self.z[i] for i in instance.V) == instance.k, name='total_nodes')
        self.solver.addConstr(gp.quicksum(self.x[(0, i)] for i in instance.V) == 1, name='root_out')

        for i in instance.V:
            self.solver.addConstr(gp.quicksum(self.x[e] for e in instance.A if i in e) <= (instance.k - 1) * self.z[i], name=f'node_1_{i}')
            self.solver.addConstr(self.z[i] <= gp.quicksum(self.x[e] for e in instance.A if e[1] == i) +
                                  gp.quicksum(self.x[e] for e in instance.A if e[0] == i), name=f'node_2_{i}')
            self.solver.addConstr(self.x[(i, 0)] == 0, name=f'root_in_{i}')
        for e in instance.Eext:
            if e[0] > e[1]:
                continue
            self.solver.addConstr(self.x[e] + self.x[e[1], e[0]] <= 1, name=f'edge_{e}')

    def define_variables_mtz(self, instance: Instance):
        """
        Define the variables for the MTZ formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        for (i, j) in instance.A:
            self.x[i, j] = self.solver.IntVar(0, 1, f'x_{i}_{j}')
        for i in instance.V:
            self.u[i] = self.solver.IntVar(0, instance.k, f'u_{i}')
            self.z[i] = self.solver.IntVar(0, 1, f'z_{i}')

    def define_constraints_mtz(self, instance: Instance):
        """
        Define the constraints for the MTZ formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        self.solver.Add(sum(self.x[i, j] + self.x[j, i] for (i, j) in instance.E) == instance.k - 1)
        self.solver.Add(sum(self.z[i] for i in instance.V) == instance.k)
        for (i, j) in instance.A:
            if self.tighten:
                self.solver.Add(self.x[i, j] + self.x[j, i] <= 1)
            self.solver.Add(self.u[i] + self.x[i, j] - self.u[j] <= instance.k * (1 - self.x[i, j]))
            # Can be rewritten with z_i instead of x_ij
        for i in instance.V:
            if self.tighten:
                self.solver.Add(self.z[i] <= self.u[i])
                self.solver.Add(self.u[i] <= instance.k * self.z[i])
                self.solver.Add(self.z[i] <= sum(self.x[i, j] for j in instance.V if (i, j) in instance.A) + sum(
                    self.x[j, i] for j in instance.V if (j, i) in instance.A))

            self.solver.Add((sum(self.x[i, j] for j in instance.V if (i, j) in instance.A) + sum(
                self.x[j, i] for j in instance.V if (j, i) in instance.A)) <= (instance.k - 1) * self.z[i])
            self.solver.Add(sum(self.x[j, i] for j in instance.V if (j, i) in instance.A) <= self.z[i])

    def define_variables_scf(self, instance: Instance):
        """
        Define the variables for the SCF formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        for e in instance.Eext:
            self.x[e] = self.solver.IntVar(0, 1, f'x_{e}')
        for i, j in instance.Aext:
            self.f[i, j] = self.solver.IntVar(0, instance.k, f'f_{i}_{j}')
        for i in instance.V:
            self.z[i] = self.solver.IntVar(0, 1, f'z_{i}')

    def define_constraints_scf(self, instance: Instance):
        """
        Define the constraints for the SCF formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        if self.tighten:
            self.solver.Add(sum(self.x[e] for e in instance.E) == instance.k - 1)
        self.solver.Add(sum(self.z[i] for i in instance.V) == instance.k)
        self.solver.Add(instance.k == sum(self.f[0, i] for i in instance.V))
        for i in instance.V:
            self.solver.Add(self.f[0, i] == instance.k * self.x[0, i])
            if self.tighten:
                self.solver.Add(self.f[i, 0] == 0)
            self.solver.Add(sum(self.f[i, j] for j in instance.Vext if (i, j) in instance.Aext) - sum(
                self.f[j, i] for j in instance.Vext if (j, i) in instance.Aext) == - self.z[i])
        for i, j in instance.E:
            self.solver.Add(self.f[i, j] + self.f[j, i] <= (instance.k - 1) * self.x[i, j])
            if self.tighten:
                self.solver.Add(self.f[i, j] <= (instance.k - 1) * self.x[i, j])
                self.solver.Add(self.f[j, i] <= (instance.k - 1) * self.x[i, j])

    def define_variables_mcf(self, instance: Instance):
        """
        Define the variables for the MCF formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        for e in instance.Eext:
            self.x[e] = self.solver.IntVar(0, 1, f'x_{e}')
        for i, j in instance.Aext:
            for l in instance.V:
                self.f[i, j, l] = self.solver.IntVar(0, 1, f'f_{i}_{j}_{l}')
        for i in instance.V:
            self.z[i] = self.solver.IntVar(0, 1, f'z_{i}')

    def define_constraints_mcf(self, instance: Instance):
        """
        Define the constraints for the MCF formulation.

        Args:
            instance: Instance object with the problem instance.
        """
        if self.tighten:
            self.solver.Add(sum(self.x[e] for e in instance.E) == instance.k - 1)
        self.solver.Add(sum(self.z[i] for i in instance.V) == instance.k)
        self.solver.Add(sum(self.f[0, i, l] for i in instance.V for l in instance.V) == instance.k)
        for i in instance.V:
            self.solver.Add(instance.k * self.x[0, i] == sum(self.f[0, i, l] for l in instance.V))
            if self.tighten:
                for l in instance.V:
                    self.solver.Add(self.f[i, 0, l] == 0)

        for (i, j) in instance.E:
            for l in instance.V:
                self.solver.Add(self.f[i, j, l] + self.f[j, i, l] <= self.x[i, j])
                if self.tighten:
                    self.solver.Add(self.f[i, j, l] <= self.x[i, j])
                    self.solver.Add(self.f[j, i, l] <= self.x[i, j])
        for i in instance.Vext:
            for l in instance.V:
                if i == 0:
                    self.solver.Add(sum(self.f[i, j, l] for j in instance.Vext if (i, j) in instance.Aext) - sum(
                        self.f[j, i, l] for j in instance.Vext if (j, i) in instance.Aext) == self.z[l])
                elif i == l:
                    self.solver.Add(sum(self.f[i, j, l] for j in instance.Vext if (i, j) in instance.Aext) - sum(
                        self.f[j, i, l] for j in instance.Vext if (j, i) in instance.Aext) == - self.z[l])
                else:
                    self.solver.Add(sum(self.f[i, j, l] for j in instance.Vext if (i, j) in instance.Aext) - sum(
                        self.f[j, i, l] for j in instance.Vext if (j, i) in instance.Aext) == 0)

    def define_objective(self, instance: Instance):
        """
        Define the objective function for each formulation.
        Args:
            instance: Instance object with the problem instance.
        """
        if self.formulation == 'MTZ':
            self.solver.Minimize(sum(self.x[i, j] * instance.weights[i, j] for (i, j) in instance.A))
        elif self.formulation in ['SCF', 'MCF']:
            self.solver.Minimize(sum(self.x[e] * instance.weights[e] for e in instance.E))
        elif self.formulation in self.exp_form:
            # Already defined as the coefficients of the variables
            pass
        else:
            raise ValueError('Formulation not recognized')

    def solve(self):
        """
        Solve the problem instance. The solver is called here. Print the status of the solution and the time it took

        Returns:
            4-tuple: (True iff optimally solved, solve_time, status integer code, objective value)
        """
        a = time.time()
        if self.formulation in self.exp_form:
            objective, optimal, status = self.solve_with_cuts()
        else:
            status = self.solver.Solve()
            objective = self.solver.Objective().Value()
            optimal = status == pywraplp.Solver.OPTIMAL

        solve_time = time.time() - a
        if optimal:
            print(f' - solved with objective {objective}')
            print(f' - solved in {solve_time} s')
            return True, solve_time, status, objective
        print(f' - not solved. Status: {status} (solver = {self.solver_name})')
        return False, solve_time, status, objective

    def solve_with_cuts(self):
        """
        Solve method for the formulations with exponential number of variables (CEC, DCC). The solver is called here.
        Redirect the correct method depending on the formulation and the cut specification (fractional, integral, both).

        Returns:
            4-tuple: (True iff optimally solved, solve_time, status integer code, objective value)
        """
        func_dict = {
            ('CEC', 'fractional'): cycle_elimination_constraint_fractional,
            ('CEC', 'integral'): cycle_elimination_constraint,
            ('CEC', 'both'): cycle_elimination_constraint_both,
            ('DCC', 'fractional'): directed_cutset_constraint_fractional,
            ('DCC', 'integral'): directed_cutset_constraint,
            ('DCC', 'both'): directed_cutset_constraint_both
        }
        self.solver.optimize(func_dict[(self.formulation, self.cuts)])
        status = self.solver.status
        objective = self.solver.ObjVal
        optimal = status == GRB.OPTIMAL
        return objective, optimal, status

    def validate(self, instance: Instance, objective: float):
        """
        Validate the solution obtained by the solver. If not self.validate_solution, do nothing.
        Otherwise, check if objective is equal to the theorital objective.
        If self.validate_solution == 'hard', check the solution subgraph structure: nodes, edges, connectivity...

        Args:
            instance: Instance object with the problem instance.
            objective: Objective value obtained by the solver.
        """
        if not self.validate_solution:
            return
        # Check if the optimal value is the same as the one obtained by the solver
        equal = round(objective) == OPTIMAL_SOLUTIONS[instance.name]
        print(f' - solved with optimal value: {equal}')
        if self.validate_solution != 'hard':
            return

        print('\n', '-' * 5, 'Validating solution', '-' * 5)
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
        elif self.formulation in ['CEC']:
            x_vals = self.solver.getAttr('X', self.x)
            for e in instance.E:
                if x_vals[e] == 1:
                    G.add_edge(*e)
        elif self.formulation in ['DCC']:
            x_vals = self.solver.getAttr('X', self.x)
            for e in instance.A:
                if x_vals[e] == 1:
                    G.add_edge(*e)

        # Check if the solution is a tree
        print(f'{G.number_of_nodes()} subgraph nodes: {G.nodes()}')
        print(f'{G.number_of_edges()} subgraph edges: {G.edges()}')
        print(f'Subgraph a tree: {nx.is_tree(G)}')
        print(f'Subgraph is connected: {nx.is_connected(G)}')
        if self.formulation == 'MTZ':
            # Check u values for the solution nodes
            u_values = {v: self.u[v].solution_value() for v in G.nodes()}
            print(f'Subgraph u values: {u_values}')
        if self.formulation == 'SCF':
            # Check f values for the solution edges
            f_values = {e: self.f[e].solution_value() for e in instance.Aext if self.f[e].solution_value()}
            print(f'Subgraph f values: {f_values}')
        if self.formulation == 'MCF':
            # Check f values for the solution edges
            f_values = {(i, j, l): self.f[i, j, l].solution_value() for (i, j) in instance.Aext for l in instance.V
                        if self.f[i, j, l].solution_value()}
            print(f'Subgraph f values: {f_values}')
            node_sums = {(i, l): sum(self.f[i, j, l].solution_value() for j in instance.Vext if (i, j) in instance.Aext) - sum(
                self.f[j, i, l].solution_value() for j in instance.Vext if (j, i) in instance.Aext)
                         for i in instance.Vext for l in instance.V}
            print(f"Subgraph node sums: { {k: v for k, v in node_sums.items() if v != 0} }")

        print('-' * 5, 'End validation', '-' * 5)

    def write_results(self, instance: Instance, solve_time: float, status: int, objective: float):
        """
        Create a dataframe with the results of the execution, and append it to the results list
        (will later be concatenated)

        Args:
            instance: Instance object with the problem instance.
            solve_time: Time it took to solve the instance (seconds).
            status: Integer code of the status of the solution.
            objective: Objective value obtained by the solver.

        """
        if (self.formulation in self.exp_form and status != GRB.OPTIMAL) or \
                (self.formulation not in self.exp_form and status == pywraplp.Solver.NOT_SOLVED):
            opt_gap = float('nan')
        else:
            opt_gap = round((objective - OPTIMAL_SOLUTIONS[instance.name]) / OPTIMAL_SOLUTIONS[instance.name], 4)
        d = {
            'instance': instance.name,
            'n': instance.n,
            'm': instance.m,
            'k': instance.k,
            'formulation': self.formulation,
            'define_hints': self.hint_solution,
            'solver': self.solver_name,
            'objective': objective,
            'time': time.time() - self.start_time,
            'nodes': self.solver.nodes() if self.formulation not in self.exp_form else self.solver.NodeCount,
            'iterations': self.solver.iterations() if self.formulation not in self.exp_form else self.solver.IterCount,
            'time_limit': self.time_limit,
            'best_bound': self.solver.Objective().BestBound() if self.formulation not in self.exp_form else self.solver.ObjBound,
            'status': status,
            'theoretical_optimal': OPTIMAL_SOLUTIONS[instance.name],
            'opt_gap': opt_gap,
            'tighten': self.tighten,
            'cuts': self.cuts,
            'solve_time': solve_time
        }
        if not len(self.results):
            self.results = pd.DataFrame(columns=list(d.keys())).astype({'define_hints': bool, 'tighten': bool})
        run_results = pd.DataFrame(d, index=[0])
        self.results = pd.concat([self.results, run_results], ignore_index=True)

    def define_hints(self, instance: Instance):
        """
        WARNING: Under development. Define hints for the solver to start with a solution.

        Args:
            instance: Instance object with the problem instance.
        """
        if self.formulation not in self.exp_form:
            raise Warning('Hints are only implemented for Gurobi -> no hints will be defined')

        # Find a tree of size k
        tree = instance.find_tree()

        if self.formulation == 'CEC':
            for e in instance.E:
                self.x[e].Start = 1 if e in tree.edges() else 0
            for i in instance.V:
                self.z[i].Start = 1 if i in tree.nodes() else 0

        else: # self.formulation == 'DCC'
            tree.add_edge(0, list(tree.nodes())[0])        # Add edge from root to some node in the tree
            # Traverse the tree from node zero and define Start values for the variables
            tree_edges = set(nx.dfs_edges(tree, source=0))
            for e in instance.Aext:
                self.x[e].Start = 1 if e in tree_edges else 0
            for i in instance.V:
                self.z[i].Start = 1 if i in tree.nodes() else 0

    def run(self):
        """
        Load the instances, defne the model, variables, constraints, objective, and solve the problem. Validate the
        solution and write the results to a dataframe. Prints the results to the console.
        """

        self.load_instances()
        for instance_name, instance in self.instances.items():
            print(f'\nRunning instance {instance_name} with {self.formulation}. Solver: {self.solver_name}, tighten: {self.tighten}')
            if self.formulation == 'CEC':
                print(f'Cutting planes: {self.cuts}. Hint solution: {self.hint_solution}')

            self.define_model()
            self.define_variables(instance)
            self.define_constraints(instance)
            self.define_objective(instance)
            print(f' - building model took {time.time() - self.start_time} seconds')
            if self.hint_solution:
                self.define_hints(instance)
            feasible, solve_time, status, objective = self.solve()
            if feasible:
                self.validate(instance, objective)
            self.write_results(instance, solve_time, status, objective)
        print(self.results)

