from ortools.linear_solver import pywraplp
import networkx as nx
import matplotlib.pyplot as plt
import itertools


class Steiner:
    def __init__(self):
        self.solver = pywraplp.Solver.CreateSolver('SCIP')

        # Define connected graph as set of edges and weights and nodes
        self.nodes = [1, 2, 3, 4, 5, 6, 7, 8]
        self.edges = [(1, 2), (1, 5), (4, 5), (2, 3), (2, 4), (3, 4), (4, 6), (5, 6), (5, 7), (6, 7), (6, 8), (7, 8)]

        self.weights = {e: 1 for e in self.edges}

        # Define terminals
        self.terminals = [1, 3, 7, 8]
        self.steiner = list(set(self.nodes) - set(self.terminals))
        self.r = self.terminals[0]

        # Define neighborhood
        self.neighborhood = {n: [e for e in self.edges if n in e] for n in self.nodes}

        # Define variables
        self.x = {}
        self.y = {}
        for e in self.edges:
            self.x[e] = self.solver.IntVar(0, 1, f'x_{e}')
            i, j = e
            self.y[i, j] = self.solver.IntVar(0, 1, f'y_{i}_{j}')
            self.y[j, i] = self.solver.IntVar(0, 1, f'y_{j}_{i}')
        self.dir_edges = self.edges + [(j, i) for i, j in self.edges]

    @staticmethod
    def subsets(S):
        return itertools.chain.from_iterable(itertools.combinations(S, r) for r in range(len(S)+1))

    def define_constraints(self):
        # self.solver.Add(sum(self.x[e] for e in self.edges) >= len(self.terminals) - 1)
        for e in self.edges:
            i, j = e
            self.solver.Add(self.y[i, j] + self.y[j, i] == self.x[e])
        # for j in self.terminals:
        #     if j == self.r:
        #         continue
        #     self.solver.Add(sum(self.y[i, j] for i in self.nodes if (i, j) in self.dir_edges) == 1)
        for S in self.subsets(self.nodes):
            if self.r in S and len(S) < len(self.nodes) and set(self.terminals) - set(S):
                # print(f'Subset: {S}')
                # print(f'Constraint: {sum(self.y[i, j] for i in S for j in set(self.nodes) \
                # - set(S) if (i, j) in self.dir_edges) >= 1}')
                self.solver.Add(sum(self.y[i, j] for (i, j) in self.dir_edges if i in S and j not in S) >= 1)

    def define_objective(self):
        self.solver.Minimize(sum(self.weights[e] * self.x[e] for e in self.edges))

    def validate_solution(self, draw=True):
        # Build networkx graph from solution
        G = nx.Graph()
        for e in self.edges:
            if self.x[e].solution_value() == 1:
                # Add nodes to graph
                G.add_node(e[0])
                G.add_node(e[1])
                # Add edge to graph
                G.add_edge(*e)

        if draw:
            # Draw graph
            nx.draw(G, with_labels=True)
            # Show graph
            plt.show()
        # Check if graph is connected
        print(f'Is connected: {nx.is_connected(G)}')
        # Check if terminals are in graph
        print(f'Are terminals in graph: {all(t in G.nodes for t in self.terminals)}')
        # Check if terminals are connected
        print(f'Are terminals connected: {all(nx.has_path(G, self.r, self.t) for self.t in self.terminals if self.t != self.r)}')

        # Build original graph:
        print('-' * 20)
        G0 = nx.Graph()
        G0.add_nodes_from(self.nodes)
        G0.add_edges_from(self.edges)
        print('Steiner tree according to NX: ')
        print(nx.approximation.steiner_tree(G0, self.terminals))
        print('Our solution: ')
        print(G)
        if draw:
            # Draw graph
            nx.draw(nx.approximation.steiner_tree(G0, self.terminals), with_labels=True)
            # Show graph
            plt.show()

    def run(self, draw=True):
        self.define_objective()
        self.define_constraints()
        # print(self.solver.ExportModelAsLpFormat(False))
        self.print_solution()
        # If solution is feasible, validate it
        if self.solver.VerifySolution(1e-7, True):
            self.validate_solution(draw=draw)

    def print_solution(self):
        print('Number of variables =', self.solver.NumVariables())
        print('Number of constraints =', self.solver.NumConstraints())
        print('Solving the problem ...')
        status = self.solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            print('Objective value =', self.solver.Objective().Value())
            for e in self.edges:
                if self.x[e].solution_value() == 1:
                    print(f'x_{e} = {self.x[e].solution_value()}')
                    i, j = e
                    print(f'y_{i}_{j} = {self.y[i, j].solution_value()}')
                    print(f'y_{j}_{i} = {self.y[j, i].solution_value()}')


if __name__ == '__main__':
    steiner = Steiner()
    steiner.run(draw=False)
