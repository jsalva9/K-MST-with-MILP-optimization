from dataclasses import dataclass
from typing import Set, Tuple, Dict
import pandas as pd
import networkx as nx


@dataclass
class Instance:
    """
    Instance of the k-MST problem.
    """
    name: str
    n: int
    m: int
    Vext: Set[int]
    Eext: Set[Tuple[int, int]]
    weights: Dict[Tuple[int, int], int]
    V: Set[int] = None
    E: Set[Tuple[int, int]] = None
    A: Set[Tuple[int, int]] = None
    Aext: Set[Tuple[int, int]] = None
    k: int = None

    def __post_init__(self):
        """
        Define the original graph from the extended graph.
        """
        self.V = self.Vext - {0}
        self.E = self.Eext - {(0, v) for v in self.V} - {(v, 0) for v in self.V}
        self.A = {(u, v) for (u, v) in self.E}.union({(v, u) for (u, v) in self.E})
        self.Aext = {(u, v) for (u, v) in self.Eext}.union({(v, u) for (u, v) in self.Eext})
        self.weights = {**self.weights, **{(v, u): w for (u, v), w in self.weights.items()}}

    def find_tree(self):
        """
        Find a tree with k vertices in the graph. We consider as many trees as the number of vertices in the graph.
        For each vertex, we run DFS to find a connected k-subset of vertices containing it. Then, we find the minimum
         spanning tree of the subset using networkx. Finally return the minimum spanning tree with the smallest cost.

        Returns:
            tree (nx.Graph): k-subtree
        """
        g = nx.Graph()
        g.add_nodes_from(self.V)
        g.add_weighted_edges_from([(e[0], e[1], self.weights[e]) for e in self.E])

        best_tree, best_tree_cost = None, float('inf')
        for i in self.V:
            dfs_nodes = list(nx.dfs_preorder_nodes(g, source=i))[:self.k]
            subgraph = g.subgraph(nodes=dfs_nodes)
            tree = nx.minimum_spanning_tree(subgraph)
            if tree.size(weight='weight') < best_tree_cost:
                best_tree = tree
                best_tree_cost = tree.size(weight='weight')

        print(f'Weight of the tree: {best_tree.size(weight="weight")}')
        assert nx.is_tree(best_tree), 'Found something that is not a tree'
        assert best_tree.number_of_nodes() == self.k, 'Found a tree with the wrong number of nodes'
        assert best_tree.number_of_edges() == self.k - 1, 'Found a tree with the wrong number of edges'

        return best_tree


OPTIMAL_SOLUTIONS = {
    '1_0': 46, '1_1': 477, '2_0': 373, '2_1': 1390,
    '3_0': 725, '3_1': 3074, '4_0': 909, '4_1': 3292,
    '5_0': 1235, '5_1': 4898, '6_0': 2068, '6_1': 6705,
    '7_0': 1335, '7_1': 4534, '8_0': 1620, '8_1': 5787,
    '9_0': 2289, '9_1': 7595, '10_0': 4182, '10_1': 14991,
}


def compute_results_statistics():
    """Compute statistics for results.csv"""

    df = pd.read_csv('../kmst_output/results.csv')
    df.loc[(df.solve_time >= 3599) | df.opt_gap.isna(), 'solve_time'] = float('nan')
    df_wide = pd.pivot(df, index='instance', values='solve_time', columns=['solver', 'formulation', 'tighten'])
    # Compute relative difference tighten vs. no tighten
    df_wide['relative_difference_tighten_gurobi_mcf'] = (df_wide['gurobi']['MCF'][True] - df_wide['gurobi']['MCF'][
        False]) / df_wide['gurobi']['MCF'][False]
    df_wide['relative_difference_tighten_gurobi_mtz'] = (df_wide['gurobi']['MTZ'][True] - df_wide['gurobi']['MTZ'][
        False]) / df_wide['gurobi']['MTZ'][False]
    df_wide['relative_difference_tighten_gurobi_scf'] = (df_wide['gurobi']['SCF'][True] - df_wide['gurobi']['SCF'][
        False]) / df_wide['gurobi']['SCF'][False]
    # Compute relative difference gurobi vs. ortools
    df_wide['relative_difference_solver_mcf'] = (df_wide['gurobi']['MCF'][True] - df_wide['ortools']['MCF'][True]) / \
                                                df_wide['ortools']['MCF'][True]
    df_wide['relative_difference_solver_mtz'] = (df_wide['gurobi']['MTZ'][True] - df_wide['ortools']['MTZ'][True]) / \
                                                df_wide['ortools']['MTZ'][True]
    df_wide['relative_difference_solver_scf'] = (df_wide['gurobi']['SCF'][True] - df_wide['ortools']['SCF'][True]) / \
                                                df_wide['ortools']['SCF'][True]
    # Compute relative difference formulations
    df_wide['relative_difference_formulation_gurobi_mcf'] = (df_wide['gurobi']['MCF'][True] - df_wide['gurobi']['MTZ'][
        True]) / df_wide['gurobi']['MTZ'][True]
    df_wide['relative_difference_formulation_gurobi_scf'] = (df_wide['gurobi']['SCF'][True] - df_wide['gurobi']['MTZ'][
        True]) / df_wide['gurobi']['MTZ'][True]
    df_wide['relative_difference_formulation_ortools_mcf'] = (df_wide['ortools']['MCF'][True] -
                                                              df_wide['ortools']['MTZ'][True]) / \
                                                             df_wide['ortools']['MTZ'][True]
    df_wide['relative_difference_formulation_ortools_scf'] = (df_wide['ortools']['SCF'][True] -
                                                              df_wide['ortools']['MTZ'][True]) / \
                                                             df_wide['ortools']['MTZ'][True]

    df_wide.to_csv(f'../kmst_output/results_agg.csv', index=True)
