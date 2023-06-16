import networkx as nx
from gurobipy import GRB
import gurobipy as gp


def find_cycle(x):
    # Find a cycle in a list of edges using networkx
    G = nx.Graph()
    G.add_edges_from([k for k, v in x.items() if v > 0.5])
    try:
        cycle = nx.find_cycle(G)
        cycle = [e if e[0] < e[1] else (e[1], e[0]) for e in cycle]
        return cycle
    except nx.NetworkXNoCycle:
        return None


def find_partition(x, z):
    # Find some vertex i and a partition of the nodes V = W_0 + W_1 such that:
    #   - W_0 contains the root node 0
    #   - W_1 contains node i
    #   - outgoing edges W_0 -> W_1 satisfy (sum(x[e] for e in outgoing) == 0 < z[i] == 1)

    # Create graph with edges in the current solution
    G = nx.DiGraph([k for k, v in x.items() if v > 0.5])

    # For every node i such that y_i = 1, check if there is a path from 0 to i
    for i in [i for i, v in z.items() if v > 0.5]:
        assert 0 in G.nodes, 'Root node 0 is not in the graph'
        assert i in G.nodes, 'Node i is not in the graph'
        if nx.has_path(G, 0, i):
            continue

        # DFS from 0
        nodes = set(nx.dfs_postorder_nodes(G, 0))
        assert i not in nodes, 'Node i is reachable from 0!'
        W_1 = [v for v in z.keys() if v not in nodes]
        W_0 = [v for v in z.keys() if v not in W_1] + [0]

        outgoing_e = [(u, v) for u, v in x.keys() if u in W_0 and v in W_1]
        assert round(sum([x[e] for e in outgoing_e]), 4) == 0, 'Cut value is not zero'

        return outgoing_e, i
    return None, None


def find_min_cut(x, z):
    G = nx.DiGraph()
    G.add_weighted_edges_from([(k[0], k[1], v) for k, v in x.items()])
    # Print min and max weight of edges by accessing G attributes
    for i in z.keys():
        flow = nx.maximum_flow_value(G, 0, i, capacity='weight')
        if flow < z[i] - 0.0001:
            val, partition = nx.minimum_cut(G, 0, i, capacity='weight')
            outgoing_e = [(u, v) for u, v in x.keys() if u in partition[0] and v in partition[1]]
            assert round(sum([x[e] for e in outgoing_e]), 4) == 0, 'Cut value is not zero'
            assert round(z[i], 4) == 1, 'z value is not one'
            return outgoing_e, i
    return None, None


def cycle_elimination_constraint_fractional(model: gp.Model, where):
    if (where == GRB.Callback.MIPNODE) and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        x = model.cbGetNodeRel(model._x)
    else:
        return
    add_cec(model, x)


def cycle_elimination_constraint(model: gp.Model, where):
    if where != GRB.Callback.MIPSOL:
        return
    x = model.cbGetSolution(model._x)
    add_cec(model, x)


def cycle_elimination_constraint_both(model: gp.Model, where):
    cycle_elimination_constraint(model, where)
    cycle_elimination_constraint_fractional(model, where)


def add_cec(model: gp.Model, x):
    # Find the shortest cycle in the selected edges
    cycle = find_cycle(x)
    if cycle is not None:
        # Add lazy constraint to the model
        model.cbLazy(gp.quicksum(model._x[e] for e in cycle) <= len(cycle) - 1)


def directed_cutset_constraint(model: gp.Model, where):
    if where != GRB.Callback.MIPSOL:
        return
    x = model.cbGetSolution(model._x)
    z = model.cbGetSolution(model._z)

    add_cutset(model, x, z)


def directed_cutset_constraint_fractional(model: gp.Model, where):
    if (where == GRB.Callback.MIPNODE) and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        x = model.cbGetNodeRel(model._x)
        z = model.cbGetNodeRel(model._z)
    else:
        return
    add_cutset(model, x, z)


def directed_cutset_constraint_both(model: gp.Model, where):
    directed_cutset_constraint(model, where)
    directed_cutset_constraint_fractional(model, where)


def add_cutset(model: gp.Model, x, z):
    outgoing_e, i = find_min_cut(x, z)
    if i is not None:
        # Add lazy constraint to the model
        model.cbLazy(gp.quicksum(model._x[e] for e in outgoing_e) >= model._z[i])
