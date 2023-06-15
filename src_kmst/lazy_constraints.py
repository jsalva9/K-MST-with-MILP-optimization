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
    # Find some vertex i and a partition of the nodes V = W + V \ W such that:
    #   - W contains the root node 0
    #   - W does not contain the vertex i
    #   - outgoing edges e of W satisfy (sum(x[e] for e in outgoing) < z[i])
    # Returns:
    #   outgoing_e (list): list of outgoing edges of the set W
    #   i (node): one such vertex

    # Create graph
    G = nx.DiGraph()

    # Add edges that are in the solution
    G.add_edges_from([k for k, v in x.items() if v > 0.5])

    # For every node i such that y_i = 1, check if there is a path from 0 to i
    for i in [i for i, v in z.items() if v > 0.5]:
        if nx.has_path(G, 0, i):
            continue
        # Found a node i such that there is no path from 0 to i -> there is a partition whose edges values are 0
        Gaux = nx.Graph()
        Gaux.add_weighted_edges_from([(k[0], k[1], v) for k, v in x.items()])

        # Find a 0->i min-cut. We know that must have value = len(outgoing_e) because we added 1 to all edges
        cut_value, partition = nx.minimum_cut(Gaux, 0, i)
        # Finally find all edges that are in the cut
        outgoing_e = [(u, v) for u, v in Gaux.edges if u in partition[0] and v in partition[1]]
        assert sum([x[e] for e in outgoing_e]) == len(outgoing_e), 'CUT VALUE IS NOT ZERO'
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
    outgoing_e, i = find_partition(x, z)
    if i is not None:
        # Add lazy constraint to the model
        model.cbLazy(gp.quicksum(model._x[e] for e in outgoing_e) >= model._z[i])
