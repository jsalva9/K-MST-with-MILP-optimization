import networkx as nx
from gurobipy import GRB
import gurobipy as gp


def find_cycle(x: dict):
    """
    Find a cycle in the graph defined by the current solution (either integral or fractional).

    Args:
        x (dict): map edge to its value in the current solution

    Returns:
        cycle (list): list of edges in the cycle. If no cycle is found, return None.
    """
    # Find a cycle in a list of edges using networkx
    G = nx.Graph()
    G.add_edges_from([k for k, v in x.items() if v > 0.0001])
    try:
        cycle = nx.find_cycle(G)
        cycle = [e if e[0] < e[1] else (e[1], e[0]) for e in cycle]
        return cycle
    except nx.NetworkXNoCycle:
        return None


def find_min_cut(x: dict, z: dict):
    """
    Find a minimum cut in the graph defined by the current solution (either integral or fractional).

    Args:
        x (dict): map edge to its value in the current solution
        z (dict): map vertex to its value in the current solution

    Returns:
        outgoing_e (list): list of edges in the minimum cut (direction 0 -> i). If no cut is found, return None.
        i (int): vertex for which the valid inequality has been found. If no cut is found, return None.
    """
    G = nx.DiGraph()
    G.add_weighted_edges_from([(k[0], k[1], v) for k, v in x.items()])
    # Print min and max weight of edges by accessing G attributes
    for i in z.keys():
        flow = nx.maximum_flow_value(G, 0, i, capacity='weight')
        if flow < z[i] - 0.0001:
            val, partition = nx.minimum_cut(G, 0, i, capacity='weight')
            outgoing_e = [(u, v) for u, v in x.keys() if u in partition[0] and v in partition[1]]
            return outgoing_e, i
    return None, None


def cycle_elimination_constraint_fractional(model: gp.Model, where):
    """
    Callback function for fractional solution. If a fractional solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurobi model
        where (int): callback location
    """
    if (where == GRB.Callback.MIPNODE) and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        x = model.cbGetNodeRel(model._x)
    else:
        return
    add_cec(model, x)


def cycle_elimination_constraint(model: gp.Model, where):
    """
    Callback function for integral solution. If an integral solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurobi model
        where (int): callback location
    """
    if where != GRB.Callback.MIPSOL:
        return
    x = model.cbGetSolution(model._x)
    add_cec(model, x)


def cycle_elimination_constraint_both(model: gp.Model, where):
    """
    Callback function for both fractional and integral solution.
    If a solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurobi model
        where (int): callback location
    """
    cycle_elimination_constraint(model, where)
    cycle_elimination_constraint_fractional(model, where)


def add_cec(model: gp.Model, x):
    """
    Search for violated cycle elimination constraint in the current solution.
    If found, add constraint to the model in the form of a lazy constraint.

    Args:
        model (gp.Model): Gurobi model
        x (dict): map edge to its value in the current solution
    """
    # Find the shortest cycle in the selected edges
    cycle = find_cycle(x)
    if cycle is not None:
        # Add lazy constraint to the model
        model.cbLazy(gp.quicksum(model._x[e] for e in cycle) <= len(cycle) - 1)


def directed_cutset_constraint(model: gp.Model, where):
    """
    Callback function for integral solution. If an integral solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurobi model
        where (int): callback location
    """
    if where != GRB.Callback.MIPSOL:
        return
    x = model.cbGetSolution(model._x)
    z = model.cbGetSolution(model._z)

    add_cutset(model, x, z)


def directed_cutset_constraint_fractional(model: gp.Model, where):
    """
    Callback function for fractional solution. If a fractional solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurobi model
        where (int): callback location
    """
    if (where == GRB.Callback.MIPNODE) and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        x = model.cbGetNodeRel(model._x)
        z = model.cbGetNodeRel(model._z)
    else:
        return
    add_cutset(model, x, z)


def directed_cutset_constraint_both(model: gp.Model, where):
    """
    Callback function for both fractional and integral solution.
    If a solution is found, add a lazy constraint to the model.

    Args:
        model (gp.Model): Gurboi model
        where (int): callback location
    """
    directed_cutset_constraint(model, where)
    directed_cutset_constraint_fractional(model, where)


def add_cutset(model: gp.Model, x: dict, z: dict):
    """
    Search for violated directed cutset constraint in the current solution.
    If found, add constraint to the model in the form of a lazy constraint.

    Args:
        model (gp.Model): Gurobi model
        x (dict): map edge to its value in the current solution
        z (dict): map vertex to its value in the current solution
    """
    outgoing_e, i = find_min_cut(x, z)
    if i is not None:
        # Add lazy constraint to the model
        model.cbLazy(gp.quicksum(model._x[e] for e in outgoing_e) >= model._z[i])
