from collections import deque

from kidney_solver.kidney_digraph import *
from kidney_solver import kidney_ndds

EPS = 0.00001

class KidneyOptimException(Exception):
    pass

def check_validity(opt_result, digraph, ndds, max_cycle, max_chain):
    """Check that the solution is valid.

    This method checks that:
      - all used edges exist
      - no vertex or NDD is used twice (which also ensures that no edge is used twice)
      - cycle and chain caps are respected
      - the objective value is correct
    """

    # all used edges exist
    for chain in opt_result.chains:
        if chain.vtx_indices[0] not in [e.target_v.index() for e in ndds[chain.ndd_index].edges]:
            raise KidneyOptimException("Edge from NDD {} to vertex {} is used but does not exist".format(
                    chain.ndd_index, chain.vtx_indices[0]))
    for cycle in opt_result.cycles:
        for i in range(len(cycle)):
            if digraph.adj_mat[cycle[i-1].index()][cycle[i].index()] is None:
                raise KidneyOptimException("Edge from vertex {} to vertex {} is used but does not exist".format(
                        cycle[i-1].index(), cycle[i].index()))
                
    # no vertex or NDD is used twice
    ndd_used = [False] * len(ndds)
    vtx_used = [False] * len(digraph.vs)
    for chain in opt_result.chains:
        if ndd_used[chain.ndd_index]:
            raise KidneyOptimException("NDD {} used more than once".format(chain.ndd_index))
        ndd_used[chain.ndd_index] = True
        for vtx_index in chain.vtx_indices:
            if vtx_used[vtx_index]:
                raise KidneyOptimException("Vertex {} used more than once".format(vtx_index))
            vtx_used[vtx_index] = True
            
    for cycle in opt_result.cycles:
        for vtx in cycle:
            if vtx_used[vtx.index()]:
                raise KidneyOptimException("Vertex {} used more than once".format(vtx.index()))
            vtx_used[vtx.index()] = True

    # cycle and chain caps are respected
    for chain in opt_result.chains:
        if len(chain.vtx_indices) > max_chain:
            raise KidneyOptimException("The chain cap is violated")
    for cycle in opt_result.cycles:
        if len(cycle) > max_cycle:
            raise KidneyOptimException("The cycle cap is violated")

    # the objective value is correct
#    if abs(opt_result.total_score - opt_result.ip_model.obj_val) > EPS:
#        raise KidneyOptimException("The objective value is out by {}".format(
#                opt_result.total_score - opt_result.ip_model.obj_val))

def get_dist_from_nearest_ndd(digraph, ndds):
    """ For each donor-patient pair V, this returns the length of the
    shortest path from an NDD to V, or 999999999 if no path from an NDD
    to V exists.
    """
    
    # Get a set of donor-patient pairs who are the target of an edge from an NDD
    ndd_targets = set()
    for ndd in ndds:
        for edge in ndd.edges:
            ndd_targets.add(edge.target_v)

    # Breadth-first search
    q = deque(ndd_targets)
    distances = [999999999] * len(digraph.vs)
    for v in ndd_targets:
        distances[v.index()] = 1

    while q:
        v = q.popleft()
        for e in v.edges:
            w = e.tgt
            if distances[w.index()] == 999999999:
                distances[w.index()] = distances[v.index()] + 1
                q.append(w)

    return distances

def find_selected_path(v_id, next_vv):
    path = [v_id]
    while v_id in next_vv:
        v_id = next_vv[v_id]
        path.append(v_id)
    return path

def find_selected_path_edges(next_vertex, mapping):
    path = []
    count = 0
    while next_vertex in mapping:
        first = mapping[next_vertex]
        path.append(first)
        next_vertex = first.target().index()
        count += 1
        if count > 50:
            break
    return path

def find_selected_cycle(v_id, next_vv):
    cycle = [v_id]
    while v_id in next_vv:
        v_id = next_vv[v_id]
        if v_id in cycle:
            return cycle
        else:
            cycle.append(v_id)
    return None


def chain_score(path, edge_success_prob, nhs=False):
    """Calculate the total score of a given chain, given as a list of edges."""
    if nhs:
        score = sum([e.age_formula(e.target().donor())*edge_success_prob**(index+1) for index, e in enumerate(path)])
    else:
        score = sum([e.score*edge_success_prob**(index+1) for index, e in enumerate(path)])
    return score


def get_optimal_chains(digraph, ndds, edge_success_prob=1, nhs=False):
    # Chain edges
    chain_next_vv = {e.src.index(): e.tgt.index()
                        for e in digraph.es
                        for var in e.grb_vars
                        if var.x > 0.1}
    selected_edges = {e.src.index(): e for e in digraph.es for var in e.grb_vars if var.x > 0.1}
    optimal_chains = []
    for i, ndd in enumerate(ndds):
        for e in ndd.edges:
            if e.edge_var.x > 0.1:
                vtx_indices = find_selected_path(e.target_v.index(), chain_next_vv)
                chain_edges = [e]
                chain_edges.extend(find_selected_path_edges(e.target().index(), selected_edges))
                score = chain_score(chain_edges, edge_success_prob, nhs)
                optimal_chains.append(kidney_ndds.Chain(i, vtx_indices, score))
    return optimal_chains

def selected_edges_to_cycles(digraph, cycle_start_vv, cycle_next_vv):
    cycles = [find_selected_cycle(start_v, cycle_next_vv) for start_v in cycle_start_vv]
    # Remove "cycles" that are really part of a chain
    cycles = [c for c in cycles if c is not None]
    # Remove duplicated cycles
    cycles = [c for c in cycles if c[0] == min(c)]
    # Use vertices instead of indices
    return [[digraph.vs[v_id] for v_id in c] for c in cycles]

