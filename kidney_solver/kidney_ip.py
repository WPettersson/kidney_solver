"""Solving the kidney-exchange problem using the Gurobi IP solver."""

import copy
import logging
import sys
from xml.etree import ElementTree as ET

from kidney_solver.kidney_digraph import cycle_score, cycle_contains
from kidney_solver.kidney_ndds import create_relabelled_ndds, Chain, find_chains
from kidney_solver import kidney_utils

from gurobipy import quicksum, Model, GRB


LOGGER = logging.getLogger(__name__)

###################################################################################################
#                                                                                                 #
#                                  Code used by all formulations                                  #
#                                                                                                 #
###################################################################################################

class OptConfig(object):
    """The inputs (problem instance and parameters) for an optimisation run

    Data members:
        digraph
        ndds
        max_cycle
        max_chain
        verbose: True if and only if Gurobi output should be writtent to screen and log file
        timelimit
        edge_success_prob
        eef_alt_constraints: True if and only if alternative EEF constraints should be used
        lp_file: The name of a .lp file to write, or None if the file should not be written
        relax: True if and only if the LP relaxation should be solved also
        only_size: True if we only maximise size of matching, not weights
    """

    def __init__(self, digraph, ndds=None, max_cycle=3, max_chain=0, verbose=False,
                 timelimit=None, edge_success_prob=1, eef_alt_constraints=False,
                 lp_file=None, relax=False, size=False, details=None):
        self.digraph = digraph
        self.ndds = ndds
        self.max_cycle = max_cycle
        self.max_chain = max_chain
        self.verbose = verbose
        self.timelimit = timelimit
        self.edge_success_prob = edge_success_prob
        self.eef_alt_constraints = eef_alt_constraints
        self.lp_file = lp_file
        self.relax = relax
        self._only_size = size
        self.details = details
        self._constrain_maximal = False
        if size:
            LOGGER.info("Only optimising for size.")

    def only_size(self):
        """Do we only care about the size of the matching found?"""
        return self._only_size

class OptSolution(object):
    """An optimal solution for a kidney-exchange problem instance.
    
    Data members:
        ip_model: The Gurobi Model object
        cycles: A list of cycles in the optimal solution, each represented
            as a list of vertices
        chains: A list of chains in the optimal solution, each represented
            as a Chain object
        total_score: The total score of the solution
    """

    def __init__(self, ip_model, cycles, chains, digraph, ndds=None, edge_success_prob=1, nhs=False):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self._ndds = ndds
        self.total_score = (sum(c.score for c in chains) +
                sum(cycle_score(c, edge_success_prob, nhs) for c in cycles))
        self.edge_success_prob = edge_success_prob
        self.nhs = nhs

    def display(self):
        """Print the optimal cycles and chains to standard output."""
        print(str(self))

    def __str__(self):
        """Get a string representation of the solution."""
        total = 0
        string = ""
        string += "cycle_count: {}\n".format(len(self.cycles))
        string += "chain_count: {}\n".format(len(self.chains))
        string += "cycles:\n"
        for cycle in self.cycles:
            # cs is a list of cycles, with each cycle represented as a list of vertex IDs
            cycle_verts = [v.index() for v in cycle]
            min_index_pos = cycle_verts.index(min(cycle_verts))
            cycle_verts = cycle_verts[min_index_pos:] + cycle_verts[:min_index_pos]
            # Sort the cycles
            cycle_verts.sort()
            score = cycle_score(cycle, self.edge_success_prob, self.nhs)
            string += "\t".join(str(v_id) for v_id in cycle_verts) + "\t%f" % score + "\n"
            total += len(cycle)
        string += "chains:\n"
        for chain in self.chains:
            string += str(chain.ndd_index) + "\t" + "\t".join(str(v) for v in chain.vtx_indices)
            string += "\n"
            total += len(chain)+1
        string += "total number of transplants: %d" % total
        return string

    def as_xml(self, filename):
        """Write an XML representation of this solution to a file.

        :param filename: The file to write to.
        """
        root = ET.Element("data")
        ET.SubElement(root, "algorithm").text = "kidney_solver"
        output = ET.SubElement(root, "output")
        all_cycles = ET.SubElement(output, "all_cycles")
        index, cycle_weight = self._xml_add_cycles(all_cycles)
        index, chain_weight = self._xml_add_chains(all_cycles, index)
        exchange_data = ET.SubElement(output, "exchange_data")
        entry = ET.SubElement(exchange_data, "entry")
        entry.set("weight", "%f" % (cycle_weight + chain_weight))
        two_ways = len([x for x in self.cycles if len(x) == 2])
        three_ways = len([x for x in self.cycles if len(x) == 3])
        entry.set("two_way_exchanges", "%d" % two_ways)
        entry.set("three_way_exchanges", "%d" % three_ways)
        num_transplants = sum([len(x) for x in self.chains])
        num_transplants += 2 * two_ways + 3 * three_ways
        entry.set("total_transplants", "%d" % num_transplants)
        ET.SubElement(entry, "description").text = "kidney_solver using __"
        exchanges = ET.SubElement(entry, "exchanges")
        for count in range(1, index):
            ET.SubElement(exchanges, "cycle").text = "%d" % count
        ET.ElementTree(element=root).write(filename)


    def _xml_add_cycles(self, cycles):
        """Private function to add cycles to an XML element.
        """
        cycle_index = 1
        total_weight = 0
        for cycle in self.cycles:
            cycle_weight = 0
            backarcs = 0
            cycle_elt = ET.SubElement(cycles, "cycle")
            for patient_index, donor in enumerate(cycle):
                weight = 1
                next_patient = cycle[(patient_index + 1) % len(cycle)]
                found = False
                for edge in donor.edges:
                    if edge.tgt == next_patient:
                        weight = edge.score
                        found = True
                        break
                    if found:
                        break
                cycle_weight += weight
                for edge in next_patient.edges:
                    if edge.tgt == donor:
                        backarcs += 1
                pair = ET.SubElement(cycle_elt, "pair")
                # The pair is the donor/patient pair who joined the exchange
                # together, not the donor + the patient to whom they are
                # donating
                ET.SubElement(pair, "p").text = "%s" % donor.patient_id()
                ET.SubElement(pair, "d").text = "%s" % donor.donor_id()
            total_weight += cycle_weight
            cycle_elt.set("id", "%d" % cycle_index)
            cycle_elt.set("backarcs", "%d" % backarcs)
            cycle_elt.set("weight", "%f" % cycle_weight)
            cycle_index += 1
        return cycle_index, total_weight

    def _xml_add_chains(self, cycles, cycle_index):
        """Private function to add chains to an XML element.
        """
        total_weight = 0
        if not self._ndds:
            return cycle_index, total_weight
        for chain in self.chains:
            chain_weight = 0
            cycle_elt = ET.SubElement(cycles, "cycle")
            donor = self._ndds[chain.ndd_index]
            first = True
            for patient_index in chain.vertices():
                weight = 1
                patient = self.digraph.vs[patient_index]
                donor_index = False
                for edge in donor.edges:
                    if edge.target() == patient:
                        weight = edge.score
                        donor_index = edge.donor().index()
                        break
                    if donor_index:
                        break
                chain_weight += weight
                pair = ET.SubElement(cycle_elt, "pair")
                if first:
                    ET.SubElement(pair, "a").text = "true"
                    first = False
                else:
                    ET.SubElement(pair, "p").text = "%s" % donor.patient_id()
                ET.SubElement(pair, "d").text = "%s" % donor_index
                donor = self.digraph.vs[patient_index]
            # donor now points to the last donor/patient in the chain.
            pair = ET.SubElement(cycle_elt, "pair")
            ET.SubElement(pair, "p").text = "%s" % donor.patient_id()
            ET.SubElement(pair, "d").text = "%s" % donor.donor_id()

            total_weight += chain_weight
            cycle_elt.set("id", "%d" % cycle_index)
            cycle_elt.set("weight", "%f" % chain_weight)
            cycle_elt.set("altruistic", "true")
            cycle_index += 1
        return cycle_index, total_weight

    def relabelled_copy(self, old_to_new_vertices, new_digraph):
        """Create a copy of the solution with vertices relabelled.

        If the solution was found on a relabelled copy of the instance digraph, this
        method can be used to transform the solution back to the original digraph. Each
        Vertex v in the OptSolution on which this method is called is replaced in the
        returned copy by old_to_new_vertices[v.index()].
        """

        relabelled_cycles = [[old_to_new_vertices[v.index()] for v in c] for c in self.cycles]
        relabelled_chains = [Chain(c.ndd_index,
                                   [old_to_new_vertices[i].index() for i in c.vtx_indices],
                                   c.score)
                             for c in self.chains]
        return OptSolution(self.ip_model, relabelled_cycles, relabelled_chains,
                           new_digraph, edge_success_prob=self.edge_success_prob)

def optimise(model, cfg):
    """Optimise a model."""
    model.update()
    LOGGER.debug("Number of variables: %d", model.numVars)
    LOGGER.debug("Number of constraints: %d", model.numConstrs)
    LOGGER.debug("Number of non-zero coefficients: %d", model.numNZs)
    if cfg.details:
        with open(cfg.details, "a") as outfile:
            outfile.write("%d,%d,%d,%d\n" % (cfg.numCycles, model.numVars,
                                             model.numConstrs, model.numNZs))
    elif cfg.lp_file:
        model.write(cfg.lp_file)
    elif cfg.relax:
        relaxed = model.relax()
        relaxed.optimize()
        print("lp_relax_obj_val:", relaxed.obj_val)
        print("lp_relax_solver_status:", relaxed.status)
        sys.exit(0)
    else:
        LOGGER.debug("Starting optimiser")
        model.optimize()
        LOGGER.debug("Optimizer finished")

def optimise_relabelled(formulation_fun, cfg):
    """Optimise on a relabelled graph such that vertices are sorted in descending
        order of (indegree + outdegree)"""

    in_degs = [0] * cfg.digraph.n
    for edge in cfg.digraph.es:
        in_degs[edge.tgt.index()] += 1

    sorted_vertices = sorted(cfg.digraph.vs,
                             key=lambda v: len(v.edges) + in_degs[v.index()],
                             reverse=True)

    relabelled_digraph = cfg.digraph.induced_subgraph(sorted_vertices)

    # old_to_new_vtx[i] is the vertex in the new graph corresponding to vertex
    # i in the original digraph
    old_to_new_vtx = [None] * cfg.digraph.n
    for index, vertex in enumerate(sorted_vertices):
        old_to_new_vtx[vertex.index()] = relabelled_digraph.vs[index]
    relabelled_ndds = create_relabelled_ndds(cfg.ndds, old_to_new_vtx)
    relabelled_cfg = copy.copy(cfg)
    relabelled_cfg.digraph = relabelled_digraph
    relabelled_cfg.ndds = relabelled_ndds

    opt_result = formulation_fun(relabelled_cfg)
    return opt_result.relabelled_copy(sorted_vertices, cfg.digraph)

def create_ip_model(time_limit, verbose):
    """Create a Gurobi Model."""

    model = Model("kidney-mip")
    if not verbose:
        model.params.outputflag = 0
    model.params.mipGap = 0
    if time_limit is not None:
        model.params.timelimit = time_limit
    return model

###################################################################################################
#                                                                                                 #
#                                       Uncapped formulation                                      #
#                                                                                                 #
###################################################################################################

def add_unlimited_vars_and_constraints(digraph, ndds, m):
    """Add the IP variables and constraints for chains in the uncapped edge formulation. 

    Args:
        digraph: the instance digraph
        ndds: a list of NDDs in the instance
        m: The Gurobi model
    """

    for v in digraph.vs:
        v.grb_vars_in  = []
        v.grb_vars_out = []

    for ndd in ndds:
        ndd_edge_vars = []
        for e in ndd.edges:
            edge_var = m.addVar(vtype=GRB.BINARY)
            e.grb_vars = [edge_var]
            e.grb_var_positions = [0]
            e.edge_var = edge_var
            ndd_edge_vars.append(edge_var)
            e.target_v.grb_vars_in.append(edge_var)
        m.update()
        m.addConstr(quicksum(ndd_edge_vars) <= 1)

    # Add pair->pair edge variables
    for e in digraph.es:
        e.grb_vars = []
        edge_var = m.addVar(vtype=GRB.BINARY)
        e.grb_vars.append(edge_var)
        e.src.grb_vars_out.append(edge_var)
        e.tgt.grb_vars_in.append(edge_var)

    m.update()

    for v in digraph.vs:
        if len(v.grb_vars_in) > 1:
            m.addConstr(quicksum(v.grb_vars_in) <= 1)

    # Sum of edges into a vertex must be >= sum of edges out
    for v in digraph.vs:
        m.addConstr(quicksum(v.grb_vars_in) >= quicksum(v.grb_vars_out))

def optimise_uuef(cfg):
    """Optimise using the uncapped edge formulation.

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    if cfg.edge_success_prob != 1:
        raise ValueError("This formulation does not support failure-aware matching.")

    m = create_ip_model(cfg.timelimit, cfg.verbose)

    add_unlimited_vars_and_constraints(cfg.digraph, cfg.ndds, m)

    obj_expr = ( quicksum(e.score * e.edge_var for ndd in cfg.ndds for e in ndd.edges) +
                 quicksum(e.score * var for e in cfg.digraph.es for var in e.grb_vars) )

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    # Try all possible cycle start positions
    cycle_start_vv = range(cfg.digraph.n)

    cycle_next_vv = {}
    for e in cfg.digraph.es:
        for var in e.grb_vars:
            if var.x > 0.1:
                cycle_next_vv[e.src.index()] = e.tgt.index()

    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(
                                    cfg.digraph, cycle_start_vv, cycle_next_vv),
                       chains=kidney_utils.get_optimal_chains(cfg.digraph, cfg.ndds),
                       digraph=cfg.digraph)
        
###################################################################################################
#                                                                                                 #
#                  Chain vars and constraints (used by HPIEF', HPIEF'' and PICEF)                 #
#                                                                                                 #
###################################################################################################

def add_chain_vars_and_constraints(digraph, ndds, max_chain, m, vtx_to_vars,
                                   store_edge_positions=False,
                                   store_collapse_chains=False):
    """Add the IP variables and constraints for chains in PICEF and HPIEF'.

    Args:
        ndds: a list of NDDs in the instance
        max_chain: the chain cap
        m: The Gurobi model
        vtx_to_vars: A list such that for each Vertex v in the Digraph,
            vtx_to_vars[v.index()] will contain the Gurobi variables representing
            edges pointing to v.
        store_edge_positions: if this is True, then an attribute grb_edge_positions
            will be added to edges that have associated Gurobi variables.
            edge.grb_edge_positions[i] will indicate the position of the edge respresented
            by edge.grb_vars[i]. (default: False)
        store_collapse_chains: If True, create the variables (and corresponding
            constraints) which will indicate when a vertex is the last vertex in a
            chain, and the chain can be collapsed by removing the first
            donor/patient pair. Note that this only works if max_chain == 3.
            The interesting variable will be attached to each vertex as
            collapse_var, as a binary variable, and will be 1 if and only if
            this vertex is the last vertex of a chain selected that can
            collapse.
    """

    if max_chain == 0:
        return
    if store_collapse_chains and max_chain != 3:
        raise Exception("store_collapse_chains can only be used if max_chain == 3")

    for vert in digraph.vs:
        vert.grb_vars_in = [[] for i in range(max_chain-1)]
        vert.grb_vars_out = [[] for i in range(max_chain-1)]
        if store_collapse_chains:
            vert.collapse_var = m.addVar(vtype=GRB.BINARY)
            vert.collapse_vars = []

    for index, ndd in enumerate(ndds):
        ndd_edge_vars = []
        for edge in ndd.edges:
            name = "ce_alt_%s" % edge.donor().index()
            name = name.replace("-", "_")
            edge_var = m.addVar(vtype=GRB.BINARY, name=name)
            edge.edge_var = edge_var
            edge.grb_vars = [edge_var]
            edge.grb_var_positions = [0]
            ndd_edge_vars.append(edge_var)
            vtx_to_vars[edge.target_v.index()].append(edge_var)
            if max_chain > 1:
                edge.target_v.grb_vars_in[0].append(edge_var)
        m.update()
        m.addConstr(quicksum(ndd_edge_vars) <= 1, "ndd_%d" % index)

    dists_from_ndd = kidney_utils.get_dist_from_nearest_ndd(digraph, ndds)

    # Add pair->pair edge variables, indexed by position in chain
    for edge in digraph.es:
        edge.grb_vars = []
        if store_edge_positions:
            edge.grb_var_positions = []
        for i in range(max_chain-1):
            if dists_from_ndd[edge.src.index()] <= i+1:
                name = "ce_%s_%s" % (edge, i)
                name = name.replace("-", "_")
                edge_var = m.addVar(vtype=GRB.BINARY, name=name)
                edge.grb_vars.append(edge_var)
                if store_edge_positions:
                    edge.grb_var_positions.append(i+1)
                vtx_to_vars[edge.tgt.index()].append(edge_var)
                edge.src.grb_vars_out[i].append(edge_var)
                if i < max_chain-2:
                    edge.tgt.grb_vars_in[i+1].append(edge_var)

    if store_collapse_chains:
        for ndd in ndds:
            for first_edge in ndd.edges:
                try:
                    first_edge.edge_var
                except AttributeError:
                    continue
                for second_edge in first_edge.target_v.edges:
                    try:
                        second_edge.edge_var
                    except AttributeError:
                        continue
                    for third_edge in second_edge.tgt.edges:
                        try:
                            third_edge.edge_var
                        except AttributeError:
                            continue
                        collapse_var = m.addVar(vtype=GRB.BINARY)
                        third_edge.tgt.collapse_vars.append(collapse_var)
                        m.addConstr(3*collapse_var <= first_edge.edge_var +
                                    second_edge.edge_var + third_edge.edge_var)
        for vert in digraph.vs:
            if vert.collapse_vars:
                m.addConstr(quicksum(vert.collapse_vars) >= vert.collapse_var)
    m.update()

    # At each chain position, sum of edges into a vertex must be >= sum of edges out
    for i in range(max_chain-1):
        for vert in digraph.vs:
            m.addConstr(quicksum(vert.grb_vars_in[i]) >= quicksum(vert.grb_vars_out[i]), "chain_%d_%d" % (vert.index(), i))

###################################################################################################
#                                                                                                 #
#                                Code shared by HPIEF' and HPIEF''                                #
#                                                                                                 #
###################################################################################################

def add_hpief_prime_vars_partial_red(max_cycle, digraph, m, hpief_2_prime=False):
    vars_and_edges = [] # A list of (gurobi_var, position, edge, low_vertex) tuples

    # max_pos is the maximum edge position for which variables may be created
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1
    
    # Index i is in the list edge_vars_in[pos][v][low_v] if and only if
    # vars_and_edges[i] corresponds to an edge at position pos, pointing to vertex
    # v, in low_v's graph copy 
    edge_vars_in = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    # Index i is in the list edge_vars_out[pos][v][low_v] if and only if
    # vars_and_edges[i] corresponds to an edge at position pos, leaving vertex
    # v, in low_v's graph copy 
    edge_vars_out = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    for low_vtx in range(digraph.n-1):
        # Length of shortest path from low vertex to each vertex with a higher index
        # Default value is 999999999 (which represents infinity)
        shortest_path_from_lv = digraph.get_shortest_path_from_low_vtx(low_vtx, max_cycle-1) 
        shortest_path_to_lv = digraph.get_shortest_path_to_low_vtx(low_vtx, max_cycle-1) 

        for v1 in digraph.vs[low_vtx+1:]:
            for e in v1.edges:
                if e.tgt.index() >=low_vtx:
                    for pos in range(1, max_pos + 1):
                        if (shortest_path_from_lv[e.src.index()] <= pos and
                                    shortest_path_to_lv[e.tgt.index()] < max_cycle - pos):
                            new_var = m.addVar(vtype=GRB.BINARY)
                            vars_and_edges.append((new_var, pos, e, low_vtx))
                            idx = len(vars_and_edges) - 1 # Index of tuple just added
                            edge_vars_in[pos][e.tgt.index()][low_vtx].append(idx)
                            edge_vars_out[pos][e.src.index()][low_vtx].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_hpief_prime_vars_full_red(max_cycle, digraph, m, hpief_2_prime=False):
    vars_and_edges = [] # A list of (gurobi_var, position, edge, low_vertex) tuples

    # max_pos is the maximum edge position for which variables may be created
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1
    
    edge_vars_in = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]
    edge_vars_out = [[[[] for __ in range(digraph.n)] for __ in range(digraph.n)] for __ in range(max_pos + 1)]

    edges_seen = set()  # (low_v_id, src_v_id, tgt_v_id, pos) tuples
    for cycle in digraph.generate_cycles(max_cycle):
        for i in range(1, len(cycle)-1):
            edges_seen.add((cycle[0].index(), cycle[i].index(), cycle[i+1].index(), i))
        if not hpief_2_prime or len(cycle) < max_cycle:
            edges_seen.add((cycle[0].index(), cycle[-1].index(), cycle[0].index(), len(cycle)-1))
            
    for low_v, src_v, tgt_v, pos in edges_seen:
        new_var = m.addVar(vtype=GRB.BINARY)
        e = digraph.adj_mat[src_v][tgt_v]
        vars_and_edges.append((new_var, pos, e, low_v))
        idx = len(vars_and_edges) - 1 # Index of tuple just added
        edge_vars_in[pos][tgt_v][low_v].append(idx)
        edge_vars_out[pos][src_v][low_v].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_hpief_prime_vars_and_constraints(max_cycle, digraph, vtx_to_in_edges, m, full_red, hpief_2_prime=False):
    max_pos = max_cycle-2 if hpief_2_prime else max_cycle-1

    if full_red:
        vars_and_edges, edge_vars_in, edge_vars_out = add_hpief_prime_vars_full_red(max_cycle, digraph, m, hpief_2_prime)
    else:
        vars_and_edges, edge_vars_in, edge_vars_out = add_hpief_prime_vars_partial_red(max_cycle, digraph, m, hpief_2_prime)
    
    for grb_var, pos, edge, low_vtx in vars_and_edges:
        vtx_to_in_edges[edge.tgt.index()].append(grb_var)
        if pos==1:
            vtx_to_in_edges[edge.src.index()].append(grb_var)
        if hpief_2_prime and pos == max_cycle - 2 and edge.tgt.index() != low_vtx:
            vtx_to_in_edges[low_vtx].append(grb_var)
        
    # Capacity constraint for vertices
    for l in vtx_to_in_edges:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)
    
    # Cycle flow-conservation constraint for vertices
    for pos in range(1, max_pos):
        for v in range(digraph.n):
            for low_v_id in range(v):
                in_vars  = [vars_and_edges[i][0] for i in edge_vars_in[pos][v][low_v_id]]
                out_vars = [vars_and_edges[i][0] for i in edge_vars_out[pos+1][v][low_v_id]]
                if len(in_vars) > 0 or len(out_vars) > 0:
                    m.addConstr(quicksum(in_vars) == quicksum(out_vars))

    return vars_and_edges

def optimise_hpief_prime(cfg, full_red=False, hpief_2_prime=False):
    """Optimise using the HPIEF' or HPIEF'' formulation.

    The HPIEF' model is based on HPIEF, but does not include cycle-edge variables at position zero.
    HPIEF'' also removes variables corresponding to edges at the last possible position of a cycle. 

    Args:
        cfg: an OptConfig object
        full_red: True if cycles should be generated in order to reduce number of variables further
        hpief_2_prime: Use HPIEF''? Default: HPIEF'

    Returns:
        an OptSolution object
    """

    if cfg.edge_success_prob != 1:
        raise ValueError("This formulation does not support failure-aware matching.")

    if cfg.max_cycle < 3:
        hpief_2_prime = False

    m = create_ip_model(cfg.timelimit, cfg.verbose)
    m.params.method = 2
    m.params.presolve = 0

    # For each vertex v, a list of variables corresponding to in-edges to v
    vtx_to_in_edges = [[] for __ in cfg.digraph.vs]

    add_chain_vars_and_constraints(cfg.digraph, cfg.ndds, cfg.max_chain, m, vtx_to_in_edges)

    vars_and_edges = add_hpief_prime_vars_and_constraints(
            cfg.max_cycle, cfg.digraph, vtx_to_in_edges, m, full_red, hpief_2_prime)

    obj_terms = []
    for var, pos, edge, low_v_id in vars_and_edges:
        score = edge.score
        if pos==1:
            score += cfg.digraph.adj_mat[low_v_id][edge.src.index()].score
        if hpief_2_prime and pos==cfg.max_cycle - 2 and edge.tgt.index() != low_v_id:
            score += cfg.digraph.adj_mat[edge.tgt.index()][low_v_id].score
        obj_terms.append(score * var)

    obj_expr = quicksum(obj_terms)
   
    if cfg.max_chain > 0:
        obj_expr += quicksum(e.score * e.edge_var for ndd in cfg.ndds for e in ndd.edges) 
        obj_expr += quicksum(e.score * var for e in cfg.digraph.es for var in e.grb_vars)
    
    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    cycle_start_vv = []
    cycle_next_vv = {}
    
    for var, pos, edge, low_v_id in vars_and_edges:
        if var.x > 0.1:
            cycle_next_vv[edge.src.index()] = edge.tgt.index()
            if pos == 1:
                cycle_start_vv.append(low_v_id)
                cycle_next_vv[low_v_id] = edge.src.index()
            if hpief_2_prime and pos == cfg.max_cycle - 2 and edge.tgt.index() != low_v_id:
                cycle_next_vv[edge.tgt.index()] = low_v_id
        
    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(
                                    cfg.digraph, cycle_start_vv, cycle_next_vv),
                       chains=[] if cfg.max_chain==0 else kidney_utils.get_optimal_chains(cfg.digraph, cfg.ndds),
                       digraph=cfg.digraph)

###################################################################################################
#                                                                                                 #
#                                               HPIEF'                                            #
#                                                                                                 #
###################################################################################################

def optimise_hpief_prime_full_red(cfg):
    return optimise_hpief_prime(cfg, True)

###################################################################################################
#                                                                                                 #
#                                             HPIEF''                                             #
#                                                                                                 #
###################################################################################################

def optimise_hpief_2prime(cfg, full_red=False):
    return optimise_hpief_prime(cfg, full_red, hpief_2_prime=True)

def optimise_hpief_2prime_full_red(cfg):
    return optimise_hpief_2prime(cfg, full_red=True)

###################################################################################################
#                                                                                                 #
#                                              PICEF                                              #
#                                                                                                 #
###################################################################################################

def optimise_picef(cfg):
    """Optimise using the PICEF formulation.

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """
    LOGGER.info("Using PICEF formulation")
    model = create_ip_model(cfg.timelimit, cfg.verbose)
    model.params.method = 2

    vtx_to_vars = [[] for __ in cfg.digraph.vs]

    LOGGER.info("Adding chain variables and constraints")
    add_chain_vars_and_constraints(cfg.digraph, cfg.ndds, cfg.max_chain, model, vtx_to_vars,
                                   store_edge_positions=cfg.edge_success_prob != 1)

    # vtx_to_vars has a list for each vertex, that list contains the cycles
    # incident with said vertex
    LOGGER.info("Finding cycles")
    cycle_vars = []
    for cycle in cfg.digraph.generate_cycles(cfg.max_cycle):
        var = model.addVar(vtype=GRB.BINARY)
        cycle_vars.append((var, cycle))
        for vert in cycle:
            vtx_to_vars[vert.index()].append(var)
    LOGGER.info("Number of cycles: %d", len(cycle_vars))
    cfg.numCycles = len(cycle_vars)
    LOGGER.info("Adding vertex constraints")
    for var_list in vtx_to_vars:
        if var_list:
            model.addConstr(quicksum(var_list) <= 1)

    LOGGER.info("Calculating objective")
    if cfg.max_chain == 0:
        if cfg.only_size():
            obj_expr = quicksum(len(cycle) * cfg.edge_success_prob ** len(cycle) * var
                                for var, cycle in cycle_vars)
        else:
            obj_expr = quicksum(cycle_score(cycle, cfg.edge_success_prob) * var
                                for var, cycle in cycle_vars)
    elif cfg.edge_success_prob == 1:
        if cfg.only_size():
            obj_expr = (quicksum(len(cycle) * var for var, cycle in cycle_vars) +
                        quicksum(e.edge_var for ndd in cfg.ndds for e in ndd.edges) +
                        quicksum(var for e in cfg.digraph.es for var in e.grb_vars))
        else:
            obj_expr = (quicksum(cycle_score(cycle, cfg.edge_success_prob) * var
                                 for var, cycle in cycle_vars) +
                        quicksum(e.score * e.edge_var for ndd in cfg.ndds
                                 for e in ndd.edges) +
                        quicksum(e.score * var for e in cfg.digraph.es
                                 for var in e.grb_vars))
    else:
        if cfg.only_size():
            obj_expr = (quicksum(len(cycle) * cfg.edge_success_prob ** len(cycle)  * var
                                 for var, cycle in cycle_vars) +
                        quicksum(cfg.edge_success_prob * e.edge_var
                                 for ndd in cfg.ndds for e in ndd.edges) +
                        quicksum(cfg.edge_success_prob**(pos+1) * var
                                 for e in cfg.digraph.es
                                 for var, pos in zip(e.grb_vars, e.grb_var_positions)))
        else:
            obj_expr = (quicksum(cycle_score(cycle, cfg.edge_success_prob) * var
                                 for var, cycle in cycle_vars) +
                        quicksum(e.score*cfg.edge_success_prob * e.edge_var
                                 for ndd in cfg.ndds for e in ndd.edges) +
                        quicksum(e.score*cfg.edge_success_prob**(pos+1) * var
                                 for e in cfg.digraph.es
                                 for var, pos in zip(e.grb_vars, e.grb_var_positions)))

    model.setObjective(obj_expr, GRB.MAXIMIZE)

    optimise(model, cfg)

    if model.status != GRB.OPTIMAL:
        return None

    return OptSolution(ip_model=model,
                       cycles=[cycle for var, cycle in cycle_vars if var.x > 0.5],
                       chains=[] if cfg.max_chain == 0 else
                       kidney_utils.get_optimal_chains(cfg.digraph, cfg.ndds,
                                                       cfg.edge_success_prob),
                       digraph=cfg.digraph, ndds=cfg.ndds,
                       edge_success_prob=cfg.edge_success_prob)

def optimise_picef_nhs(cfg):
    """Optimise using the PICEF formulation, but use the NHS objectives.

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """
    LOGGER.info("Using PICEF formulation with NHS objectives")
    model = create_ip_model(cfg.timelimit, cfg.verbose)
    model.params.method = 2
    LOGGER.debug("Adding cycle variables")
    cycle_vars = []
    two_cycles = []
    three_cycles = []
    three_cycles_with_backarcs = []
    vtx_to_vars = [[] for __ in cfg.digraph.vs]
    matchable_constraints = {}
    if cfg._constrain_maximal:
        LOGGER.debug("Getting max matchable edges")
        max_matchable = list(cfg.digraph.get_max_matchable_edges_networkx(cfg.ndds))
        LOGGER.debug([str(s) for s in max_matchable])
        for e in max_matchable:
            matchable_constraints[e] = []
            LOGGER.debug("%s is necessary", e)
        LOGGER.info("Constraining %d maximally matchable edges" % len(max_matchable))
    for index, cycle in enumerate(cfg.digraph.generate_cycles(cfg.max_cycle)):
        name = "_".join(str(v) for v in cycle)
        var = model.addVar(vtype=GRB.BINARY, name=name)
        cycle_vars.append((var, cycle))
        for vert in cycle:
            vtx_to_vars[vert.index()].append(var)
        # Find 2-cycles and 3-cycles with back-arcs
        if len(cycle) == 2:
            two_cycles.append(var)
        elif len(cycle) == 3:
            three_cycles.append(var)
            backarc = False
            for index, vert1 in enumerate(cycle):
                vert2 = cycle[index - 1 % len(cycle)]
                for edge in vert1.edges:
                    if edge.tgt == vert2:
                        backarc = True
                        break
                if backarc:
                    three_cycles_with_backarcs.append(var)
                    break
        if cfg._constrain_maximal:
            for edge in max_matchable:
                if cycle_contains(cycle, edge):
                    matchable_constraints[edge].append(var)

    LOGGER.info("Number of cycles: %d", len(cycle_vars))
    cfg.numCycles = len(cycle_vars)
    model.update()


    LOGGER.debug("Adding chain variables and constraints")
    add_chain_vars_and_constraints(cfg.digraph, cfg.ndds, cfg.max_chain, model, vtx_to_vars,
                                   store_edge_positions=True, store_collapse_chains=False)
    if cfg._constrain_maximal:
        for edge in max_matchable:
            matchable_constraints[edge].extend(edge.grb_vars)
        for key, list_of_vars in matchable_constraints.items():
            if list_of_vars:
                name = "must_use_%s" % key
                name = name.replace("-","_")
                model.addConstr(quicksum(set(list_of_vars)) == 1, name=name)

    LOGGER.debug("Adding vertex constraints")
    for index, var_list in enumerate(vtx_to_vars):
        if var_list:
            model.addConstr(quicksum(var_list) <= 1, "vert_%d" % index)


    first_edges = []
    second_edges = []
    for ndd in cfg.ndds:
        for edge in ndd.edges:
            first_edges.append(edge.edge_var)
    for vert in cfg.digraph.vs:
        if cfg.max_chain > 1:
            for var in vert.grb_vars_out[0]:
                second_edges.append(var)
    short_chains = quicksum(first_edges) - quicksum(second_edges)
    long_chains = quicksum(second_edges)


    # First just maximise 2-cycles, 3-cycles with backarcs, short
    # chains and collapsible long chains
    obj = (quicksum(two_cycles) +
        quicksum(three_cycles_with_backarcs) + short_chains + long_chains)
    model.setObjective(obj, GRB.MAXIMIZE)
    optimise(model, cfg)
    if model.status != GRB.OPTIMAL:
        LOGGER.warning("First round not found optimal: %s.", model.status)
        return None
    obj1_value = model.getObjective().getValue()
    LOGGER.debug("First round complete, objective is %d", obj1_value)
    # Add constraint keeping the first objective at its optimal value
    model.addConstr(obj >= obj1_value, "obj1_limit")

    # Next maximise number of transplants
    obj = (quicksum(len(cycle) * var for (var, cycle) in cycle_vars) +
        quicksum(e.edge_var for ndd in cfg.ndds for e in ndd.edges) +
        quicksum(var for e in cfg.digraph.es for var in e.grb_vars))
    model.setObjective(obj, GRB.MAXIMIZE)
    if cfg._constrain_maximal:
        model.write("matchable.mps")
        model.write("matchable.lp")
    else:
        model.write("default.mps")
        model.write("default.lp")
    optimise(model, cfg)
    if model.status != GRB.OPTIMAL:
        return None
    obj2_value = model.getObjective().getValue()
    LOGGER.debug("Second round complete, size is %d", obj2_value)
    # Add constraint on second objective
    model.addConstr(obj >= obj2_value)

    # Third objective, minimize three-ways
    obj = quicksum(three_cycles)
    model.setObjective(obj, GRB.MINIMIZE)
    optimise(model, cfg)
    if model.status != GRB.OPTIMAL:
        return None
    obj3_value = model.getObjective().getValue()
    LOGGER.debug("Third round complete, size is %d", obj3_value)
    model.addConstr(obj <= obj3_value)

    # Fourth objective, maximise three-ways with back-arcs and collapsible
    # long chains
    obj = quicksum(three_cycles_with_backarcs) + quicksum(two_cycles)
    model.setObjective(obj, GRB.MAXIMIZE)
    optimise(model, cfg)
    if model.status != GRB.OPTIMAL:
        return None
    obj4_value = model.getObjective().getValue()
    LOGGER.debug("Fourth round complete, size is %d", obj4_value)
    model.addConstr(obj >= obj4_value)

    # This must now be the last iteration
    obj_expr = quicksum(cycle_score(cycle, cfg.edge_success_prob, nhs=True)
                        * var for (var, cycle) in cycle_vars)
    if cfg.max_chain != 0:
        for vertex in cfg.ndds + cfg.digraph.vs:
            for chain_one in vertex.edges:
                next_vert = chain_one.target()
                next_donor = chain_one.target().donor()
                for chain_var_one, chain_one_index in zip(chain_one.grb_vars, chain_one.grb_var_positions):
                    obj_expr += chain_one.age_formula(next_donor) * cfg.edge_success_prob**(chain_one_index+1) * chain_var_one
                    for chain_two in next_vert.edges:
                        for chain_var_two, chain_two_index in zip(chain_two.grb_vars, chain_two.grb_var_positions):
                            # Ensure that chain_two can follow chain_one
                            if chain_one_index + 1 != chain_two_index:
                                continue
                            combined = model.addVar(vtype=GRB.BINARY)
                            model.addConstr(-2 * combined + chain_var_one + chain_var_two >= 0)

    model.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(model, cfg)

    if model.status != GRB.OPTIMAL:
        return None
    obj5_value = model.getObjective().getValue()
    LOGGER.debug("Fifth round complete, score is %d", obj5_value)

    return OptSolution(ip_model=model,
                       cycles=[cycle for (var, cycle) in cycle_vars if var.x > 0.5],
                       chains=[] if cfg.max_chain==0 else kidney_utils.get_optimal_chains(
                            cfg.digraph, cfg.ndds, cfg.edge_success_prob, nhs=True),
                       digraph=cfg.digraph, ndds=cfg.ndds,
                       edge_success_prob=cfg.edge_success_prob, nhs=True)


#def optimise_picef_nhs_chains(cfg):
#    """Optimise using the PICEF formulation, but use the NHS objectives.
#    Uses alternative definition of fail-able chain.
#
#    Args:
#        cfg: an OptConfig object
#
#    Returns:
#        an OptSolution object
#    """
#    LOGGER.info("Using PICEF formulation with NHS objectives")
#    LOGGER.info("Finding cycles")
#    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
#    model = create_ip_model(cfg.timelimit, cfg.verbose)
#    model.params.method = 2
#
#    LOGGER.info("Adding cycle variables")
#    cycle_vars = [model.addVar(vtype=GRB.BINARY) for __ in cycles]
#    model.update()
#
#    vtx_to_vars = [[] for __ in cfg.digraph.vs]
#
#    LOGGER.info("Adding chain variables and constraints")
#    add_chain_vars_and_constraints(cfg.digraph, cfg.ndds, cfg.max_chain, model, vtx_to_vars,
#                                   store_edge_positions=cfg.edge_success_prob != 1,
#                                   store_collapse_chains=True)
#
#    # vtx_to_vars has a list for each vertex, that list contains the cycles
#    # incident with said vertex
#    for index, cycle in enumerate(cycles):
#        for vert in cycle:
#            vtx_to_vars[vert.index()].append(cycle_vars[index])
#
#    LOGGER.info("Adding vertex constraints")
#    for var_list in vtx_to_vars:
#        if var_list:
#            model.addConstr(quicksum(var_list) <= 1)
#
#    # Note that for 2-cycles and 3-cycles, we store a list of tuples, each
#    # tuple containing (cycle, gurobi_variable_for_cycle).
#    # Find 2-cycles
#    two_cycles = [(cycle, var) for cycle, var in zip(cycles, cycle_vars)
#                  if len(cycle) == 2]
#    # Find 3-cycles with back-arcs
#    three_cycles_with_backarcs = []
#    for cycle, var in zip(cycles, cycle_vars):
#        if len(cycle) != 3:
#            continue
#        backarc = False
#        for index, vert1 in enumerate(cycle):
#            vert2 = cycle[index - 1 % len(cycle)]
#            for edge in vert1.edges:
#                if edge.tgt == vert2:
#                    backarc = True
#                    break
#            if backarc:
#                three_cycles_with_backarcs.append((cycle, var))
#                break
#    second_edges = []
#    third_edges = []
#    long_chain_collapse_vars = []
#    for vert in cfg.digraph.vs:
#        for edge in vert.edges:
#            try:
#                edge.grb_vars_out
#            except AttributeError:
#                continue
#            for var in edge.grb_vars_out[2]:
#                second_edges.append(var)
#            for var in edge.grb_vars_out[3]:
#                third_edges.append(var)
#        long_chain_collapse_vars.append(vert.collapse_var)
#    short_chains = quicksum(second_edges) - quicksum(third_edges)
#    long_collapse_chains = quicksum(long_chain_collapse_vars)
#    long_chains = quicksum(third_edges)
#
#    for iteration in range(5):
#        LOGGER.info("Calculating objective in round %d", iteration)
#        # First just maximise 2-cycles, 3-cycles with backarcs, short
#        # chains and collapsible long chains
#        obj = (quicksum([var for (cycle, var) in two_cycles]) +
#               quicksum([var for (cycle, var) in three_cycles_with_backarcs]) +
#               short_chains + long_collapse_chains)
#        if iteration == 0:
#            model.setObjective(obj, GRB.MAXIMIZE)
#            optimise(model, cfg)
#            if model.status != GRB.OPTIMAL:
#                return None
#            obj1_value = model.getObjective().getValue()
#            LOGGER.info("First round complete, objective is %d", obj1_value)
#            continue
#        # Add constraint keeping the first objective at its optimal value
#        model.addConstr(obj == obj1_value)
#
#        # Next maximise number of transplants
#        if cfg.max_chain == 0:
#            obj = quicksum(len(c) * cfg.edge_success_prob ** len(c) * var for
#                           c, var in zip(cycles, cycle_vars))
#        elif cfg.edge_success_prob == 1:
#            obj = (quicksum(len(c) * var
#                            for c, var in zip(cycles, cycle_vars)) +
#                   quicksum(e.edge_var for ndd in cfg.ndds
#                            for e in ndd.edges) +
#                   quicksum(var for e in cfg.digraph.es for var in e.grb_vars))
#        else:
#            obj = (quicksum(len(c) * cfg.edge_success_prob ** len(c)  * var for
#                            c, var in zip(cycles, cycle_vars)) +
#                   quicksum(cfg.edge_success_prob * e.edge_var for ndd in
#                            cfg.ndds for e in ndd.edges) +
#                   quicksum(cfg.edge_success_prob**(pos+1) * var for e in
#                            cfg.digraph.es for var, pos in zip(e.grb_vars,
#                                                               e.grb_var_positions)))
#        if iteration == 1:
#            model.setObjective(obj, GRB.MAXIMIZE)
#            optimise(model, cfg)
#            if model.status != GRB.OPTIMAL:
#                return None
#            obj2_value = model.getObjective().getValue()
#            LOGGER.info("Second round complete, size is %d", obj2_value)
#            continue
#        # Add constraint on second objective
#        model.addConstr(obj == obj2_value)
#
#        # Third objective, minimize three-ways + long chains
#        obj = long_chains + quicksum(var for c, var in zip(cycles, cycle_vars)
#                                     if len(c) == 3)
#        if iteration == 2:
#            model.setObjective(obj, GRB.MINIMIZE)
#            optimise(model, cfg)
#            if model.status != GRB.OPTIMAL:
#                return None
#            obj3_value = model.getObjective().getValue()
#            LOGGER.info("Third round complete, size is %d", obj3_value)
#            continue
#
#        model.addConstr(obj == obj3_value)
#
#        # Fourth objective, maximise three-ways with back-arcs and collapsible
#        # long chains
#        obj = (quicksum([var for (cycle, var) in three_cycles_with_backarcs]) +
#               long_collapse_chains)
#        if iteration == 3:
#            model.setObjective(obj, GRB.MAXIMIZE)
#            optimise(model, cfg)
#            if model.status != GRB.OPTIMAL:
#                return None
#            obj4_value = model.getObjective().getValue()
#            LOGGER.info("Fourth round complete, size is %d", obj4_value)
#            continue
#        model.addConstr(obj == obj4_value)
#        # This must now be the last iteration
#        if cfg.max_chain == 0:
#            obj_expr = quicksum(cycle_score(c, cfg.edge_success_prob) * var
#                                for c, var in zip(cycles, cycle_vars))
#        else:
#            obj_expr = (quicksum(cycle_score(c, cfg.edge_success_prob) * var
#                                 for c, var in zip(cycles, cycle_vars)) +
#                        quicksum(e.score*cfg.edge_success_prob * e.edge_var
#                                 for ndd in cfg.ndds for e in ndd.edges) +
#                        quicksum(e.score*cfg.edge_success_prob**(pos+1) * var
#                                 for e in cfg.digraph.es
#                                 for var, pos in zip(e.grb_vars, e.grb_var_positions)))
#
#        model.setObjective(obj_expr, GRB.MAXIMIZE)
#        optimise(model, cfg)
#
#        if model.status != GRB.OPTIMAL:
#            return None
#        obj5_value = model.getObjective().getValue()
#        LOGGER.info("Fifth round complete, size is %d", obj5_value)
#
#    return OptSolution(ip_model=model,
#                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
#                       chains=[] if cfg.max_chain==0 else kidney_utils.get_optimal_chains(
#                            cfg.digraph, cfg.ndds, cfg.edge_success_prob),
#                       digraph=cfg.digraph, ndds=cfg.ndds,
#                       edge_success_prob=cfg.edge_success_prob)

###################################################################################################
#                                                                                                 #
#                                        Cycle formulation                                        #
#                                                                                                 #
###################################################################################################

def optimise_ccf(cfg):
    """Optimise using the cycle formulation (with one var per cycle and one var per chain).

    Args:
        cfg: an OptConfig object

    Returns:
        an OptSolution object
    """

    LOGGER.info("Finding cycles")
    cycles = cfg.digraph.find_cycles(cfg.max_cycle)
    LOGGER.info("Finding chains")
    chains = find_chains(cfg.digraph, cfg.ndds, cfg.max_chain, cfg.edge_success_prob)

    m = create_ip_model(cfg.timelimit, cfg.verbose)
    m.params.method = 2

    cycle_vars = [m.addVar(vtype=GRB.BINARY, name="cycle_%d"%index) for index, __ in enumerate(cycles)]
    chain_vars = [m.addVar(vtype=GRB.BINARY, name="chain_%d"%index) for index, __ in enumerate(chains)]
    m.update()

    ndd_to_vars = [[] for __ in cfg.ndds]
    vtx_to_vars = [[] for __ in cfg.digraph.vs]

    for var, c in zip(cycle_vars, cycles):
        for v in c:
            vtx_to_vars[v.index()].append(var)

    for var, c in zip(chain_vars, chains):
        ndd_to_vars[c.ndd_index].append(var)
        for v in c.vtx_indices:
            vtx_to_vars[v].append(var)

    LOGGER.info("Adding vertex constraints")
    # Each donor-patient pair and each each NDD is in at most one chosen cycle or chain
    for l in vtx_to_vars + ndd_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    LOGGER.info("Calculating objective")
    if cfg.only_size():
        obj_expr = (quicksum(len(c) * cfg.edge_success_prob ** len(c) * var
                             for (c, var) in zip(cycles, cycle_vars)) +
                    quicksum(len(c) * var for (c, var) in zip(chains, chain_vars)))
    else:
        obj_expr = (quicksum(cycle_score(c, cfg.edge_success_prob) * var
                             for (c, var) in zip(cycles, cycle_vars)) +
                    quicksum(c.score * var for (c, var) in zip(chains, chain_vars)))

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    return OptSolution(ip_model=m,
                       cycles=[c for c, v in zip(cycles, cycle_vars) if v.x > 0.5],
                       chains=[c for c, v in zip(chains, chain_vars) if v.x > 0.5],
                       digraph=cfg.digraph, ndds=cfg.ndds,
                       edge_success_prob=cfg.edge_success_prob)

###################################################################################################
#                                                                                                 #
#                                    Extended Edge Formulation                                    # 
#                                                                                                 #
###################################################################################################

def add_eef_vars_partial_red(max_cycle, digraph, m):
    vars_and_edges = [] # A list of (gurobi_var, edge, low_vertex) tuples

    # Index i is in the list edge_vars_in[low_v][v] if and only if
    # vars_and_edges[i] corresponds to an edge pointing to vertex v, in low_v's graph copy 
    edge_vars_in = [[[] for __ in range(digraph.n)] for __ in range(digraph.n)]

    # Index i is in the list edge_vars_out[low_v][v] if and only if
    # vars_and_edges[i] corresponds to an edge leaving vertex v, in low_v's graph copy 
    edge_vars_out = [[[] for __ in range(digraph.n)] for __ in range(digraph.n)]

    for low_vtx in range(digraph.n-1):
        # Length of shortest path from low vertex to each vertex with a higher index
        # Default value is 999999999 (which represents infinity)
        shortest_path_from_lv = digraph.get_shortest_path_from_low_vtx(low_vtx, max_cycle-1) 
        shortest_path_to_lv = digraph.get_shortest_path_to_low_vtx(low_vtx, max_cycle-1) 

        for v1 in digraph.vs[low_vtx:]:
            for e in v1.edges:
                if e.tgt.index() >=low_vtx:
                    if (shortest_path_from_lv[e.src.index()] +
                                shortest_path_to_lv[e.tgt.index()] < max_cycle):
                        new_var = m.addVar(vtype=GRB.BINARY)
                        vars_and_edges.append((new_var, e, low_vtx))
                        idx = len(vars_and_edges) - 1 # Index of tuple just added
                        edge_vars_in[low_vtx][e.tgt.index()].append(idx)
                        edge_vars_out[low_vtx][e.src.index()].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_eef_vars_full_red(max_cycle, digraph, m):
    vars_and_edges = [] # A list of (gurobi_var, edge, low_vertex) tuples

    edge_vars_in = [[[] for __ in range(digraph.n)] for __ in range(digraph.n)]
    edge_vars_out = [[[] for __ in range(digraph.n)] for __ in range(digraph.n)]

    edges_seen = set()  # (low_v_id, src_v_id, tgt_v_id) tuples
    for cycle in digraph.generate_cycles(max_cycle):
        for i in range(len(cycle)):
            edges_seen.add((cycle[0].index(), cycle[i-1].index(), cycle[i].index()))
            
    for low_v, src_v, tgt_v in edges_seen:
        new_var = m.addVar(vtype=GRB.BINARY)
        e = digraph.adj_mat[src_v][tgt_v]
        vars_and_edges.append((new_var, e, low_v))
        idx = len(vars_and_edges) - 1 # Index of tuple just added
        edge_vars_in[low_v][tgt_v].append(idx)
        edge_vars_out[low_v][src_v].append(idx)
    m.update()
    return vars_and_edges, edge_vars_in, edge_vars_out

def add_eef_vars_and_constraints(max_cycle, digraph, m, full_red, eef_alt_constraints, vtx_to_in_edges):
    if full_red:
        vars_and_edges, edge_vars_in, edge_vars_out = add_eef_vars_full_red(max_cycle, digraph, m)
    else:
        vars_and_edges, edge_vars_in, edge_vars_out = add_eef_vars_partial_red(max_cycle, digraph, m)
    
    for grb_var, edge, low_vtx in vars_and_edges:
        vtx_to_in_edges[edge.tgt.index()].append(grb_var)
        
    # Capacity constraint for vertices
    for l in vtx_to_in_edges:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)
    
    # Cycle flow-conservation constraint for vertices
    for v in range(digraph.n):
        for low_v_id in range(v):
            in_vars  = [vars_and_edges[i][0] for i in edge_vars_in[low_v_id][v]]
            out_vars = [vars_and_edges[i][0] for i in edge_vars_out[low_v_id][v]]
            if len(in_vars) > 0 or len(out_vars) > 0:
                m.addConstr(quicksum(in_vars) == quicksum(out_vars))

    if eef_alt_constraints:
        for low_v_id in range(v):
            edge_indices_in_graph_copy = [i for indices in edge_vars_in[low_v_id] for i in indices]

            edge_vars_leaving_l = []
            edge_vars_not_involving_l = []    # Edge vars where low_v_id is neither the src nor the tgt
            for i in edge_indices_in_graph_copy:
                var, edge, _ = vars_and_edges[i]
                if edge.src.index()==low_v_id:
                    edge_vars_leaving_l.append(var)
                elif edge.tgt.index()!=low_v_id:
                    edge_vars_not_involving_l.append(var)
            
            # Number of edges constraint for each graph copy
            # (Note that this is redundant, but removing it seems to slow the program down
            # quite a bit.)
            m.addConstr(quicksum(edge_vars_not_involving_l) <= max_cycle-2) 

            # In each graph copy, if any edge is selected then an edge is selected
            # that leaves the low-numbered vertex in the graph copy
            # Note: this differs from (9e) in Constantino et al.
            m.addConstr(quicksum(edge_vars_not_involving_l) <=
                        (max_cycle-2) * quicksum(edge_vars_leaving_l))

    else:
        for low_v_id in range(v):
            edge_indices_in_graph_copy = [i for indices in edge_vars_in[low_v_id] for i in indices]

            # Number of edges constraint for each graph copy
            m.addConstr(quicksum(vars_and_edges[i][0] for i in edge_indices_in_graph_copy) <= max_cycle)        

            # Constraint (9e) from Constantino et al.
            sum_of_edge_vars_leaving_l = quicksum(
                    vars_and_edges[i][0] for i in edge_vars_out[low_v_id][low_v_id])
            for i in range(low_v_id+1, digraph.n):
                vars_leaving_i = [vars_and_edges[j][0] for j in edge_vars_out[low_v_id][i]]
                if len(vars_leaving_i):
                    m.addConstr(quicksum(vars_leaving_i) <= sum_of_edge_vars_leaving_l)


    return vars_and_edges

def optimise_eef(cfg, full_red=False):
    """Optimise using the reduced extended edge formulation (Constantino et al., EJOR, 2013).

    Note that this implementation does not yet include chains, and throws an exception
    if a chain cap greater than zero is used.

    Args:
        cfg: an OptConfig object
        full_red: True if cycles should be generated in order to reduce number of variables further

    Returns:
        an OptSolution object
    """

    if cfg.edge_success_prob != 1:
        raise ValueError("This formulation does not support failure-aware matching.")

    m = create_ip_model(cfg.timelimit, cfg.verbose)
    m.params.method = 2
    m.params.presolve = 0

    # For each vertex v, a list of variables corresponding to in-edges to v
    vtx_to_in_edges = [[] for __ in cfg.digraph.vs]

    add_chain_vars_and_constraints(cfg.digraph, cfg.ndds, cfg.max_chain, m, vtx_to_in_edges)

    vars_and_edges = add_eef_vars_and_constraints(cfg.max_cycle, cfg.digraph, m, full_red,
                                                  cfg.eef_alt_constraints, vtx_to_in_edges)

    obj_expr = quicksum(edge.score * var for var, edge, low_v_id in vars_and_edges)
    if cfg.max_chain > 0:
        obj_expr += quicksum(e.score * e.edge_var for ndd in cfg.ndds for e in ndd.edges) 
        obj_expr += quicksum(e.score * var for e in cfg.digraph.es for var in e.grb_vars)

    m.setObjective(obj_expr, GRB.MAXIMIZE)
    optimise(m, cfg)

    cycle_start_vv = []
    cycle_next_vv = {}
    
    for var, edge, low_v_id in vars_and_edges:
        if var.x > 0.1:
            cycle_next_vv[edge.src.index()] = edge.tgt.index()
            cycle_start_vv.append(edge.src.index())
        
    return OptSolution(ip_model=m,
                       cycles=kidney_utils.selected_edges_to_cycles(
                                    cfg.digraph, cycle_start_vv, cycle_next_vv),
                       chains=[] if cfg.max_chain==0 else kidney_utils.get_optimal_chains(cfg.digraph, cfg.ndds),
                       digraph=cfg.digraph)

def optimise_eef_full_red(cfg):
    return optimise_eef(cfg, full_red=True)
