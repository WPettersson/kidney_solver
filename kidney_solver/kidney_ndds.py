"""This module has objects for non-directed donors (NDDs) and chains initiated by NDDs.

In contrast to many kidney-exchange formulations, we do not include in the main directed
graph a vertex for each NDD. Rather, the program maintains a separate array
of Ndd objects. Each of these Ndd objects maintains a list of outgoing edges, and
each of these edges points towards a vertex (which represents a donor-patient pair)
in the directed graph.
"""

class Ndd(object):
    """A non-directed donor"""
    def __init__(self):
        self.edges = []
        self._donor_id = None

    def add_edge(self, target, score, explanation):
        """Add an edge representing compatibility with a patient who appears as a
        vertex in the directed graph."""
        self.edges.append(NddEdge(target, score, explanation))


class NddEdge:
    """An edge pointing from an NDD to a vertex in the directed graph"""
    def __init__(self, target_v, score, donor):
        self.target_v = target_v
        self.score = score
        self._donor = donor
        self.grb_vars = []
        self.grb_var_positions = []

    def donor(self):
        """The donor corresponding to this transplant."""
        return self._donor

    def target(self):
        """The target of this edge."""
        return self.target_v

def create_relabelled_ndds(ndds, old_to_new_vtx):
    """Creates a copy of a n array of NDDs, with target vertices changed.

    If a target vertex in an original NDD had ID i, then the new target vertex
    will be old_to_new_vtx[i].
    """
    new_ndds = [Ndd() for ndd in ndds]
    for i, ndd in enumerate(ndds):
        for edge in ndd.edges:
            new_ndds[i].add_edge(old_to_new_vtx[edge.target().index()], edge.score, edge.donor())

    return new_ndds

class Chain(object):
    """A chain initiated by an NDD.
    
    Data members:
        ndd_index: The index of the NDD
        vtx_indices: The indices of the vertices in the chain, in order
        score: the chain's score
    """

    def __init__(self, ndd_index, vtx_indices, score):
        self.ndd_index = ndd_index
        self.vtx_indices = vtx_indices
        self.score = score

    def __repr__(self):
        return ("Chain NDD{} ".format(self.ndd_index) +
                        " ".join(str(v) for v in self.vtx_indices) +
                        " with score " + str(self.score))

    def __cmp__(self, other):
        # Compare on NDD ID, then chain length, then on score, then
        # lexicographically on vtx indices
        if self.ndd_index < other.ndd_index:
            return -1
        elif self.ndd_index > other.ndd_index:
            return 1
        elif len(self.vtx_indices) < len(other.vtx_indices):
            return -1
        elif len(self.vtx_indices) > len(other.vtx_indices):
            return 1
        elif self.score < other.score:
            return -1
        elif self.score > other.score:
            return 1
        else:
            for i, j in zip(self.vtx_indices, other.vtx_indices):
                if i < j:
                    return -1
                elif i > j:
                    return 1
        return 0

    def vertices(self):
        """The list of vertices in this chain"""
        return self.vtx_indices

    def __len__(self):
        return len(self.vtx_indices)
            
def find_chains(digraph, ndds, max_chain, edge_success_prob=1):
    """Generate all chains with up to max_chain edges."""

    def find_chains_recurse(vertices, score):
        chains.append(Chain(ndd_idx, vertices[:], score))
        if len(vertices) < max_chain:
            for e in digraph.vs[vertices[-1]].edges:
                if e.tgt.index() not in vertices:
                    vertices.append(e.tgt.index())
                    find_chains_recurse(vertices, score+e.score*edge_success_prob**len(vertices))
                    del vertices[-1]
    chains = []
    if max_chain == 0:
        return chains
    for ndd_idx, ndd in enumerate(ndds):
        for e in ndd.edges:
            vertices = [e.target_v.index()]
            find_chains_recurse(vertices, e.score*edge_success_prob)
    return chains

