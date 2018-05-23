"""A Digraph class which can be used for representing donor-patient pairs
(as vertices) and their compatibilities (as weighted edges), along with
some related methods.
"""

from collections import deque
import logging

from networkx.algorithms.matching import max_weight_matching
from networkx import Graph as nx_Graph


from kidney_solver.edmonds import Graph


LOGGER = logging.getLogger(__name__)


# Constants for the age formula
AGE_SELECTOR = 20
AGE_WEIGHT = 3
AGE_FORMULA_FACTOR = 1e-5
AGE_FORMULA_MAX_AGE_DIFF = 70


def cycle_score(cycle, edge_success_prob=1, nhs=False):
    """Calculate the sum of a cycle's edge scores.

    Args:
        cycle: A list of Vertex objects in the cycle, with the first Vertex not repeated.
        digraph: The digraph in which this cycle appears.
    """
    total = 0
    for index, current_vertex in enumerate(cycle):
        next_vertex = cycle[(index + 1) % len(cycle)]
        # Find edge
        edge = None
        for edge in current_vertex.edges:
            if edge.target() == next_vertex:
                break
        if not nhs:
            total += edge.score * edge_success_prob
            continue
        # NHS scoring needs next donation, to get the donor-donor age
        # difference, which is part of the weight calculations
        next_plus_one_vertex = cycle[(index + 2) % len(cycle)]
        # Find edge
        next_edge = None
        for next_edge in next_vertex.edges:
            if next_edge.target() == next_plus_one_vertex:
                break
        total += edge.age_formula(next_edge.donor()) * edge_success_prob
    return total


def cycle_contains(cycle, edge):
    """Return True if the given cycle contains the given edge. The cycle is a
    list of vertices as given by generate_cycles and the edge is an Edge
    object.
    """
    for index, v in enumerate(cycle):
        if v.index() == edge.source().index() and cycle[(index+1) % len(cycle)].index() == edge.target().index():
            return True
        if v.index() == edge.source().index() and cycle[(index-1) % len(cycle)].index() == edge.target().index():
            return True
    return False


class Vertex(object):
    """A vertex in a directed graph (see the Digraph class)."""

    def __init__(self, index):
        self._index = index
        self.edges = []
        self._donor_id = -1
        self._patient_id = -1

    def set_patient_id(self, patient_id):
        """Set the patient ID of this vertex."""
        self._patient_id = patient_id

    def patient_id(self):
        """Get the patient ID of the patient-patient-pair represented by this
        vertex.
        """
        return self._patient_id

    def set_donor_id(self, donor_id):
        """Set the donor ID of this vertex."""
        self._donor_id = donor_id

    def donor_id(self):
        """Get the donor ID of the donor-patient-pair represented by this
        vertex.
        """
        return self._donor_id

    def __str__(self):
        return "V{}".format(self._index)

    def index(self):
        """The index of the vertex."""
        return self._index

class Edge(object):
    """An edge in a directed graph (see the Digraph class)."""

    def __init__(self, id, score, src, tgt, donor):
        self.id = id
        self.score = score
        self.src = src   # source vertex
        self.tgt = tgt # target vertex
        self._donor = donor
    def __str__(self):
        return ("V" + str(self.src.index()) + "-V" + str(self.tgt.index()))

    def source(self):
        """The source of this edge."""
        return self.src

    def target(self):
        """The target of this edge."""
        return self.tgt

    def donor(self):
        """Which donor is donating for this edge."""
        return self._donor

    def age_formula(self, next_donor):
        """Calculate the weight of this edge, taking into account the next
        donor and the age difference between the two donors.
        """
        weight = AGE_WEIGHT
        age_diff = abs(self.donor().age() - next_donor.age())
        if age_diff > AGE_SELECTOR:
            weight = 0
        age_formula = (AGE_FORMULA_MAX_AGE_DIFF - age_diff) ** 2
        return self.score + weight + age_formula * AGE_FORMULA_FACTOR

class Digraph(object):
    """A directed graph, in which each edge has a numeric score.

    Data members:
        n: the number of vertices in the digraph
        vs: an array of Vertex objects, such that vs[i].index() == i
        es: an array of Edge objects, such that es[i].index() = i
    """

    def __init__(self, n):
        """Create a Digraph with n vertices"""
        self.n = n
        self.vs = [Vertex(i) for i in range(n)]
        self.adj_mat = [[None for x in range(n)] for x in range(n)]
        self.es = []

    def add_edge(self, score, source, tgt, donor):
        """Add an edge to the digraph

        Args:
            score: the edge's score, as a float
            source: the source Vertex
            tgt: the edge's target Vertex
            donor: the Donor for this edge
        """

        id = len(self.es)
        e = Edge(id, score, source, tgt, donor)
        self.es.append(e)
        source.edges.append(e)
        self.adj_mat[source.index()][tgt.index()] = e
    
    def find_cycles(self, max_length):
        """Find cycles of length up to max_length in the digraph.

        Returns:
            a list of cycles. Each cycle is represented as a list of
            vertices, with the first vertex _not_ repeated at the end.
        """

        return [cycle for cycle in self.generate_cycles(max_length)]

    def generate_cycles(self, max_length):
        """Generate cycles of length up to max_length in the digraph.

        Each cycle yielded by this generator is represented as a list of
        vertices, with the first vertex _not_ repeated at the end.
        """

        vtx_used = [False] * len(self.vs)  # vtx_used[i]==True iff vertex i is in current path

        def cycle(current_path):
            last_vtx = current_path[-1]
            if self.edge_exists(last_vtx, current_path[0]):
                yield current_path[:]
            if len(current_path) < max_length:
                for edge in last_vtx.edges:
                    v = edge.tgt
                    if (len(current_path) + shortest_paths_to_low_vtx[v.index()] <= max_length
                            and not vtx_used[v.index()]):
                        current_path.append(v)
                        vtx_used[v.index()] = True
                        for c in cycle(current_path):
                            yield c
                        vtx_used[v.index()] = False
                        del current_path[-1]

        # Adjacency lists for transpose graph
        transp_adj_lists = [[] for v in self.vs]
        for edge in self.es:
            transp_adj_lists[edge.tgt.index()].append(edge.src)

        for v in self.vs:
            shortest_paths_to_low_vtx = self.calculate_shortest_path_lengths(
                    v, max_length - 1,
                    lambda u: (w for w in transp_adj_lists[u.index()] if w.index() > v.index()))
            vtx_used[v.index()] = True
            for c in cycle([v]):
                yield c
            vtx_used[v.index()] = False

    def get_shortest_path_from_low_vtx(self, low_vtx, max_path):
        """ Returns an array of path lengths. For each v > low_vtx, if the shortest
            path from low_vtx to v is shorter than max_path, then element v of the array
            will be the length of this shortest path. Otherwise, element v will be
            999999999."""
        return self.calculate_shortest_path_lengths(self.vs[low_vtx], max_path,
                                                    adj_list_accessor=lambda v:
                                                    (e.tgt for e in v.edges if
                                                     e.tgt.index() >= low_vtx))

    def get_shortest_path_to_low_vtx(self, low_vtx, max_path):
        """ Returns an array of path lengths. For each v > low_vtx, if the shortest
            path to low_vtx from v is shorter than max_path, then element v of the array
            will be the length of this shortest path. Otherwise, element v will be
            999999999."""
        def adj_list_accessor(v):
            for i in range(low_vtx, len(self.vs)):
                if self.adj_mat[i][v.index()]:
                    yield self.vs[i]

        return self.calculate_shortest_path_lengths(self.vs[low_vtx], max_path,
                                                    adj_list_accessor=adj_list_accessor)

    def calculate_shortest_path_lengths(self, from_v, max_dist,
                                        adj_list_accessor=lambda v: (e.tgt for
                                                                     e in
                                                                     v.edges)):
        """Calculate the length of the shortest path from vertex from_v to each
        vertex with a greater or equal index, using paths containing
        only vertices indexed greater than or equal to from_v.

        Return value: a list of distances of length equal to the number of vertices.
        If the shortest path to a vertex is greater than max_dist, the list element
        will be 999999999.

        Args:
            from_v: The starting vertex
            max_dist: The maximum distance we're interested in
            adj_list_accessor: A function taking a vertex and returning an
                iterable of out-edge targets
        """
        # Breadth-first search
        q = deque([from_v])
        distances = [999999999] * len(self.vs)
        distances[from_v.index()] = 0

        while q:
            v = q.popleft()
            #Note: >= is used instead of == on next line in case max_dist<0
            if distances[v.index()] >= max_dist:
                break
            for w in adj_list_accessor(v):
                if distances[w.index()] == 999999999:
                    distances[w.index()] = distances[v.index()] + 1
                    q.append(w)

        return distances

    def num_edges(self):
        """The number of edges (arcs) in the graph."""
        return len(self.es)

    def edge_exists(self, v1, v2):
        """Returns true if and only if an edge exists from Vertex v1 to Vertex v2."""

        return self.adj_mat[v1.index()][v2.index()] is not None

    def update_edge(self, v1, v2, new_score, donor):
        """If an edge exists, update its score, but only if the new score is
        higher."""
        if not self.edge_exists(v1, v2):
            self.add_edge(v1, v2, new_score, donor)
        if self.adj_mat[v1.index()][v2.index()].score < new_score:
            self.adj_mat[v1.index()][v2.index()].score = new_score
            self.adj_mat[v1.index()][v2.index()]._donor = donor

    def induced_subgraph(self, vertices):
        """Returns the subgraph indiced by a given list of vertices."""

        subgraph = Digraph(len(vertices))
        for i, v in enumerate(vertices):
            for j, w in enumerate(vertices):
                e = self.adj_mat[v.index()][w.index()]
                if e is not None:
                    new_src = subgraph.vs[i]
                    new_tgt = subgraph.vs[j]
                    subgraph.add_edge(e.score, new_src, new_tgt, e.donor())
        return subgraph

    def __str__(self):
        return "\n".join([str(v) for v in self.vs])

    def get_max_matchable_edges_networkx(self, ndds=None):
        """Find the set of edges which must be in every maximum matching.
        """
        g = nx_Graph()
        translate = {}
        for edge in self.es:
            v1 = edge.source()
            v2 = edge.target()
            if v2.index() < v1.index():
                continue
            new_edge = (v1.index(), v2.index())
            if self.edge_exists(v2, v1) or edge.donor().is_altruistic():
                g.add_node(new_edge[0])
                g.add_node(new_edge[1])
                g.add_edge(v1.index(), v2.index())
                translate[new_edge] = edge
        count = len(translate)
        for ndd in ndds:
            for edge in ndd.edges:
                v1 = edge.donor()
                v2 = edge.target()
                new_edge = (v2.index(), count+v1.index())
                g.add_node(new_edge[0])
                g.add_node(new_edge[1])
                g.add_edge(v2.index(), count+v1.index())
                translate[new_edge] = edge
        # TODO Add NDD edges to this graph!
        largest = max_weight_matching(g)
        LOGGER.debug("Largest matching has size %d", len(largest))
        edges = []
        for v1, v2 in largest.items():
            if v1 < v2:
                edges.append([v1, v2])
        matchable = []
        while edges:
            v1, v2 = edges.pop()
            LOGGER.debug("Testing [%s, %s]", v1, v2)
            g.remove_edge(v1, v2)
            new_max = max_weight_matching(g)
            if len(new_max) < len(largest):
                LOGGER.debug("new matching has size %d", len(new_max))
                edges = list(filter(lambda x: x[0] in new_max, edges))
                matchable.append((v1, v2))
                LOGGER.debug("[%s, %s] is matchable", v1, v2)
            g.add_edge(v1, v2)
        LOGGER.info("Found %s maximally matchable edges" % len(matchable))
        return (translate[e] for e in matchable)


    def get_max_matchable_edges_self(self, ndds=None):
        """Find the set of edges which must be in every maximum matching.
        """
        # Build the simple graph
        g = Graph()
        translate = {}
        for edge in self.es:
            v1 = edge.source()
            v2 = edge.target()
            if v2.index() < v1.index():
                continue
            new_edge = (v1.index(), v2.index())
            if self.edge_exists(v2, v1) or edge.donor().is_altruistic():
                g.add_node(new_edge[0])
                g.add_node(new_edge[1])
                g.add_edge(new_edge)
                translate[new_edge] = edge
        count = len(translate)
        for ndd in ndds:
            for edge in ndd.edges:
                v1 = edge.donor()
                v2 = edge.target()
                new_edge = (v2.index(), count+v1.index())
                g.add_node(new_edge[0])
                g.add_node(new_edge[1])
                g.add_edge(new_edge)
                translate[new_edge] = edge
        # TODO Add NDD edges to this graph!
        matchable = g.find_maximally_matchable_edges()
        LOGGER.info("Found %s maximally matchable edges" % len(matchable))
        return (translate[e] for e in matchable)
