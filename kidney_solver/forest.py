"""An undirected forest, and associated classes. This is specifically designed
for Edmond's Blossom algorithm for finding maximum weight matchings in graphs,
and may not be suitable for other uses.
"""

import logging


class InvalidEdgeException(Exception):
    pass


LOGGER = logging.getLogger(__name__)


class Forest(object):
    """An undirected graph."""
    def __init__(self):
        self._nodes = {}
        self._parent = {}
        self._root = {}
        self._distances = {}

    def add_node(self, node):
        """Add a node to the forest."""
        self._nodes[node] = node
        self._parent[node] = None

    def add_edge(self, node_one, node_two):
        """Add an edge to the forest. Note that when adding the edge (a, b),
        the forest must already contain the node a.
        """
        if node_one not in self._nodes:
            raise InvalidEdgeException()
        self._nodes[node_two] = node_two
        self._parent[node_two] = node_one

    def root(self, node):
        """Return the root of the given node."""
        parent = self._parent[node]
        the_root = node
        while parent is not None:
            the_root = parent
            parent = self._parent[parent]
        return the_root

    def has_node(self, node):
        """Does a node exist in this forest?"""
        return node in self._nodes

    def even_height_nodes(self):
        """Return a generator for even-height nodes.
        """
        return (x for x in self._nodes if self.distance(x, self.root(x)) % 2 == 0)

    def path_between(self, one, two):
        """Find the path between two nodes."""
        # First just go up to the root from node one, storing the height along
        # the way
        parent = one
        if one == two:
            return []
        paths = {}
        path = []
        while parent:
            path.append(parent)
            paths[parent] = list(path)
            parent = self._parent[parent]
        # Now go up the parents of node_two until we find a parent of node_one
        parent = two
        path = []
        while parent:
            path.append(parent)
            if parent in paths:
                path.reverse()
                return paths[parent] + path
            parent = self._parent[parent]
        return None

    def distance(self, node_one, node_two):
        """Find the distance between two nodes."""
        # First just go up to the root from node one, storing the height along
        # the way
        parent = node_one
        if node_one == node_two:
            return 0
        parents = {}
        height = 0
        while parent:
            parents[parent] = height
            height += 1
            parent = self._parent[parent]
        # Now go up the parents of node_two until we find a parent of node_one
        parent = node_two
        height = 0
        while parent:
            if parent in parents:
                return height + parents[parent]
            height += 1
            parent = self._parent[parent]
        return 99999999999

    def path_to_root(self, node):
        """Return the path to the root from this node, as a list of edges."""
        now = node
        path = []
        parent = self._parent[now]
        while parent:
            path.append((now, parent))
            now = parent
            parent = self._parent[now]
        return path
