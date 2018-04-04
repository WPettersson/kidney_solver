"""An undirected forest, and associated classes. This is specifically designed
for Edmond's Blossom algorithm for finding maximum weight matchings in graphs,
and may not be suitable for other uses.
"""

class InvalidEdgeException(Exception):
    pass


class Forest(object):
    """An undirected graph."""
    def __init__(self):
        self._nodes = {}
        self._parent = {}
        self._root = {}
        self._distances = {}

    def add_node(self, node, height=0):
        """Add a node to the forest."""
        self._nodes[node] = height
        self._root[node] = node

    def add_edge(self, node_one, node_two):
        """Add an edge to the forest. Note that when adding the edge (a, b),
        the forest must already contain the node a.
        """
        if node_one not in self._nodes:
            raise InvalidEdgeException()
        self._nodes[node_two] = self._nodes[node_one] + 1
        self._root[node_two] = self._root[node_one]
        self._parent[node_two] = node_one

    def root(self, node):
        """Return the root of the given node."""
        return self._root[node]

    def has_node(self, node):
        """Does a node exist in this forest?"""
        return node in self._nodes

    def even_height_nodes(self):
        """Return a generator for even-height nodes.
        """
        return x for x in self._nodes if x % 2 == 0

    def height(node):
        """Get the height (distance to root) for a node.
        """
        return self._nodes[node]

    def distance(node_one, node_two):
        """Find the distance between two nodes."""
        # First just go up to the root from node one, storing the height along
        # the way
        parent = node_one
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
            parent = self._parent[parent]
        return 99999999999
