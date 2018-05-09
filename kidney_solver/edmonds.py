"""Edmond's Blossom algorithm.
"""


import logging

from kidney_solver.forest import Forest


LOGGER = logging.getLogger(__name__)


def get_blossom_path(blossom, entry, leave):
    """Given a blossom (odd length cycle) and an leave and entry point, find the
    even length route around the cycle.

    :param blossom: the cycle, as a list of vertices
    :param entry, leave: the entry and leave points around the cycle

    :return: the even length path from entry to leave around the cycle
    """

    # Build up the backwards path in two stages, one from the "first" vertex in
    # the list, and one from the "last" vertex in the list
    back = []
    backtwo = []
    fore = []
    # Stage = 1 means before first of leave/entry, stage 2 means between the
    # two, stage 3 means after the final leave/entry
    stage = 1
    for vert in blossom:
        if stage == 1:
            back.append(vert)
            if vert == entry or vert == leave:
                fore.append(vert)
                stage = 2
                continue
        if stage == 2:
            fore.append(vert)
            if vert == entry or vert == leave:
                backtwo.append(vert)
                stage = 3
                continue
        if stage == 3:
            backtwo.append(vert)
    if len(fore) % 2 == 0:
        return fore
    back.reverse()
    backtwo.reverse()
    return back + backtwo


class Graph(object):
    """A graph which we can use to find a maximum matching."""
    def __init__(self):
        self._nodes = {}
        self._edges = []

    def add_node(self, node):
        """Add a node."""
        # List of adjacent vertices.
        self._nodes[node] = []

    def add_edge(self, edge):
        """Add an edge."""
        if not isinstance(edge, tuple):
            raise TypeError
        self._edges.append(edge)
        self._nodes[edge[0]].append(edge[1])
        self._nodes[edge[1]].append(edge[0])

    def remove_edge(self, edge):
        """Remove an edge."""
        self._edges.remove(edge)
        self._nodes[edge[0]].remove(edge[1])
        self._nodes[edge[1]].remove(edge[0])

    def __str__(self):
        s = "V: " + " ".join(str(x) for x in self._nodes)
        s += "\n"
        s += "E: " + " ".join(str(x) for x in self._edges)
        return s

    def degree(self, node):
        """Return the degree of a node."""
        return len(self._nodes[node])

    def find_a_maximum_matching(self):
        """Find a maximum matching"""
        matching = []
        while True:
            path, forced = self.find_augmenting_path(matching)
            if not path:
                break
            LOGGER.debug("Found path " + str(path))
            # Follow augmenting path, which starts and ends with open edges not in
            # the matching, and alternates with edges that currently are in the
            # matching.
            for edge in path:
                if edge in matching:
                    matching.remove(edge)
                else:
                    matching.append(edge)
        return matching

    def find_maximally_matchable_edges(self):
        """Find all edges which are in every maximum sized matching."""
        largest = len(self.find_a_maximum_matching())
        matchable = []
        for e in list(self._edges):
            self.remove_edge(e)
            if len(self.find_a_maximum_matching()) < largest:
                matchable.append(e)
            self.add_edge(e)
        return matchable

    def find_augmenting_path(self, matching=None):
        """Finds an augmenting path in the graph.

        :param matching: A matching as a list of edges
        :returns: An augmenting path as a list of edges, or an empty list if no
        augmenting path exists.
        """
        forest = Forest()
        marked = {}
        marked_v = {}
        exposed = {}
        for vertex in self._nodes.keys():
            exposed[vertex] = True
        for m in matching:
            marked[m] = True
            if m[0] in exposed:
                del exposed[m[0]]
            if m[1] in exposed:
                del exposed[m[1]]
        for vertex in exposed.keys():
            LOGGER.debug("" + str(vertex) + " is exposed")
            forest.add_node(vertex)
        v = None
        for v in forest.even_height_nodes():
            if forest.distance(v, forest.root(v)) % 2 != 0:
                continue
            if v in marked_v:
                continue
            break
        keep_going = v is not None
        LOGGER.debug("Found " + str(v))
        while keep_going:
            for w in self._nodes.keys():
                LOGGER.debug("Looking at " + str(v) + ", " + str(w))
                if w == v:
                    continue
                if (v, w) not in self._edges and (w, v) not in self._edges:
                    LOGGER.debug("not edge")
                    continue
                if (v, w) in marked or (w, v) in marked:
                    LOGGER.debug("marked")
                    continue
                LOGGER.debug("Interesting " + str(v) + ", " + str(w))
                if not forest.has_node(w):
                    LOGGER.debug("w is matched")
                    # w is matched in M
                    for m in matching:
                        if m[0] == w:
                            forest.add_edge(v, w)
                            forest.add_edge(w, m[1])
                            break
                        if m[1] == w:
                            forest.add_edge(v, w)
                            forest.add_edge(w, m[0])
                            break
                else:
                    LOGGER.debug("w is not matched")
                    if forest.distance(w, forest.root(w)) % 2 == 1:
                        # Do nothing
                        pass
                    else:
                        if forest.root(v) != forest.root(w):
                            # Report augmenting path
                            forced = self.degree(v) == 1 and self.degree(w) == 1
                            path = forest.path_to_root(v) + [(v, w)] + forest.path_to_root(w)
                            return path, forced
                        else:
                            # Blossom on e + path from v to w
                            blossom = [(v, w)] + forest.path_between(v, w)
                            # Convert blossom to list of vertices only
                            # Shrink graph (creating new?)
                            g = Graph()
                            new_node = len(self._nodes)
                            translates = {}
                            for n in self._nodes:
                                g.add_node(n)
                            g.add_node(new_node)
                            for edge in self._edges:
                                if edge[0] in blossom and edge[1] in blossom:
                                    continue
                                if edge[0] in blossom:
                                    g.add_edge((edge[1], new_node))
                                    translates[(edge[1], new_node)] = edge
                                elif edge[1] in blossom:
                                    g.add_edge((edge[0], new_node))
                                    translates[(edge[0], new_node)] = edge
                                else:
                                    g.add_edge(edge)
                            # Also convert matching!
                            new_matching = []
                            for edge in matching:
                                if edge[0] in blossom and edge[1] in blossom:
                                    continue
                                if edge[0] in blossom:
                                    new_matching.append((edge[1], new_node))
                                elif edge[1] in blossom:
                                    new_matching.append((edge[0], new_node))
                                else:
                                    new_matching.append(edge)
                            # Find augmenting path in shrunken graph
                            path, forced = g.find_augmenting_path(new_matching)
                            # Expand path onto "this" graph
                            real_path = []
                            entry = None
                            if path is not None:
                                for edge in path:
                                    if edge in translates:
                                        translated = translates[edge]
                                        # Add blossom bit

                                        # If entry is None, we don't know where
                                        # we "enter" the blossom, so just mark
                                        # that.
                                        if entry is None:
                                            real_path.append(translated)
                                            if translated[0] == new_node:
                                                entry = translated[1]
                                            else:
                                                entry = translated[0]
                                        else:
                                            # We know where we entered the blossom,
                                            # so we either go "forwards" or
                                            # "backwards" around the cycle,
                                            # whichever is even length (knowing
                                            # that the cycle must have odd length).
                                            if translated[0] == new_node:
                                                leave = translated[1]
                                            else:
                                                leave = translated[0]
                                            bpath = get_blossom_path(blossom,
                                                                    leave, entry)
                                            real_path.extend(bpath)
                                            real_path.append(translated)
                                    else:
                                        real_path.append(edge)
                            return real_path, forced
                marked[(v, w)] = True
            marked_v[v] = True
            keep_going = False
            for v in forest.even_height_nodes():
                if forest.distance(v, forest.root(v)) % 2 != 0:
                    continue
                if v in marked_v:
                    continue
                keep_going = True
                break
        return None, True
