"""Edmond's Blossom algorithm.
"""


import logging

from kidney_solver.forest import Forest

LOGGER = logging.getLogger(__name__)
#LOGGER.setLevel(logging.INFO)

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
        self._biggest = -1

    def add_node(self, node):
        """Add a node."""
        # List of adjacent vertices.
        if node not in self._nodes:
            self._nodes[node] = []
            if node > self._biggest:
                self._biggest = node

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

    def find_a_maximum_matching(self, matching=None):
        """Find a maximum matching"""
        if matching is None:
            matching = []
            # Initialise with a greedy matching
            known = []
            for edge in self._edges:
                if edge[0] not in known and edge[1] not in known:
                    matching.append(edge)
                    known.extend([edge[0], edge[1]])
        while True:
            path, forced = self.find_augmenting_path(matching)
            if not path:
                break
            for index, edge in enumerate(path):
                # Follow augmenting path, which starts and ends with open edges not in
                # the matching, and alternates with edges that currently are in the
                # matching.
                if index % 2 == 1:
                    # Edge to be removed from matching
                    try:
                        matching.remove(edge)
                    except ValueError as e:
                        LOGGER.error("Tried to remove %s from %s", edge, sorted(matching, key=lambda x:x[0]))
                        raise e
                else:
                    matching.append(edge)
        return matching

    def find_maximally_matchable_edges(self):
        """Find all edges which are in every maximum sized matching."""
        largest = self.find_a_maximum_matching()
        LOGGER.info("Largest matching has size %d", len(largest))
        matchable = []
        interesting = list(largest)
        while interesting:
            e = interesting.pop()
            if e not in self._edges:
                e = (e[1], e[0])
            LOGGER.debug("Testing %s", e)
            self.remove_edge(e)
            new_matching = list(largest)
            if e in new_matching:
                new_matching.remove(e)
            new_max = self.find_a_maximum_matching(new_matching)
            if len(new_max) < len(largest):
                LOGGER.debug("new matching has size %d", len(new_max))
                matchable.append(e)
                LOGGER.debug("%s is matchable", e)
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
            #marked_v[m[0]] = True
            #marked_v[m[1]] = True
            if m[0] in exposed:
                del exposed[m[0]]
            if m[1] in exposed:
                del exposed[m[1]]
        for vertex in exposed.keys():
            forest.add_node(vertex)
        v = None
        for v in forest.even_height_nodes():
            if forest.distance(v, forest.root(v)) % 2 != 0:
                continue
            if v in marked_v:
                continue
            break
        keep_going = v is not None
        LOGGER.debug("Looking for edge to " + str(v))
        while keep_going:
            for w in self._nodes.keys():
                if w == v:
                    continue
                if (v, w) not in self._edges and (w, v) not in self._edges:
                    continue
                if (v, w) in marked or (w, v) in marked:
                    continue
                if not forest.has_node(w):
                    #LOGGER.debug("w is matched")
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
                    #LOGGER.debug("w is not matched")
                    if forest.distance(w, forest.root(w)) % 2 == 1:
                        # Do nothing
                        pass
                    else:
                        if forest.root(v) != forest.root(w):
                            # Report augmenting path
                            forced = self.degree(v) == 1 and self.degree(w) == 1
                            if v < w:
                                pathA = forest.path_to_root(v)
                                pathA.reverse()
                                path = pathA + [(v, w)] + forest.path_to_root(w)
                            else:
                                pathA = forest.path_to_root(v)
                                pathA.reverse()
                                path = pathA + [(w, v)] + forest.path_to_root(w)
                            LOGGER.debug("Found a path %s", path)
                            return path, forced
                        else:
                            # Blossom on e + path from v to w
                            blossom = [(v, w)] + forest.path_between(v, w)
                            # Convert blossom to list of vertices only
                            # Shrink graph (creating new?)
                            g = Graph()
                            LOGGER.debug("Shrinking on %s", blossom)
                            new_node = self._biggest + 1
                            LOGGER.debug("Adding new node %d" % new_node)
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
                                    if edge[1] == 490:
                                        LOGGER.debug("Adding %s", (edge[1], new_node))
                                elif edge[1] in blossom:
                                    g.add_edge((edge[0], new_node))
                                    translates[(edge[0], new_node)] = edge
                                    if edge[0] == 490:
                                        LOGGER.debug("Adding %s", (edge[0], new_node))
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
                            LOGGER.debug("Found untranslated path %s", path)
                            if path is not None:
                                for edge in path:
                                    if edge in translates:
                                        translated = translates[edge]
                                        LOGGER.debug("%s is translated from %s", edge, translated)
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
                            LOGGER.debug("Found translated path %s", real_path)
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
