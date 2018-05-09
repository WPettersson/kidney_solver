
from nose.tools import raises, eq_, ok_
from kidney_solver.edmonds import Graph

def test_edmonds_one():
    g = Graph()
    g.add_node(0)
    g.add_node(1)
    g.add_edge((0, 1))
    eq_(len(g.find_a_maximum_matching()), 1)


def test_edmonds_two():
    g = Graph()
    g.add_node(0)
    g.add_node(1)
    g.add_node(2)
    g.add_edge((0, 1))
    g.add_edge((1, 2))
    g.add_edge((0, 2))
    eq_(len(g.find_a_maximum_matching()), 1)


def test_edmonds_three():
    g = Graph()
    g.add_node(0)
    g.add_node(1)
    g.add_node(2)
    g.add_node(3)
    g.add_edge((0, 1))
    g.add_edge((1, 2))
    g.add_edge((2, 3))
    eq_(len(g.find_a_maximum_matching()), 2)


def test_edmonds_four():
    g = Graph()
    g.add_node(0)
    g.add_node(1)
    g.add_node(2)
    g.add_node(3)
    g.add_edge((0, 1))
    g.add_edge((1, 2))
    g.add_edge((2, 3))
    g.add_edge((0, 3))
    eq_(len(g.find_a_maximum_matching()), 2)


def test_edmonds_five():
    g = Graph()
    nodes = [0, 1, 2, 3, 4, 5, 6, 7]
    edges = [(0, 3), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (4, 7)]
    for n in nodes:
        g.add_node(n)
    for e in edges:
        g.add_edge(e)
    eq_(len(g.find_a_maximum_matching()), 2)


def test_edmonds_six():
    g = Graph()
    nodes = [0, 1, 2, 3]
    edges = [(0, 1), (1, 2), (2, 3), (3, 1)]
    for n in nodes:
        g.add_node(n)
    for e in edges:
        g.add_edge(e)
    eq_(len(g.find_a_maximum_matching()), 2)


def test_edmonds_seven():
    g = Graph()
    nodes = [0, 1, 2, 3, 4]
    edges = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]
    for n in nodes:
        g.add_node(n)
    for e in edges:
        g.add_edge(e)
    eq_(len(g.find_a_maximum_matching()), 2)


def test_matchable_one():
    g = Graph()
    g.add_node(0)
    g.add_node(1)
    g.add_node(2)
    g.add_node(3)
    g.add_edge((0, 1))
    g.add_edge((1, 2))
    g.add_edge((2, 3))
    ll = g.find_maximally_matchable_edges()
    eq_(len(ll), 2)
    ok_((0,1) in ll)
    ok_((2,3) in ll)


def test_edmonds_seven():
    g = Graph()
    nodes = [0, 1, 2, 3]
    edges = [(0, 1), (1, 2), (2, 3), (3, 1)]
    for n in nodes:
        g.add_node(n)
    for e in edges:
        g.add_edge(e)
    ll = g.find_maximally_matchable_edges()
    eq_(len(ll), 2)
    ok_((0, 1) in ll)
    ok_((2, 3) in ll)
