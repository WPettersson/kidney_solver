import time

from nose.tools import raises, eq_, ok_

from kidney_solver.kidney_digraph import *
from kidney_solver.readers import KidneyReadException, read_digraph
from kidney_solver.agents import Donor
from tests.simple_find_cycles import simple_find_cycles

def read(filename):
    with open(filename) as f:
        lines = f.readlines()
    return read_digraph(lines)

def test_cycle_score():
    d = Digraph(4)
    d.add_edge(1.5, d.vs[0], d.vs[1], Donor(1, False, 25))
    d.add_edge(1, d.vs[1], d.vs[2], Donor(2, False, 25))
    d.add_edge(1, d.vs[2], d.vs[3], Donor(3, False, 25))
    d.add_edge(1, d.vs[3], d.vs[0], Donor(4, False, 25))
    c = [d.vs[i] for i in range(4)]
    eq_(cycle_score(c), 4.5)
    eq_(cycle_score(c, nhs=True), 16.696)

def test_shortest_path():
    d = read("test-fixtures/small1.input")
    spl = d.calculate_shortest_path_lengths(d.vs[1], 4)
    assert spl == [1, 0, 2, 2, 1, 2, 3, 4, 999999999]
    print(spl)

    min_vtx = 1
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: (e.tgt for e in v.edges if e.tgt.index() > min_vtx))
    assert spl == [999999999, 0, 999999999, 2, 1, 2, 3, 4, 999999999]
    print(spl)

    min_vtx = 7
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: (e.tgt for e in v.edges if e.tgt.index() > min_vtx))
    print(spl)
    assert spl == [999999999, 999999999, 999999999, 999999999, 999999999, 999999999, 999999999, 0, 1]

    min_vtx = 1
    spl = d.calculate_shortest_path_lengths(d.vs[min_vtx], 4,
            lambda v: [d.vs[v.index()+1]] if v.index()<4 else [])
    assert spl == [999999999, 0, 1, 2, 3, 999999999, 999999999, 999999999, 999999999]
    print(spl)


def test_find_cycles():
    d = read("test-fixtures/100.input")
    max_cycle = 5
    start = time.time()
    cycles = d.find_cycles(max_cycle)
    print(time.time() - start)
    print(len(cycles))
    start = time.time()
    slow_cycles = simple_find_cycles(d, max_cycle)
    print(time.time() - start)
    ok_(len(cycles) > 100)
    eq_(len(slow_cycles), len(cycles))

@raises(KidneyReadException)
def test_raises_exception_on_self_loop():
    d = read("test-fixtures/self-loop.input")

@raises(KidneyReadException)
def test_raises_exception_on_duplicate_edge():
    d = read("test-fixtures/duplicate-edge.input")

@raises(KidneyReadException)
def test_raises_exception_on_incorrect_edge_count_1():
    d = read("test-fixtures/incorrect_edge_count1.input")

@raises(KidneyReadException)
def test_raises_exception_on_incorrect_edge_count_2():
    d = read("test-fixtures/incorrect_edge_count2.input")

@raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_1():
    d = read("test-fixtures/out-of-range1.input")

@raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_2():
    d = read("test-fixtures/out-of-range2.input")

@raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_3():
    d = read("test-fixtures/out-of-range3.input")

@raises(KidneyReadException)
def test_raises_exception_on_index_out_of_range_4():
    d = read("test-fixtures/out-of-range4.input")
