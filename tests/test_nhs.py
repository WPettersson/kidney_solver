import logging
from nose.tools import eq_, assert_almost_equal
from kidney_solver.readers import read_digraph_file
import kidney_solver.kidney_ip as k_ip

LOGGER = logging.getLogger(__name__)

def test_instance_with_twocycle_and_threecycle():
    graph, altruists = read_digraph_file("test-fixtures/test1.json")
    opt_result = k_ip.optimise_picef_nhs(k_ip.OptConfig(graph, ndds=altruists,
                                                        max_cycle=3,
                                                        max_chain=2))
    eq_(len(opt_result.cycles), 1)
    eq_(opt_result.total_score, 8.098)


def test_instance_with_twocycle_and_threecycle_max_matchable():
    graph, altruists = read_digraph_file("test-fixtures/test1.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    cfg._constrain_maximal = True
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 1)
    eq_(opt_result.total_score, 8.098)


def test_instance_with_two_threeway_or_three_twoway():
    graph, altruists = read_digraph_file("test-fixtures/test2.json")
    opt_result = k_ip.optimise_picef_nhs(k_ip.OptConfig(graph, ndds=altruists,
                                                        max_cycle=3,
                                                        max_chain=2))
    eq_(len(opt_result.cycles), 3)
    assert_almost_equal(opt_result.total_score, 24.294)


def test_instance_with_two_threeway_or_three_twoway_max_matchable():
    graph, altruists = read_digraph_file("test-fixtures/test2.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    cfg._constrain_maximal = True
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 3)
    assert_almost_equal(opt_result.total_score, 24.294)


def test_instance_with_threeway_with_backarc():
    graph, altruists = read_digraph_file("test-fixtures/test3.json")
    opt_result = k_ip.optimise_picef_nhs(k_ip.OptConfig(graph, ndds=altruists,
                                                        max_cycle=3,
                                                        max_chain=2))
    eq_(len(opt_result.cycles), 1)
    assert_almost_equal(opt_result.total_score, 12.147)


def test_instance_with_threeway_with_backarc_max_matchable():
    graph, altruists = read_digraph_file("test-fixtures/test3.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    cfg._constrain_maximal = True
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 1)
    assert_almost_equal(opt_result.total_score, 12.147)


def test_instance_with_threeway_with_backarc_or_threeway_without():
    graph, altruists = read_digraph_file("test-fixtures/test4.json")
    opt_result = k_ip.optimise_picef_nhs(k_ip.OptConfig(graph, ndds=altruists,
                                                        max_cycle=3,
                                                        max_chain=2))
    eq_(len(opt_result.cycles), 1)
    assert_almost_equal(opt_result.total_score, 12.147)
    cycle = [v.patient_id() for v in opt_result.cycles[0]]
    eq_(cycle, [1, 2, 4])


def test_instance_with_threeway_with_backarc_or_threeway_without_max_matchable():
    graph, altruists = read_digraph_file("test-fixtures/test4.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    cfg._constrain_maximal = True
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 1)
    assert_almost_equal(opt_result.total_score, 12.147)
    cycle = [v.patient_id() for v in opt_result.cycles[0]]
    eq_(cycle, [1, 2, 4])


def test_instance_with_one_chain():
    graph, altruists = read_digraph_file("test-fixtures/test-alt1.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 0)
    eq_(len(opt_result.chains), 1)
    assert_almost_equal(opt_result.total_score, 8.098)


def test_instance_with_chain_with_option():
    graph, altruists = read_digraph_file("test-fixtures/test-alt2.json")
    cfg = k_ip.OptConfig(graph, ndds=altruists, max_cycle=3, max_chain=2)
    opt_result = k_ip.optimise_picef_nhs(cfg)
    eq_(len(opt_result.cycles), 0)
    eq_(len(opt_result.chains), 1)
    assert_almost_equal(opt_result.total_score, 9.098)
