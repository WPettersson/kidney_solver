"""Functions to read in data from input, usually stdin."""


from kidney_solver.kidney_digraph import Digraph
from kidney_solver.kidney_ndds import Ndd


class KidneyReadException(Exception):
    """Thrown on failure when reading from input."""
    pass


def read_digraph(lines):
    """Reads a digraph from an array of strings in the input format."""

    vtx_count, edge_count = [int(x) for x in lines[0].split()]
    digraph = Digraph(vtx_count)
    for line in lines[1:edge_count+1]:
        tokens = [x for x in line.split()]
        src_id = int(tokens[0])
        tgt_id = int(tokens[1])
        explanation = "Donor %s to patient %s" % (src_id, tgt_id)
        if src_id < 0 or src_id >= vtx_count:
            raise KidneyReadException("Vertex index {} out of range.".format(src_id))
        if tgt_id < 0 or tgt_id >= vtx_count:
            raise KidneyReadException("Vertex index {} out of range.".format(tgt_id))
        if src_id == tgt_id:
            raise KidneyReadException("Self-loop from {0} to {0} not permitted".format(src_id))
        if digraph.edge_exists(digraph.vs[src_id], digraph.vs[tgt_id]):
            raise KidneyReadException("Duplicate edge from {} to {}".format(src_id, tgt_id))
        score = float(tokens[2])
        digraph.add_edge(score, digraph.vs[src_id], digraph.vs[tgt_id], explanation)

    if lines[edge_count+1].split()[0] != "-1" or len(lines) < edge_count+2:
        raise KidneyReadException("Incorrect edge count")

    return digraph


def read_digraph_xml(lines):
    """Reads a digraph from a list of strings in JSON format. Note that this
    includes altruistic donors.
    """
    from defusedxml import ElementTree as ET
    data = ET.fromstring("\n".join(lines))
    num_altruists = 0
    num_regular = 0
    # The XML doesn't necessarily use sequential integers to number
    # donors/patients. In fact, it often skips numbers. These two dictionaries
    # are used to convert from a JSON-index, to our index.
    donors = data
    donor_lookup = {}
    patient_lookup = {}
    altruist_lookup = {}
    for index, donor in enumerate(donors):
        sources = donor.find("sources")
        if not sources:
            altruist_lookup[index] = num_altruists
            num_altruists += 1
        else:
            sources = list(sources)
            if len(sources) != 1:
                raise KidneyReadException("Only donors with exactly 1 source are supported")
            donor_lookup[index] = num_regular
            patient_lookup[int(sources[0].text)] = num_regular
            num_regular += 1
    digraph = Digraph(num_regular)
    ndds = [Ndd() for _ in range(num_altruists)]
    edge_exists = [[False for _ in digraph.vs] for __ in ndds]
    for index, donor in enumerate(donors):
        sources = donor.find("sources")
        if not sources:
            source = altruist_lookup[index]
            ndds[source].set_donor_id(donor.get("donor_id"))
            matches = donor.find("matches")
            if matches:
                for match in matches:
                    explanation = ("Donor %s (altruist) to patient %s" %
                                   (donor.get("donor_id"), match.find("recipient").text))
                    target = patient_lookup[int(match.find("recipient").text)]
                    score = float(match.find("score").text)
                    if source < 0 or source >= num_altruists:
                        raise KidneyReadException("Vertex index {} out of "
                                                  "range.".format(source))
                    if target < 0 or target >= num_regular:
                        raise KidneyReadException("Vertex index {} out of "
                                                  "range.".format(target))
                    if edge_exists[source][target]:
                        raise KidneyReadException("Duplicate edge from NDD {0} "
                                                  "to vertex {1}.".format(source, target))
                    digraph.vs[target].set_patient_id(match.find("recipient").text)
                    ndds[source].add_edge(digraph.vs[target], score, explanation)
                    edge_exists[source][target] = True
        else:
            sources = list(sources)
            if len(sources) != 1:
                raise KidneyReadException("Only donors with exactly 1 source are supported")
            source = patient_lookup[int(sources[0].text)]
            matches = donor.find("matches")
            if matches:
                for match in matches:
                    explanation = ("Donor %s to patient %s" %
                                   (index, match.find("recipient").text))
                    target = patient_lookup[int(match.find("recipient").text)]
                    score = float(match.find("score").text)
                    if source < 0 or source >= num_regular:
                        raise KidneyReadException("Vertex index {} out of "
                                                  "range.".format(source))
                    if target < 0 or target >= num_regular:
                        raise KidneyReadException("Vertex index {} out of "
                                                  "range.".format(target))
                    if source == target:
                        raise KidneyReadException("Self-loop from {0} to {0} not "
                                                  "permitted".format(source))
                    digraph.vs[source].set_donor_id(donor.attrib["donor_id"])
                    digraph.vs[target].set_patient_id(match.find("recipient").text)
                    if digraph.edge_exists(digraph.vs[source], digraph.vs[target]):
                        digraph.update_edge(digraph.vs[source],
                                            digraph.vs[target], score,
                                            explanation)
                    else:
                        digraph.add_edge(score, digraph.vs[source], digraph.vs[target], explanation)
    digraph.patient_lookup = patient_lookup
    digraph.altruist_lookup = altruist_lookup
    return digraph, ndds

def read_digraph_json(lines):
    """Reads a digraph from a list of strings in JSON format. Note that this
    includes altruistic donors.
    """
    import json
    json_data = json.loads("\n".join(lines))
    num_altruists = 0
    num_regular = 0
    # The JSON doesn't necessarily use sequential integers to number
    # donors/patients. In fact, it often skips numbers. These two dictionaries
    # are used to convert from a JSON-index, to our index.
    donor_lookup = {}
    patient_lookup = {}
    altruist_lookup = {}
    for index, donor in json_data["data"].items():
        if "altruistic" in donor:
            altruist_lookup[int(index)] = num_altruists
            num_altruists += 1
        else:
            donor_lookup[int(index)] = num_regular
            patient_lookup[donor["sources"][0]] = num_regular
            num_regular += 1
    digraph = Digraph(num_regular)
    ndds = [Ndd() for _ in range(num_altruists)]
    edge_exists = [[False for _ in digraph.vs] for __ in ndds]
    for index, donor in json_data["data"].items():
        if "altruistic" in donor:
            source = altruist_lookup[int(index)]
            ndds[source].set_donor_id(int(index))
            for match in donor["matches"]:
                explanation = ("Donor %s (altruist) to patient %s" %
                               (index, match["recipient"]))
                target = patient_lookup[match["recipient"]]
                score = match["score"]
                if source < 0 or source >= num_altruists:
                    raise KidneyReadException("Vertex index {} out of "
                                              "range.".format(source))
                if target < 0 or target >= num_regular:
                    raise KidneyReadException("Vertex index {} out of "
                                              "range.".format(target))
                if edge_exists[source][target]:
                    raise KidneyReadException("Duplicate edge from NDD {0} to "
                                              "vertex {1}.".format(source,
                                                                   target))
                digraph.vs[target].set_patient_id(match["recipient"])
                ndds[source].add_edge(digraph.vs[target], score, explanation)
                edge_exists[source][target] = True
        else:
            if len(donor["sources"]) != 1:
                raise KidneyReadException("We do not support donors providing "
                                          "for more than one patient")
            source = patient_lookup[donor["sources"][0]]
            for match in donor["matches"]:
                explanation = "Donor %s to patient %s" % (index,
                                                          match["recipient"])
                target = patient_lookup[match["recipient"]]
                score = match["score"]
                if source < 0 or source >= num_regular:
                    raise KidneyReadException("Vertex index {} out of "
                                              "range.".format(source))
                if target < 0 or target >= num_regular:
                    raise KidneyReadException("Vertex index {} out of "
                                              "range.".format(target))
                if source == target:
                    raise KidneyReadException("Self-loop from {0} to {0} not "
                                              "permitted".format(source))
                digraph.vs[source].set_donor_id(donor["sources"][0])
                digraph.vs[target].set_patient_id(match["recipient"])
                if digraph.edge_exists(digraph.vs[source], digraph.vs[target]):
                    digraph.update_edge(digraph.vs[source], digraph.vs[target], score, explanation)
                else:
                    digraph.add_edge(score, digraph.vs[source], digraph.vs[target], explanation)
    digraph.patient_lookup = patient_lookup
    digraph.altruist_lookup = altruist_lookup
    return digraph, ndds


def read_ndds(lines, digraph):
    """Reads NDDs from an array of strings in the .ndd format."""

    ndds = []
    ndd_count, edge_count = [int(x) for x in lines[0].split()]
    ndds = [Ndd() for _ in range(ndd_count)]

    # Keep track of which edges have been created already so that we can
    # detect duplicates
    edge_exists = [[False for _ in digraph.vs] for __ in ndds]

    for line in lines[1:edge_count+1]:
        tokens = [t for t in line.split()]
        src_id = int(tokens[0])
        tgt_id = int(tokens[1])
        score = float(tokens[2])
        explanation = "Donor %s (altruist) to patient %s" % (src_id, tgt_id)
        if src_id < 0 or src_id >= ndd_count:
            raise KidneyReadException("NDD index {} out of range.".format(src_id))
        if tgt_id < 0 or tgt_id >= digraph.n:
            raise KidneyReadException("Vertex index {} out of range.".format(tgt_id))
        if edge_exists[src_id][tgt_id]:
            raise KidneyReadException("Duplicate edge from NDD {0} to vertex "
                                      "{1}.".format(src_id, tgt_id))
        ndds[src_id].add_edge(digraph.vs[tgt_id], score, explanation)
        edge_exists[src_id][tgt_id] = True

    if lines[edge_count+1].split()[0] != "-1" or len(lines) < edge_count+2:
        raise KidneyReadException("Incorrect edge count")

    return ndds
