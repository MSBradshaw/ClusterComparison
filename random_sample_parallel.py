import networkx as nx
import argparse
import random
import pandas as pd
import multiprocessing

"""
python snowball.py --edgelist edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt
--output snowball.infomap.String_HPO_2015.phenotypic_branch.tsv
--coms Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt
--new_edges Data/new_jenkins_edges.tsv
--reps 100
"""


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='name of file to save the community to')

    parser.add_argument('--new_edges',
                        dest='new_edges',
                        required=True,
                        help='file with new edges, tab separated')

    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='file of communities and there members, each row is a com, '
                             'first item is the com name, all others are members')

    parser.add_argument('--reps',
                        dest='reps',
                        required=True,
                        default=100,
                        type=int,
                        help='number of repitions, default = 100')

    parser.add_argument('--threads',
                        dest='threads',
                        required=True,
                        default=1,
                        type=int,
                        help='number of threads to use, default = 1')

    args = parser.parse_args()
    return args


def score_com(com, edges_set):
    """

    :param com: list of members of a community
    :param edges_set: set of all edges (in both directions) named like "gene1_gene2" or "gene2_HP:0000001"
    :return: number of edges from the dict found in the com
    """
    possible_edges = [y + '_' + x for y in com for x in com]
    count = sum(e in edges_set for e in possible_edges)

    return count



def process(fun_args):
    (G, com, com_id, reps, new_edges) = fun_args
    data = {'com_id': [], 'com_score': [], 'replicate_id': [], 'replicate_score': [], 'rep_and_com_size': []}

    com_score = score_com(G, com, new_edges)
    for i in range(reps):
        # choose a random node from the community
        start = random.sample(list(G.nodes), 1)[0]
        # snowball till we have a set of size len(com)
        snowball = set()
        used = set()
        snowball.add(start)
        current = start
        finished = set()
        while len(snowball) < len(com):
            print('\t', str(len(snowball)))
            not_used = [x for x in snowball if x not in finished]
            try:
                current = random.sample(not_used, 1)[0]
            except ValueError:
                print(com)
                print(snowball)
                print(finished)
                return None
            neighs = list(G.neighbors(current))
            deficit = len(com) - len(snowball)
            choose_x = min(deficit, len(neighs))
            choices = random.sample(neighs, choose_x)
            if len(choices) == len(neighs):
                finished.add(current)
            for x in choices:
                snowball.add(x)
        snowball_score = score_com(G, list(snowball), new_edges)
        data['com_id'].append(com_id)
        data['com_score'].append(com_score)
        data['replicate_id'].append(i)
        data['replicate_score'].append(snowball_score)
        data['rep_and_com_size'].append(len(com))

    df = pd.DataFrame(data)
    return df


args = get_args()

_el = args.edgelist
_output = args.output
_coms = args.coms
_new_edge_list = args.new_edges
_reps = args.reps
_threads = args.threads

# _el = 'Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt'
# _output = 'del.txt'
# _coms = 'Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt'
# _new_edge_list = 'Data/new_jenkins_edges.tsv'
# _reps = 5
# _threads = 10

_new_edges = set()

for line in open(_new_edge_list, 'r'):
    row = line.strip().split('\t')
    _new_edges.add(row[0] + '_' + row[1])
    _new_edges.add(row[1] + '_' + row[0])

_G = nx.read_edgelist(_el)

# create the parameters for the multi-processing
param_sets = []
for line in open(_coms, 'r'):
    _row = line.strip().split('\t')
    _com = _row[1:]
    _com_id = _row[0]
    # this is one of the proteins that is in the weird connected component of size 2, it causes errors is not ignored
    if '9606.ENSP00000411694' in _com: continue
    param_sets.append([_G, _com, _com_id, _reps, _new_edges])

_p = multiprocessing.Pool(10)
_res = _p.map(process, param_sets)
_res_df = pd.concat(_res)

_res_df.to_csv(_output, index=False, sep='\t')
