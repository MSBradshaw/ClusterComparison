import sknetwork
from sknetwork.hierarchy import Paris, cut_straight, cut_balanced
import pandas as pd
import argparse
import networkx as nx
import random
import os


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='prefix to be use din output files')

    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='file of communities and there members, each row is a com, '
                             'first item is the com name, all others are members')

    args = parser.parse_args()
    return args

args = get_args()

# load the graph
el = args.edgelist
output = args.output
coms = args.coms

# el = 'Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt'
# output = 'DELTest/infomap.String_HPO_2015.phenotypic_branch'
# coms = 'Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt'

G = nx.read_edgelist(el)
for line in open(coms,'r'):
    row = line.strip().split('\t')
    id = row[0]
    com = row[1:]
    if len(com) < 3:
        print('Warning this community with < 3 members, skipping')
        continue
    g = nx.subgraph(G, com)
    tmp_el = '.hierarchical_clustering_temp_edge_list_' + str(int(random.random()  * 100000)) + '.txt'
    with open(tmp_el,'w') as outfile:
        for edge in g.edges():
            outfile.write('\t'.join(edge))
            outfile.write('\n')
    paris = Paris()
    adjacency = sknetwork.data.load_edge_list(tmp_el)
    dendrogram = paris.fit_transform(adjacency['adjacency'])
    # write the communities to a coms file
    num_coms = min(200, dendrogram.shape[0])
    if num_coms < 2:
        print('Warning: community has less than 2 members')
        continue
    com_ids = cut_balanced(dendrogram, num_coms)
    coms_df = pd.DataFrame({'node': list(adjacency['names']), 'com': list(com_ids)})

    with open(output + '_' + id + '.paris_balanced_cut.coms.txt', 'w') as outfile:
        for com in coms_df['com'].unique():
            sub = coms_df[coms_df['com'] == com]
            outfile.write('\t'.join([str(com)] + list(sub['node'])) + '\n')
    # delete the temp edge list
    os.remove(tmp_el)
