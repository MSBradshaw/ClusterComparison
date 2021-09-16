# import pandas as pd
# import networkx as nx
#
# G = nx.Graph()
#
# nodes_neighbors = {}
# # for the edges of Jenkins
# for line in open('Data/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_12_2015.txt'):
#     row = line.strip().split()
#     if len(row) < 4: continue
#     if (row[1], row[3]) in G.edges:
#         G[row[1]][row[3]]['weight'] = G[row[1]][row[3]]['weight'] + 1
#     else:
#         G.add_edge(row[1], row[3])
#         G[row[1]][row[3]]['weight'] = 1
#     if row[1] not in nodes_neighbors:
#         nodes_neighbors[row[1]] = []
#     if row[3] not in nodes_neighbors:
#         nodes_neighbors[row[3]] = []
#     nodes_neighbors[row[3]].append(row[1])
#     nodes_neighbors[row[1]].append(row[3])


import networkx as nx
import obonet
import number_edge_list

"""
----------------------Make the Phenotypic Abnormality Only Edge List----------------------
"""
# for the edges of HPO
H = obonet.read_obo('Data/hp_2015_week_46.obo')

# make it a simple graph, not a multi or directed graph
h = nx.Graph()
for edge in H.edges:
    h.add_edge(edge[0], edge[1])
    h[edge[0]][edge[1]]['weight'] = 1

# keep only thing that are children of Phenotypic abnormality
root = 'HP:0000001'
children = list(h.neighbors(root))
pheno_ab = 'HP:0000118'


def get_kids_that_pass_through_node(g, target, required_node):
    keepers = set()
    for n in g.nodes():
        path = nx.shortest_path(g, target, n)
        if required_node in path:
            for x in path:
                keepers.add(x)
    return keepers

# get the nodes that are part of the Phenotypic abnormality tree
ks = get_kids_that_pass_through_node(h, root, pheno_ab)
bad_kiddos0 = get_kids_that_pass_through_node(h, root, children[0])
bad_kiddos2 = get_kids_that_pass_through_node(h, root, children[2])
bad_kiddos3 = get_kids_that_pass_through_node(h, root, children[3])
print([x for x in bad_kiddos0 if x in ks])
print([x for x in bad_kiddos2 if x in ks])
print([x for x in bad_kiddos3 if x in ks])

# make a sub graph of just the Phenotypic abnormality tree
G = nx.Graph(h.subgraph(ks))


# Rename to gene symbols
protein_mapping = {}
for line in open('Data/9606.protein.info.v11.0.txt'):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

# for the edges of String
error_count = 0
fine_count = 0
not_founds = set()
founds = set()
for line in open('Data/9606.protein.links.v9.1.txt'):
    row = line.strip().split()
    n1 = row[0]
    n2 = row[1]
    try:
        n1 = protein_mapping[n1]
        founds.add(n1)
    except KeyError:
        not_founds.add(n1)

    try:
        n2 = protein_mapping[n2]
        founds.add(n2)
    except KeyError:
        not_founds.add(n2)

    G.add_edge(n1, n2)
    G[n1][n2]['weight'] = 1

print(len(not_founds))
print(len(founds))

for line in open('Data/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_12_2015.txt'):
    row = line.strip().split()
    if len(row) < 4: continue
    if (row[1], row[3]) in G.edges:
        G[row[1]][row[3]]['weight'] = G[row[1]][row[3]]['weight'] + 1
    else:
        G.add_edge(row[1], row[3])
        G[row[1]][row[3]]['weight'] = 1


# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open('Edgelists/String_HPO_2015.weighted.phenotypic_branch.edgelist.txt', 'w') as outfile:
    for edge in G.edges:
        row = [edge[0], edge[1]]
        row.sort()
        if row[0] == row[1]: continue
        if str(row) in edges:
            continue
        edges.add(str(row))
        outfile.write(row[0])
        outfile.write('\t')
        outfile.write(row[1])
        outfile.write('\t')
        outfile.write(str(G[row[0]][row[1]]['weight']))
        outfile.write('\n')

number_edge_list.num_edgelist('Edgelists/String_HPO_2015.weighted.phenotypic_branch.edgelist.txt',
                              'Edgelists/String_HPO_2015.weighted.phenotypic_branch.numbered.edgelist.txt',
                              'Edgelists/String_HPO_2015.weighted.phenotypic_branch.nodenames.txt')
