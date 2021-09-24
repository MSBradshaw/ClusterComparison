import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import number_edge_list


G = nx.read_edgelist('Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt')
og_G = nx.read_edgelist('Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt')
og_cc = list(nx.algorithms.connected_components(G))
print(og_cc)
print('Number of connected components pre-pruning', str(len(og_cc)))
print('Number of HPO nodes pre-pruning', str(len([n for n in G.nodes if 'HP:' in n])))
print('Number of HPO edges pre-pruning', str(len([e for e in G.edges if 'HP:' in e[0] and 'HP:' in e[1]])))

"""
There are 3 connected components, one with virtually everything, 
one with two column labels that accidentally got added as nodes and one with two proteins that appear to not exist in 
STRING anymore and are reported as having no interactions on other sources.
"""

leaf_node_look_up = {}
distance_look_up = {}

for n in G.nodes:
    if 'HP:' not in n: continue
    neighs = nx.neighbors(G, n)
    try:
        distance_look_up[n] = len(nx.shortest_path(G, n, 'HP:0000118'))
    except nx.exception.NetworkXNoPath:
        distance_look_up[n] = 0

for n in G.nodes:
    if 'HP:' not in n: continue
    leaf_node_look_up[n] = len(
        [x for x in neighs if x in distance_look_up and distance_look_up[x] > distance_look_up[n]]) == 0

nodes_to_prune = []

for n in G.nodes:
    if 'HP:' not in n: continue
    neighs = nx.neighbors(G, n)
    # prune the node if non of it's children are leafs
    # get children
    children = [x for x in neighs if x in distance_look_up and distance_look_up[x] > distance_look_up[n]]
    # get leaf children
    leaf_children = [x for x in children if leaf_node_look_up[x]]
    if len(leaf_children) == 0:
        nodes_to_prune.append(n)

# number of nodes that will be removed from the HPO tree
print(len(nodes_to_prune))

for n in nodes_to_prune:
    G.remove_node(n)

for x in ['HP:0001470', 'HP:0100613', 'HP:0010984', 'HP:0001427', 'HP:0003811', 'HP:0003831', 'HP:0003743', 'HP:0003678', 'HP:0003584', 'HP:0001452', 'HP:0003676', 'HP:0003677', 'HP:0001450', 'HP:0001423', 'HP:0001466', 'HP:0001417', 'HP:0001522', 'HP:0012275', 'HP:0003829', 'HP:0003812', 'HP:0011463', 'HP:0003819', 'HP:0001699', 'HP:0011462', 'HP:0003674', 'HP:0003587', 'HP:0003581', 'HP:0001442', 'HP:0003680', 'HP:0001425', 'HP:0003745', 'HP:0003828', 'HP:0003623', 'HP:0003826', 'HP:0003744', 'HP:0003593', 'HP:0001426', 'HP:0001472', 'HP:0001419', 'HP:0010982', 'HP:0003577', 'HP:0001428', 'HP:0001444', 'HP:0003621', 'HP:0000007', 'HP:0000006', 'HP:0003596']:
    if x in G.nodes:
        print('FUCK', x)


edges = set()
with open('Edgelists/String_HPO_2015.pruned.edgelist.txt', 'w') as outfile:
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
        outfile.write('\n')

number_edge_list.num_edgelist('Edgelists/String_HPO_2015.pruned.edgelist.txt',
                              'Edgelists/String_HPO_2015.pruned.numbered.edgelist.txt',
                              'Edgelists/String_HPO_2015.pruned.nodenames.txt')
