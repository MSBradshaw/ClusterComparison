import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

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

# how many components does this split HPO into?
print('Number of connected components with STRING', str(nx.algorithms.number_connected_components(G)))

# sizes of these connected components
cc = list(nx.algorithms.connected_components(G))
cc_sizes = [len(x) for x in cc]
sns.histplot(cc_sizes)
plt.xlabel('Size of connected components')
plt.savefig('Figures/prune_2015_hist.png')
plt.show()

# Rename to gene symbols
protein_mapping = {}
for line in open('Data/9606.protein.info.v11.0.txt'):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

g = nx.Graph()
for line in open('Data/9606.protein.links.v9.1.txt'):
    if 'protein1' in line: continue
    row = line.strip().split()
    n1 = row[0]
    n2 = row[1]
    try:
        n1 = protein_mapping[n1]
    except KeyError:
        pass

    try:
        n2 = protein_mapping[n2]
    except KeyError:
        pass

    g.add_edge(n1, n2)

# how many of the coms of size 1 have no connections in STRING 2015
[len(list(g.neighbors(list(x)[0]))) for x in cc if len(x) == 1 and list(x)[0] in g]

# The 46 coms of size one are not in STRING 2015 at all, are they in Jenkins?

# load Jenkins

jenkins_2015_genes = set()
jenkins_2015_hpos = set()
for line in open('Data/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_12_2015.txt', 'r'):
    if '#' == line[0]: continue
    row = line.strip().split('\t')
    jenkins_2015_genes.add(row[1])
    jenkins_2015_hpos.add(row[3])

print('Number the 46 founds in Jenkins 2015',
      str(len([x for x in cc if len(x) == 1 and list(x)[0] in jenkins_2015_genes])))
print('Number of HPO nodes post-pruning', str(len([n for n in G.nodes if 'HP:' in n])))
print('Number of HPO edges post-pruning', str(len([e for e in G.edges if 'HP:' in e[0] and 'HP:' in e[1]])))
"""
There is 1 community with 25748 members, 2 communities with 2 members and 46 with 1.
The 2 coms of size 2 are the same 2 that existed in the unpruned version
The 46 of size 1 are all genes that have no connections in STRING 2015
"""

# what is the connectivity of just HPO post pruning
# get the names of the HPO nodes
hpos = [n for n in G.nodes if 'HP:' in n]

sub = nx.subgraph(G, hpos)
post_cc = list(nx.algorithms.connected_components(sub))
print(post_cc)
print('Number of connected components pre-pruning', str(len(post_cc)))

post_lens = [len(x) for x in post_cc]

sns.histplot(post_lens)
plt.xlabel('Size of connected components')
plt.savefig('Figures/prune_2015_HPO_ONLY_hist.png')
plt.show()

# What are these 47? Leafs? What is there depth?

singles = [list(x)[0] for x in post_cc if len(x) == 1]
are_singles_leafs = [leaf_node_look_up[x] for x in singles]

# are these disconnected from HPOs in the original graph?
# how many of these are fully disconnected from other HPOs in the un pruned graph?
sum(sum('HP:' in y for y in og_G.neighbors(x)) == 0 for x in singles)

# Where did these 46 come from?
# are they in the raw HPO 2015 data?
print('Number of the disconnected HPOs in the original 2015 HPO network', str(sum(x in g.nodes for x in singles)))
# are they from Jenkins?
print('Number of the disconnected HPOs in the original 2015 Jenkins list',
      str(sum(x in jenkins_2015_hpos for x in singles)))
# they are all from Jenkins

# are the 46 disconnected Genes and 47 disconnected HPOs related?
single_genes = [list(x)[0] for x in cc if len(x) == 1 and 'protein' not in list(x)[0]]
# create a dictionary to indicate if a single gene is connected any single hpos
single_genes_connected_to_a_single_hpo = {x: 0 for x in single_genes}
for hpo in singles:
    hpo_neighbors = og_G.neighbors(hpo)
    for gene in single_genes:
        if gene in hpo_neighbors:
            print(gene)
            print(hpo)
            single_genes_connected_to_a_single_hpo[gene] += 1

# The only connected pair is connected to Autosomal recessive inheritance, this node should not be in this graph

# Even though HPO was pruned, Jenkins was not. The non-phenotypic abnormality branch terms need to be removed.

