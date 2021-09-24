import networkx as nx
import obonet
import number_edge_list


"""
----------------------Make the Phenotypic Abnormality Only Edge List----------------------
"""

# for the edges of HPO
H = obonet.read_obo('Data/hp_July_8_2021.obo')

# make it a simple graph, not a multi or directed graph
h = nx.Graph()
for edge in H.edges:
    h.add_edge(edge[0], edge[1])

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
# ks = get_kids_that_pass_through_node(h, root, root)


for i in range(len(children)):
    if children[i] == pheno_ab: continue
    bad_kiddos = get_kids_that_pass_through_node(h, root, children[i])
    print([x for x in bad_kiddos if x in ks])

G = nx.Graph(h.subgraph(ks))
# use this line if you want to prune phenotypic abnormality from the HPO tree
# G.remove_node(pheno_ab)

# Rename to gene symbols
protein_mapping = {}
for line in open('Data/9606.protein.info.v11.0.txt'):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

# for the edges of String
error_count = 0
fine_count = 0
edges_added = 0
not_founds = set()
founds = set()
for line in open('Data/9606.protein.links.v11.0.txt'):
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
    edges_added += 1

print(len(not_founds))
print(len(founds))
print(edges_added)

# for the edges of Jenkins
for line in open('Data/genes_to_phenotype_july_8_2021.txt'):
    row = line.strip().split()
    if len(row) < 3: continue
    # this will skip adding edges that include edges pruned from HPO
    if row[2] not in G.nodes: continue
    G.add_edge(row[1], row[2])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open('Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt', 'w') as outfile:
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

number_edge_list.num_edgelist('Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt',
                              'Edgelists/String_HPO_2021.phenotypic_branch.numbered.edgelist.txt',
                              'Edgelists/String_HPO_2021.phenotypic_branch.nodenames.txt')


"""
----------------------Make the full HPO tree Edge List----------------------
"""

# for the edges of HPO
H = obonet.read_obo('Data/hp_July_8_2021.obo')

# make it a simple graph, not a multi or directed graph
G = nx.Graph()
for edge in H.edges:
    G.add_edge(edge[0], edge[1])

# Rename to gene symbols
protein_mapping = {}
for line in open('Data/9606.protein.info.v11.0.txt'):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

# for the edges of String
error_count = 0
fine_count = 0
edges_added = 0
not_founds = set()
founds = set()
for line in open('Data/9606.protein.links.v11.0.txt'):
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
    edges_added += 1

print(len(not_founds))
print(len(founds))
print(edges_added)

# for the edges of Jenkins
for line in open('Data/genes_to_phenotype_july_8_2021.txt'):
    row = line.strip().split()
    if len(row) < 3: continue
    G.add_edge(row[1], row[2])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open('Edgelists/String_HPO_2021.all_hpo.edgelist.txt', 'w') as outfile:
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

number_edge_list.num_edgelist('Edgelists/String_HPO_2021.all_hpo.edgelist.txt',
                              'Edgelists/String_HPO_2021.all_hpo.numbered.edgelist.txt',
                              'Edgelists/String_HPO_2021.all_hpo.nodenames.txt')
