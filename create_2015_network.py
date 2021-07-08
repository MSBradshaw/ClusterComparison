import networkx as nx
import obonet

# for the edges of HPO
G = obonet.read_obo('Data/hp_2015_week_46.obo')

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
    n2 = row[0]
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

print(len(not_founds))
print(len(founds))

# for the edges of Jenkins
for line in open('Data/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_12_2015.txt'):
    row = line.strip().split()
    if len(row) < 4: continue
    G.add_edge(row[1], row[3])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open('Data/String_HPO_2015.edgelist', 'w') as outfile:
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

# remove doubled edges

# remove circular edges