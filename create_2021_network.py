import networkx as nx
import obonet

# for the edges of HPO
G = obonet.read_obo('Data/hp_July_8_2021.obo')

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
for line in open('Data/9606.protein.links.v11.0.txt'):
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
for line in open('Data/genes_to_phenotype_july_8_2021.txt'):
    row = line.strip().split()
    if len(row) < 3: continue
    G.add_edge(row[1], row[2])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open('Data/String_HPO_2021.edgelist', 'w') as outfile:
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
