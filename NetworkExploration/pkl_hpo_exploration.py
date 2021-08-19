import networkx as nx

P = nx.read_edgelist('Data/pkl.tsv')
hpos = [n for n in P.nodes if 'http://purl.obolibrary.org/obo/HP' in n]

nx.shortest_path(P, 'http://purl.obolibrary.org/obo/HP_0000831', 'http://purl.obolibrary.org/obo/HP_0000819')

spaths = []

for x in hpos:
    for y in hpos:
        if x == y:
            continue
        try:
            p = nx.shortest_path(P, x, y)
            spaths.append(p)
        except nx.exception.NetworkXNoPath:
            continue

intermediate_node_counts = {}
for p in spaths:
    for x in p[1:-1]:
        if x in intermediate_node_counts:
            intermediate_node_counts[x] += 1
        else:
            intermediate_node_counts[x] = 1

