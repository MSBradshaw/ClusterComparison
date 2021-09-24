import networkx as nx

G = nx.read_edgelist('Data/HPO_String_edgelist_multi_edge.tsv')

nx.write_edgelist(G,'Data/HPO_String_edgelist_june_22_2021.tsv')

