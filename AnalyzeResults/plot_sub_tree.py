import matplotlib.pyplot as plt
import networkx as nx
import obonet
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--all_coms',
                        dest='all_coms',
                        required=True,
                        help='file with all communities listed, first column is algorithm type, second is community id')

    parser.add_argument('--algo',
                        dest='algo',
                        required=True,
                        help='name of algorithm of interest')

    parser.add_argument('--com',
                        dest='com',
                        required=True,
                        type=int,
                        help='community id of interest')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='community id of interest')

    args = parser.parse_args()
    return args

args = get_args()
el = args.edgelist
all_coms = args.all_coms
algo = args.algo
com_num = args.com
output = args.output

com = None

for line in open(all_coms,'r'):
    if algo not in line:
        continue
    row = line.strip().split('\t')
    if int(row[1]) != com_num:
        continue
    com = row[2:]

# el = 'Data/HPO_String_edgelist_june_22_2021.tsv'
# all_coms = 'all_june_22_coms.txt'
# algo = 'infomap'
# com_num = 5
# output = 'filename2.png'

# el = 'Data/HPO_String_edgelist_june_22_2021.tsv'
url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
hpo_g = obonet.read_obo(url)
G = nx.read_edgelist(el)
# hpos = ['HP:0012093', 'HP:0002012', 'HP:0001732', 'HP:0011620', 'HP:0012091']
# genes = ['CFAP53']
hpos = []
genes = []
for mem in com:
    if 'HP:' in mem:
        hpos.append(mem)
    else:
        genes.append(mem)

print(hpos)
print(genes)
distance_2_root = []
nodes = set()
# find the top of the tree

for h in hpos:
    try:
        path = nx.shortest_path(hpo_g,h, 'HP:0000001')
    except nx.exception.NodeNotFound:
        continue
    for p in path:
        nodes.add(p)
    for gene in genes:
        try:
            path = nx.shortest_path(G, h, gene)
        except nx.exception.NodeNotFound:
            continue
        for p in path:
            nodes.add(p)
    distance_2_root.append(len(path))


distance_2_root, hpos = zip(*sorted(zip(distance_2_root, hpos)))

g = nx.subgraph(G, list(nodes) + genes)

pos = nx.nx_agraph.graphviz_layout(g, prog="twopi")

# nx.draw(g, pos)
in_com = list(hpos)+genes
others = [n for n in nodes if n not in in_com]
fig = plt.gcf()
fig.set_size_inches(10, 10)
nx.draw_networkx_nodes(g,pos,nodelist=in_com,node_color='green')
nx.draw_networkx_nodes(g,pos,nodelist=others,node_color='blue')
nx.draw_networkx_edges(g, pos, width=1.0, alpha=0.5)
labels = {n:n+'\n'+hpo_g.nodes[n]['name'] for n in nodes if 'HP:' in n}
for gene in nodes:
    if 'HP:' not in gene:
        labels[gene] = gene
nx.draw_networkx_labels(g,pos,labels,font_size=10,font_color='black')
plt.savefig(args.output)

# plt.show()


