import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import obonet
from matplotlib.lines import Line2D
import seaborn as sns

H = obonet.read_obo('Data/hp_July_8_2021.obo')
# make it a simple graph, not a multi or directed graph
h = nx.Graph()
for edge in H.edges:
    h.add_edge(edge[0], edge[1])

hpo_to_genes = {}

hpos = []
distances = []
genes = []
file = []
is_leaf = []
root = 'HP:0000001'

hpo_to_genes = {}
for line in open('Data/phenotype_to_genes.txt'):
    if '<tab>' in line: continue
    row = line.strip().split('\t')
    hp = row[0]
    distances.append(len(nx.shortest_path(h, root, hp)))
    file.append('P_to_G')
    hpos.append(hp)
    is_leaf.append(len(list(h.neighbors(hp))) == 1)
    if hp not in hpo_to_genes:
        hpo_to_genes[hp] = set()
    hpo_to_genes[hp].add(row[3])
genes = genes + [len(hpo_to_genes[hp]) for hp in hpos]

hpo_to_genes = {}
hpos = []
for line in open('Data/genes_to_phenotype.txt'):
    if '<tab>' in line: continue
    row = line.strip().split('\t')
    hp = row[2]
    distances.append(len(nx.shortest_path(h, root, hp)))
    file.append('G_to_P')
    hpos.append(hp)
    is_leaf.append(len(list(h.neighbors(hp))) == 1)
    if hp not in hpo_to_genes:
        hpo_to_genes[hp] = set()
    hpo_to_genes[hp].add(row[1])
genes = genes + [len(hpo_to_genes[hp]) for hp in hpos]

colors = {'G_to_P': 'red', 'P_to_G': 'blue'}

pc = [colors[x] for x in file]

legend_elements = []
for key in colors.keys():
    legend_elements.append(Line2D([0], [0],
                                  marker='o',
                                  color=colors[key],
                                  label=key,
                                  markerfacecolor=colors[key],
                                  markersize=5,
                                  linewidth=0))

plt.scatter(x=distances, y=genes, c=pc)
plt.xlabel('Distance from root')
plt.ylabel('Number of genes')
plt.legend(handles=legend_elements, frameon=False)
plt.savefig('Figures/annotation_connectivity.png')
plt.show()

# how many leaf nodes are there at each depth
df = pd.DataFrame({'is_leaf': is_leaf, 'distance': distances})

sns.histplot(data=df, x='distance', hue=is_leaf, multiple="dodge", shrink=.8)
x1 = df[df['is_leaf']]['distance']
x2 = df[~df['is_leaf']]['distance']
plt.hist([x1, x2])
plt.xlabel('Distance from HPO root')
plt.ylabel('Freq')
plt.title('Is the node a leaf (blue = yes)')
plt.savefig('leaf_hist.png')
plt.show()

df = pd.DataFrame({'is_leaf': is_leaf, 'distance': distances, 'num_genes': genes, 'file': file})
df = df[df['file'] == 'G_to_P']
colors = {True: 'red', False: 'blue'}

pc = [colors[x] for x in df['is_leaf']]

legend_elements = []

legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color='red',
                              label='Leaf',
                              markerfacecolor='red',
                              markersize=5,
                              linewidth=0))
legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color='blue',
                              label='Not Leaf',
                              markerfacecolor='blue',
                              markersize=5,
                              linewidth=0))

df['adj_dist'] = [ df.iloc[i,1] + .5 if df.iloc[i,0] else df.iloc[i,1] for i in range(df.shape[0])]

plt.scatter(x=df['adj_dist'], y=df['num_genes'], c=pc)
plt.xlabel('Distance from root')
plt.ylabel('Number of genes')
plt.legend(handles=legend_elements, frameon=False)
plt.savefig('Figures/annotation_connectivity_leaf.png')
plt.show()
