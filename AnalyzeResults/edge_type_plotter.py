import matplotlib.pyplot as plt
import networkx as nx
import plotly.express as px
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
import obonet
import networkx as nx


el = 'Data/String_HPO_2015.edgelist'

G = nx.read_edgelist(el, delimiter='\t')

# count number of each type of edge
# gene gene
gg_count = 0
gh_count = 0
hh_count = 0
for edge in G.edges:
    count = 0
    if 'HP:' in edge[0]:
        count += 1
    if 'HP:' in edge[1]:
        count += 1
    if count == 0:
        gg_count += 1
    elif count == 1:
        gh_count += 1
    elif count == 2:
        hh_count += 1

print(gg_count)
print(gh_count)
print(hh_count)

ggs = []
ghs = []
hhs = []
pure = []
not_pure = []
types = []
nodes = []
for node in G.nodes:
    n_gg = 0
    n_gh = 0
    n_hh = 0
    for n in nx.neighbors(G, node):
        count = 0
        if 'HP:' in node:
            count += 1
        if 'HP:' in n:
            count += 1
        if count == 0:
            n_gg += 1
        elif count == 1:
            n_gh += 1
        elif count == 2:
            n_hh += 1
    ggs.append(n_gg)
    ghs.append(n_gh)
    hhs.append(n_hh)
    pure.append(n_gg + n_hh)
    not_pure.append(n_gh)
    t = 'Gene'
    if 'HP:' in node:
        t = 'HPO'
    types.append(t)
    nodes.append(node)

ys = ggs + ghs + hhs
xs = ['Gene-Gene\nn=' + str(gg_count)] * len(ggs) + \
     ['Gene-HPO\nn=' + str(gh_count)] * len(ghs) + \
     ['HPO-HPO\nn=' + str(hh_count)] * len(hhs)
sns.boxplot(xs, ys)
plt.yscale('log')
plt.savefig('Figures/edge_type_box_plotter15.png')
plt.show()

pure_df = pd.DataFrame({'pure': pure, 'not_pure': not_pure, 'type': types, 'node': nodes})
color_map = {'HPO': 'red', 'Gene': 'blue'}
colors = [color_map[x] for x in types]
legend_elements=[]
for key in color_map.keys():
    legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=color_map[key],
                              label=key,
                              markerfacecolor=color_map[key],
                              markersize=5,
                              linewidth=0))

plt.scatter(pure, not_pure, c=colors)
plt.legend(handles=legend_elements,frameon=False)
plt.xlabel('Self-Self Type Edge Count')
plt.ylabel('Self-Other Type Edge Count')
plt.savefig('Figures/self-other-edge_count15.png')
plt.show()


plt.scatter(pure, not_pure, c=colors)
plt.legend(handles=legend_elements,frameon=False)
plt.xlabel('Self-Self Type Edge Count')
plt.ylabel('Self-Other Type Edge Count')
plt.xscale('log')
plt.yscale('log')
plt.savefig('Figures/self-other-edge_count_log15.png')
plt.show()


# find HPOs with not pure > 250
sub = pure_df[pure_df['type'] == 'HPO']
# how far are they from the root of the HPO?
sub['distance'] = [len(nx.shortest_path(G, 'HP:0000001', x)) if 'HP:' in x else -1 for x in sub['node']]

colors = {'Over 250': 'green', 'Under 250': 'blue'}
pc = ['green' if x else 'blue' for x in sub['not_pure'] > 250]

legend_elements=[]
for key in colors.keys():
    legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=colors[key],
                              label=key,
                              markerfacecolor=colors[key],
                              markersize=5,
                              linewidth=0))


plt.scatter(sub['distance'], sub['not_pure'], c=pc)
plt.ylabel('Self-Other Type Edge Count')
plt.xlabel('Distance from HPO Root')
plt.legend(handles=legend_elements,frameon=False)
plt.savefig('Figures/jenkins_edges15.png')
plt.show()

H = obonet.read_obo('Data/hp_July_8_2021.obo')
H.nodes['HP:0000028']
ss = sub[sub['not_pure'] > 250].sort_values('not_pure', ascending=False)
for x, n in zip(ss['node'],ss['not_pure']):
    print(x, H.nodes[x]['name'], str(n))


"""
-----------------------------------------------------------------------------------
Now do all of those plots again but for the 2021 graph
-----------------------------------------------------------------------------------
"""

el = 'Data/String_HPO_2021.edgelist'

G = nx.read_edgelist(el, delimiter='\t')

# count number of each type of edge
# gene gene
gg_count = 0
gh_count = 0
hh_count = 0
for edge in G.edges:
    count = 0
    if 'HP:' in edge[0]:
        count += 1
    if 'HP:' in edge[1]:
        count += 1
    if count == 0:
        gg_count += 1
    elif count == 1:
        gh_count += 1
    elif count == 2:
        hh_count += 1

print(gg_count)
print(gh_count)
print(hh_count)

ggs = []
ghs = []
hhs = []
pure = []
not_pure = []
types = []
nodes = []
for node in G.nodes:
    n_gg = 0
    n_gh = 0
    n_hh = 0
    for n in nx.neighbors(G, node):
        count = 0
        if 'HP:' in node:
            count += 1
        if 'HP:' in n:
            count += 1
        if count == 0:
            n_gg += 1
        elif count == 1:
            n_gh += 1
        elif count == 2:
            n_hh += 1
    ggs.append(n_gg)
    ghs.append(n_gh)
    hhs.append(n_hh)
    pure.append(n_gg + n_hh)
    not_pure.append(n_gh)
    t = 'Gene'
    if 'HP:' in node:
        t = 'HPO'
    types.append(t)
    nodes.append(node)

ys = ggs + ghs + hhs
xs = ['Gene-Gene\nn=' + str(gg_count)] * len(ggs) + \
     ['Gene-HPO\nn=' + str(gh_count)] * len(ghs) + \
     ['HPO-HPO\nn=' + str(hh_count)] * len(hhs)
sns.boxplot(xs, ys)
plt.yscale('log')
plt.savefig('Figures/edge_type_box_plotter21.png')
plt.show()

pure_df = pd.DataFrame({'pure': pure, 'not_pure': not_pure, 'type': types, 'node': nodes})
color_map = {'HPO': 'red', 'Gene': 'blue'}
colors = [color_map[x] for x in types]
legend_elements=[]
for key in color_map.keys():
    legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=color_map[key],
                              label=key,
                              markerfacecolor=color_map[key],
                              markersize=5,
                              linewidth=0))

plt.scatter(pure, not_pure, c=colors)
plt.legend(handles=legend_elements,frameon=False)
plt.xlabel('Self-Self Type Edge Count')
plt.ylabel('Self-Other Type Edge Count')
plt.savefig('Figures/self-other-edge_count21.png')
plt.show()


plt.scatter(pure, not_pure, c=colors)
plt.legend(handles=legend_elements,frameon=False)
plt.xlabel('Self-Self Type Edge Count')
plt.ylabel('Self-Other Type Edge Count')
plt.xscale('log')
plt.yscale('log')
plt.savefig('Figures/self-other-edge_count_log21.png')
plt.show()


# find HPOs with not pure > 250
sub = pure_df[pure_df['type'] == 'HPO']
# how far are they from the root of the HPO?
sub['distance'] = [len(nx.shortest_path(G, 'HP:0000001', x)) if 'HP:' in x else -1 for x in sub['node']]
sub['is_leaf'] = [ sum( 'HP:' in y for y in nx.neighbors(G,x)) == 1 for x in sub['node']]

colors = {'Over 250': 'green', 'Under 250': 'blue'}
pc = ['green' if x else 'blue' for x in sub['not_pure'] > 250]

legend_elements=[]
for key in colors.keys():
    legend_elements.append(Line2D([0], [0],
                              marker='o',
                              color=colors[key],
                              label=key,
                              markerfacecolor=colors[key],
                              markersize=5,
                              linewidth=0))


plt.scatter(sub['distance'], sub['not_pure'], c=pc)
plt.ylabel('Self-Other Type Edge Count')
plt.xlabel('Distance from HPO Root')
plt.legend(handles=legend_elements,frameon=False)
plt.savefig('Figures/jenkins_edges21.png')
plt.show()


colors = {True: 'red', False: 'blue'}
pc = ['red' if x else 'blue' for x in sub['is_leaf']]

legend_elements=[]
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
                              label='Not leaf',
                              markerfacecolor='blue',
                              markersize=5,
                              linewidth=0))


plt.scatter(sub['distance'], sub['not_pure'], c=pc)
plt.ylabel('Self-Other Type Edge Count')
plt.xlabel('Distance from HPO Root')
plt.legend(handles=legend_elements,frameon=False)
plt.savefig('Figures/jenkins_edges21_leaf.png')
plt.show()

H = obonet.read_obo('Data/hp_July_8_2021.obo')
H.nodes['HP:0000028']
ss = sub[sub['not_pure'] > 250].sort_values('not_pure', ascending=False)
a = ss[ss['is_leaf']]
sss = ss.iloc[0:28,:]

for x, n in zip(sss['node'],sss['not_pure']):
    print(x, H.nodes[x]['name'], str(n))


for x, n in zip(a['node'],a['not_pure']):
    print(x, H.nodes[x]['name'], str(n))



