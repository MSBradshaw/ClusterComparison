from BOCC import BOCC
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

com_files = ['cesna_coms_june_22_2021.txt', 'greedy_coms_june_22_2021.txt', 'infomap_coms_june_22_2021.txt',
             'walktrap_coms_june_22_2021_filtered.txt']

c = BOCC.load_clusters(com_files[0])
t = c[0].mygene2_stats('Data/mygene2_gene_hpo_family.tsv')
fam = t[-1]

# x: total pairs
# y: num rediscovered
# color: aglo
total_num_pairs_per_fam = []
num_rediscovered_per_fam = []
algo = []
families = []
com_id = []
for file in com_files:
    coms = BOCC.load_clusters(file)
    for c in coms:
        t = c.mygene2_stats('Data/mygene2_gene_hpo_family.tsv')
        fam = t[-1]
        for key in fam.keys():
            x = fam[key]
            total_num_pairs_per_fam.append(len(x['pairs']) + len(x['not_pairs']))
            num_rediscovered_per_fam.append(len(x['pairs']))
            algo.append(file)
            families.append(key)
            com_id.append(c.name)

df = pd.DataFrame({'total_num_pairs': total_num_pairs_per_fam,
                   'rediscovered_pairs_count': num_rediscovered_per_fam,
                   'file': algo,
                   'family_id': families,
                   'com_id': com_id})

colors = ['#8931EF', '#F2CA19', '#FF00BD', '#0057E9', '#87E911', '#E11845']

legend_elements = []
fig, axs = plt.subplots(2, 2, gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [3, 1]})
fig.set_size_inches(8, 8)
for i, f in enumerate(com_files):
    legend_elements.append(Line2D([0], [0],
                                  marker='o',
                                  color=colors[i],
                                  label=f,
                                  markerfacecolor=colors[i],
                                  markersize=5,
                                  linewidth=0))
    sub = df[df['file'] == f]
    # add pseudo counts
    sub['rediscovered_pairs_count'] = sub['rediscovered_pairs_count'] + 1
    sub['total_num_pairs'] = sub['total_num_pairs'] + 1
    axs[1, 0].scatter(sub['total_num_pairs'],
                      sub['rediscovered_pairs_count'],
                      color=colors[i],
                      alpha=.5,
                      marker='o')

axs[0, 1].legend(handles=legend_elements, frameon=False)
vert_hist = np.histogram(df['total_num_pairs'], bins=1000)
axs[1, 1].plot(vert_hist[0], vert_hist[1][:-1])
hist = np.histogram(df['rediscovered_pairs_count'], bins=1000)
axs[0, 0].plot(hist[1][:-1], hist[0])

axs[1, 0].set_xscale('log')
axs[0, 0].set_xscale('log')
axs[1, 1].set_xscale('log')
axs[0, 0].set_yscale('log')
axs[0, 0].set_xlim(axs[1, 0].get_xlim())

for i in range(2):
    for j in range(2):
        axs[i][j].spines['top'].set_visible(False)
        axs[i][j].spines['right'].set_visible(False)
axs[0][1].spines['bottom'].set_visible(False)
axs[0][1].spines['left'].set_visible(False)
axs[0][1].set_yticks([])
axs[0][1].set_xticks([])
axs[0][0].set_xticks([])
axs[1][1].set_yticklabels([])
axs[1][0].set_xlabel('Number of pairs')
axs[1][0].set_ylabel('Number of rediscovered pairs')

flat_line_x = list(range(int(axs[1][0].get_xlim()[1])))
flat_line_y = [x for x in flat_line_x]

axs[1][0].plot(flat_line_x, flat_line_y, '--')

plt.savefig('Figures/mygene2_patient_level_scatter.png', bbox_inches='tight', dpi=600)
plt.show()

xx = df[df.total_num_pairs > 500]

# how many patients get rediscovered in a single community
# size of those distinct to a com vs not distinct
no_zeros = df[df['rediscovered_pairs_count'] != 0]
no_zeros['file_com_id'] = no_zeros['file'] + no_zeros['com_id'].astype(str)
g = no_zeros.groupby('family_id').agg(set)
gg = g.reset_index()
gg['num_coms'] = [len(x) for x in gg.file_com_id]
gg = gg.sort_values('num_coms')

# are there coms with only one patient?
b = no_zeros.groupby('file_com_id').agg(set)
bb = b.reset_index()
bb['num_patients'] = [len(x) for x in bb.family_id]
bb = bb.sort_values('num_patients')

