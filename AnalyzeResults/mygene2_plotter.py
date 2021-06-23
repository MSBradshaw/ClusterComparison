import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
import numpy as np

results_files = ['cesna_jun21_mg2.tsv', 'greedy_jun21_mg2.tsv', 'info_jun21_mg2.tsv', 'walktrap_jun21_mg2.tsv']
df = None
for f in results_files:
    tdf = pd.read_csv(f, sep='\t')
    tdf['file'] = f
    if df is None:
        df = tdf
    else:
        df = pd.concat([df, tdf])

# add pseudo counts to avoid the 0's being long in a log scale
df['mg2_pairs_count'] += 1

"""
Index(['cluster_id', 'cluster_size', 'gene_ratio', 'HPO_ratio',
       'num_sig_go_enrichment_terms', 'sig_go_enrichment_p_vals',
       'sig_go_enrichment_fdr_corrected_p_vals', 'sig_go_enrichment_terms',
       'go_sig_threshold', 'max_norm_cell_type_specificity',
       'max_norm_cell_type_comma_sep_string', 'max_norm_disease_specificity',
       'max_norm_disease_comma_sep_string', 'mg2_pairs_count',
       'mg2_not_pairs_count', 'mg2_portion_families_recovered', 'file'],
      dtype='object')
"""

colors = ['#8931EF', '#F2CA19', '#FF00BD', '#0057E9', '#87E911', '#E11845']
files = list(set(df['file']))
files.sort()
markers = ['o', '^', 's', 'P', 'D', '*']

legend_elements = []
fig, axs = plt.subplots(2, 2, gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [3, 1]})
fig.set_size_inches(8, 8)
for i, f in enumerate(files):
    legend_elements.append(Line2D([0], [0],
                                  marker=markers[i],
                                  color=colors[i],
                                  label=f,
                                  markerfacecolor=colors[i],
                                  markersize=5,
                                  linewidth=0))
    sub = df[df['file'] == f]
    axs[1, 0].scatter(sub['mg2_pairs_count'],
               sub['mg2_portion_families_recovered'],
               color=colors[i],
               alpha=.5,
               marker=markers[i])

axs[0, 1].legend(handles=legend_elements,frameon=False)
vert_hist = np.histogram(df['mg2_portion_families_recovered'], bins=30)
axs[1, 1].plot(vert_hist[0], vert_hist[1][:-1])
hist = np.histogram(df['mg2_pairs_count'], bins=10000)
axs[0, 0].plot(hist[1][:-1], hist[0])

axs[1, 0].set_xscale('log')
axs[0, 0].set_xscale('log')
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
axs[1][0].set_xlabel('Num. MyGene2 Pairs in Com')
axs[1][0].set_ylabel('Portion of MyGene2\nPatients Rediscovered')
plt.savefig('Figures/mygene2_scatter.png', bbox_inches='tight', dpi=600)
plt.show()
