import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

G = nx.read_edgelist('Data/HPO_String_edgelist.tsv')

gene_gene_deg = []
hpo_gene_deg = []
hpo_hpo_deg = []
genes_over_5000 = []
for n in G.nodes:
    is_hpo = 'HP:' in n
    # count HPO neighbors
    hpo_counts = 0
    # count gene neighbors
    gene_counts = 0
    for x in nx.neighbors(G, n):
        if 'HP:' in x:
            hpo_counts += 1
        else:
            gene_counts += 1
    if is_hpo:
        hpo_gene_deg.append(gene_counts)
        hpo_hpo_deg.append(hpo_counts)
    else:
        if gene_counts > 5000:
            genes_over_5000.append(n)
        gene_gene_deg.append(gene_counts)
        hpo_gene_deg.append(hpo_counts)

sns.boxplot(y=hpo_hpo_deg + gene_gene_deg + hpo_gene_deg, x=(['Gene-Gene'] * len(gene_gene_deg)) +
                                                            (['HPO-HPO'] * len(hpo_hpo_deg)) +
                                                            (['Gene-HPO'] * len(hpo_gene_deg)))
plt.yscale('log')
plt.show()

# TODO what are the genes with >5000 other gene connections
print(genes_over_5000)
# TODO what is the total number of Gene-HPO connections
# divising by two becuase of double counting
print(sum(hpo_gene_deg) / 2)
# TODO same plots but for PKL, or merge them with PKL

skips = ['https://reactome.org/content/detail',
         'http://purl.bioontology.org/ontology/SNOMEDCT',
         'http://purl.obolibrary.org/obo/HsapDv',
         'http://purl.obolibrary.org/obo/UBERON',
         'http://purl.obolibrary.org/obo/CLO',
         'http://purl.obolibrary.org/obo/CARO',
         'http://purl.obolibrary.org/obo/OAE',
         'http://purl.obolibrary.org/obo/NCBITaxon',
         'http://purl.obolibrary.org/obo/PR',
         'http://purl.obolibrary.org/obo/GO',
         'http://purl.obolibrary.org/obo/IAO',
         'http://purl.obolibrary.org/obo/ENVO',
         'http://purl.obolibrary.org/obo/ECO',
         'http://purl.obolibrary.org/obo/CL',
         'http://purl.obolibrary.org/obo/DOID',
         'http://purl.obolibrary.org/obo/MOD',
         'http://purl.obolibrary.org/obo/HGNC',
         'http://purl.obolibrary.org/obo/VO',
         'http://purl.obolibrary.org/obo/TRANS',
         'http://purl.obolibrary.org/obo/FBbt',
         'http://purl.obolibrary.org/obo/GENO',
         'http://purl.obolibrary.org/obo/NBO',
         'http://purl.obolibrary.org/obo/OBI',
         'http://purl.obolibrary.org/obo/MPATH',
         'http://purl.obolibrary.org/obo/FOODON',
         'http://purl.obolibrary.org/obo/CHEBI',
         'http://purl.obolibrary.org/obo/RO',
         'http://purl.obolibrary.org/obo/D96882F1-8709-49AB-BCA9-772A67EA6C33',
         'http://purl.obolibrary.org/obo/BFO',
         'http://purl.obolibrary.org/obo/OGMS',
         'http://purl.obolibrary.org/obo/FMA',
         'http://purl.obolibrary.org/obo/PW',
         'http://purl.obolibrary.org/obo/SO',
         'http://purl.obolibrary.org/obo/OMIM',
         'http://purl.obolibrary.org/obo/CP',
         'http://purl.obolibrary.org/obo/OGG',
         'http://purl.obolibrary.org/obo/NCIT']
hpo_like = ['http://purl.obolibrary.org/obo/SYMP',
            'http://purl.obolibrary.org/obo/PATO',
            'http://purl.obolibrary.org/obo/UPHENO',
            'http://purl.obolibrary.org/obo/HP']
gene_tags = ['https://uswest.ensembl.org/Homo_sapiens/Transcript', 'https://www.ncbi.nlm.nih.gov/gene',
             'http://zfin.org/action/marker/view', 'https://uswest.ensembl.org/Homo_sapiens/Transcript',
             'http://www.informatics.jax.org/marker', 'https://www.ncbi.nlm.nih.gov/gene']

P = nx.read_edgelist('Data/pkl.tsv')

pkl = nx.Graph()
for line in open('Data/pkl_raw.tsv', 'r'):
    row = line.strip().split('\t')
    # skip all the non-wanted info
    # if sum(1 for x in skips if x in line) > 0:
    #     continue
    pkl.add_edge(row[0], row[1])
    if sum(1 for x in hpo_like if x in row[0]) > 0:
        pkl.nodes[row[0]]['type'] = 'HPO_like'
    elif sum(1 for x in gene_tags if x in row[0]) > 0:
        pkl.nodes[row[0]]['type'] = 'gene'
    else:
        pkl.nodes[row[0]]['type'] = 'none'

    if sum(1 for x in hpo_like if x in row[1]) > 0:
        pkl.nodes[row[1]]['type'] = 'HPO_like'
    elif sum(1 for x in gene_tags if x in row[1]) > 0:
        pkl.nodes[row[1]]['type'] = 'gene'
    else:
        pkl.nodes[row[1]]['type'] = 'none'

pkl_gene_gene_deg = []
pkl_hpo_gene_deg = []
pkl_hpo_hpo_deg = []

for n in pkl.nodes:
    is_hpo = pkl.nodes[n]['type'] == 'HPO_like'
    # count HPO neighbors
    pkl_hpo_counts = 0
    # count gene neighbors
    pkl_gene_counts = 0
    for x in nx.neighbors(pkl, n):
        if pkl.nodes[x]['type'] == 'HPO_like':
            pkl_hpo_counts += 1
        else:
            pkl_gene_counts += 1
    if is_hpo:
        pkl_hpo_gene_deg.append(pkl_gene_counts)
        pkl_hpo_hpo_deg.append(pkl_hpo_counts)
    else:
        pkl_gene_gene_deg.append(pkl_gene_counts)
        pkl_hpo_gene_deg.append(pkl_hpo_counts)

sns.boxplot(y=hpo_hpo_deg + gene_gene_deg + hpo_gene_deg +
              pkl_gene_gene_deg + pkl_hpo_hpo_deg + pkl_hpo_gene_deg,
            x=(['Gene-Gene'] * len(gene_gene_deg)) +
              (['HPO-HPO'] * len(hpo_hpo_deg)) +
              (['Gene-HPO'] * len(hpo_gene_deg)) +
              (['PKL Gene-Gene'] * len(pkl_gene_gene_deg)) +
              (['PKL HPO-HPO'] * len(pkl_hpo_hpo_deg)) +
              (['PKL Gene-HPO'] * len(pkl_hpo_gene_deg))
            )
plt.yscale('log')
plt.xticks([0, 1, 2, 3, 4, 5],
           ['Gene-Gene\nn='+str(len(gene_gene_deg)),
            'HPO-HPO\nn='+str(len(hpo_hpo_deg)),
            'Gene-HPO\nn='+str(int(len(hpo_gene_deg)/2)),
            'PKL Gene-Gene\nn='+str(len(pkl_gene_gene_deg)),
            'PKL HPO-HPO\nn='+str(len(pkl_hpo_hpo_deg)),
            'PKL Gene-HPO\nn='+str(int(len(pkl_hpo_gene_deg)/2))],
           rotation=45)
plt.ylabel('Degree')
plt.tight_layout()
plt.savefig('Figures/pkl_stringHPO_degrees.boxplot.png',bbox_inches='tight')
plt.show()

# check if the HPO terms in B really are not connected to each other
