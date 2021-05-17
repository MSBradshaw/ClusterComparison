import json
import os

# load each file, parse as json, add to dictionary indexed by gene symbol
all_genes_dict = {}
for item in os.listdir('HumanProteinAtlas/Data/'):
    with open('HumanProteinAtlas/Data/' + item,'r') as f:
        single_gene_info = json.load(f)
    all_genes_dict[single_gene_info['Gene']] = single_gene_info

# save the dictionary as one huge json file
with open('HumanProteinAtlas/all_genes_info.json', 'w') as json_file:
    json.dump(all_genes_dict, json_file)
