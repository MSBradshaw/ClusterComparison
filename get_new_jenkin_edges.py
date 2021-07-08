string_hpo_el = "Data/HPO_String_edgelist_june_22_2021.tsv"
old_jenkins_el = "Data/jenkins_12_2015.tsv"

old_jenkins = set()
old_jenkin_genes = set()
jenkins = set()
raw_jenkins = set()
raw_jenkin_genes = set()

for line in open(old_jenkins_el,'r'):
    row = line.strip().split('\t')
    old_jenkins.add(str(row))
    old_jenkin_genes.add(row[0])
print(row)

for line in open(string_hpo_el,'r'):
    row = line.strip().split('\t')
    row = row[:2]
    is_hpo = ['HP:' in x for x in row]
    if sum(is_hpo) == 1:
        hpo = row[[i for i,x in enumerate(is_hpo) if x][0]]
        try:
            gene = row[[i for i, x in enumerate(is_hpo) if not x][0]]
        except IndexError:
            continue
        if gene in old_jenkin_genes:
            jenkins.add(str([gene, hpo]))
        raw_jenkins.add(str([gene, hpo]))
        raw_jenkin_genes.add(gene)
print(row)


print(len(jenkins))
print(len(old_jenkins))
print(len(old_jenkin_genes))
print(len(raw_jenkins))
print(len(raw_jenkin_genes))

new_edges = [x for x in jenkins if x not in old_jenkins]

with open('Data/new_jenkins_edges.tsv','w') as outfile:
    for pair in new_edges:
        pair.replace("[",'').replace("]",'').replace("'",'').split(', ')
        outfile.write()
