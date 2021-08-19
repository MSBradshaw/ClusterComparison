import sys
import obonet

"""
Usage:
as stdin give this script as many tab delimited HPOs as you like and it will print out the name and description

example:
cat all_june_22_coms.txt | grep -p 'greedy\t118\t' | python AnalyzeResults/get_hpo_info.py 
"""

url = 'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo'
hpo_g = obonet.read_obo(url)

for line in sys.stdin:
    row = line.strip().split('\t')
    print(row)
    for item in row:
        try:
            print(item + ':' + hpo_g.nodes[item]['name'] + ' - ' + hpo_g.nodes[item]['def'])
        except:
            print(item + ' Not Found')
            continue


