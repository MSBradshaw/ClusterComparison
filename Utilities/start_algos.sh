# WARNING do not start these all at once, only 5 at a time, it will crash your computer

screen python Algorithms/greedy.py --edgelist Edgelists/String_HPO_2015.all_hpo.edgelist.txt --out greedy.String_HPO_2015.all_hpo.coms.txt\n
screen python Algorithms/greedy.py --edgelist Edgelists/String_HPO_2021.all_hpo.edgelist.txt --out greedy.String_HPO_2021.all_hpo.coms.txt\n
screen python Algorithms/greedy.py --edgelist Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt --out greedy.String_HPO_2015.phenotypic_branch.coms.txt\n
screen python Algorithms/greedy.py --edgelist Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt --out greedy.String_HPO_2021.phenotypic_branch.coms.txt\n
screen python Algorithms/greedy.py --edgelist Edgelists/String_HPO_2015.pruned.edgelist.txt --out greedy.String_HPO_2015.pruned.coms.txt\n

screen python Algorithms/walktrap.py --edgelist Edgelists/String_HPO_2015.all_hpo.edgelist.txt --out walktrap.String_HPO_2015.all_hpo.coms.txt\n
screen python Algorithms/walktrap.py --edgelist Edgelists/String_HPO_2021.all_hpo.edgelist.txt --out walktrap.String_HPO_2021.all_hpo.coms.txt\n
screen python Algorithms/walktrap.py --edgelist Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt --out walktrap.String_HPO_2015.phenotypic_branch.coms.txt\n
screen python Algorithms/walktrap.py --edgelist Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt --out walktrap.String_HPO_2021.phenotypic_branch.coms.txt\n
screen python Algorithms/walktrap.py --edgelist Edgelists/String_HPO_2015.pruned.edgelist.txt --out walktrap.String_HPO_2015.pruned.coms.txt\n

screen venv/bin/python Algorithms/cluster_infomap.py --edgelist Edgelists/String_HPO_2015.all_hpo.numbered.edgelist.txt --output infomap.String_HPO_2015.all_hpo.coms.txt --nodenames Edgelists/String_HPO_2015.all_hpo.nodenames.txt\n
screen venv/bin/python Algorithms/cluster_infomap.py --edgelist Edgelists/String_HPO_2021.all_hpo.numbered.edgelist.txt --output infomap.String_HPO_2021.all_hpo.coms.txt --nodenames Edgelists/String_HPO_2021.all_hpo.nodenames.txt\n
screen venv/bin/python Algorithms/cluster_infomap.py --edgelist Edgelists/String_HPO_2015.phenotypic_branch.numbered.edgelist.txt --output infomap.String_HPO_2015.phenotypic_branch.coms.txt --nodenames Edgelists/String_HPO_2015.phenotypic_branch.nodenames.txt\n
screen venv/bin/python Algorithms/cluster_infomap.py --edgelist Edgelists/String_HPO_2021.phenotypic_branch.numbered.edgelist.txt --output infomap.String_HPO_2021.phenotypic_branch.coms.txt --nodenames Edgelists/String_HPO_2021.phenotypic_branch.nodenames.txt\n
screen venv/bin/python Algorithms/cluster_infomap.py --edgelist Edgelists/String_HPO_2015.pruned.numbered.edgelist.txt --output infomap.String_HPO_2015.pruned.coms.txt --nodenames Edgelists/String_HPO_2015.pruned.nodenames.txt\n

