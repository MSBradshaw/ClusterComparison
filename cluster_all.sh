rj -n infomap -m 30G -T 4 -c 'python Algorithms/cluster_info.py --edgelist Data/HPO_String_edgelist.tsv --out infomap_coms.txt'
rj -n belief -t 72:00:00 -m 30G -T 4 -c 'python Algorithms/belief.py --edgelist Data/HPO_String_edgelist.tsv --out belief_coms.txt'
rj -n greedy -m 30G -T 4 -c 'python Algorithms/greedy.py --edgelist Data/HPO_String_edgelist.tsv --out greedy_coms.txt'
rj -n walktrap -m 30G -T 4 -c 'python Algorithms/walktrap.py  --edgelist Data/HPO_String_edgelist.tsv --out walktrap_coms.txt'
