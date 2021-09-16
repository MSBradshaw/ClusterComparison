cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2021.all_hpo.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2021.all_hpo.nodenames.txt -c -1 -nt 1 -o String_HPO_2021.all_hpo
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2021.all_hpo.coms.txt --node_names Edgelists/String_HPO_2021.all_hpo.nodenames.txt 
ls cesna.String_HPO_2021.all_hpo.coms.txt
wc -l  cesna.String_HPO_2021.all_hpo.coms.txt

cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2021.phenotypic_branch.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2021.phenotypic_branch.nodenames.txt -c -1 -nt 1 -o String_HPO_2021.phenotypic_branch
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2021.phenotypic_branch.coms.txt --node_names Edgelists/String_HPO_2021.phenotypic_branch.nodenames.txt 
ls cesna.String_HPO_2021.phenotypic_branch.coms.txt
wc -l  cesna.String_HPO_2021.phenotypic_branch.coms.txt

cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2021.pruned.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2021.pruned.nodenames.txt -c -1 -nt 1 -o String_HPO_2021.pruned
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2021.pruned.coms.txt --node_names Edgelists/String_HPO_2021.pruned.nodenames.txt 
ls cesna.String_HPO_2021.pruned.coms.txt
wc -l  cesna.String_HPO_2021.pruned.coms.txt

cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2015.all_hpo.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2015.all_hpo.nodenames.txt -c -1 -nt 1 -o String_HPO_2015.all_hpo
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2015.all_hpo.coms.txt --node_names Edgelists/String_HPO_2015.all_hpo.nodenames.txt 
ls cesna.String_HPO_2015.all_hpo.coms.txt
wc -l  cesna.String_HPO_2015.all_hpo.coms.txt

cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2015.phenotypic_branch.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2015.phenotypic_branch.nodenames.txt -c -1 -nt 1 -o String_HPO_2015.phenotypic_branch
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2015.phenotypic_branch.coms.txt --node_names Edgelists/String_HPO_2015.phenotypic_branch.nodenames.txt 
ls cesna.String_HPO_2015.phenotypic_branch.coms.txt
wc -l  cesna.String_HPO_2015.phenotypic_branch.coms.txt

cd snap/examples/cesna/
./cesna -i ../../../Edgelists/String_HPO_2015.pruned.numbered.edgelist.txt -l ../../../Edgelists/String_HPO_2015.pruned.nodenames.txt -c -1 -nt 1 -o String_HPO_2015.pruned
cd ../../..
python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna.String_HPO_2015.pruned.coms.txt --node_names Edgelists/String_HPO_2015.pruned.nodenames.txt 
ls cesna.String_HPO_2015.pruned.coms.txt
wc -l  cesna.String_HPO_2015.pruned.coms.txt
