# ClusterComparison

![example workflow](https://github.com/MSBradshaw/ClusterComparison/actions/workflows/ci.yml/badge.svg)

Contained within this reposity are the methods and results of a study on the effectiveness of clustering on biological ontologies

## Installation & set up

Clone this repo

`git clone git@github.com:MSBradshaw/ClusterComparison.git`

Move into the repo

`cd ClusterComparison/`

Create an environment

`python3 -m venv env`

Start the enviroment
`source env/bin/activate`

Install requirements
`python3 -m pip install -r requirements.txt`

Uncompress data

`gunzip BOCC/all_genes_info.json.gz`

## CESNA

### Install SNAP

`git clone git@github.com:snap-stanford/snap.git`

`cd snap/examples/cesna`

`make`

### Run CENSA

Change `-nt` to the number of threads to be used for parallelization

`./cesna -i ../../../Data/HPO_String_edgelist.numbered.tsv -l ../../../Data/HPO_String_edgelist.nodenames.tsv -c -1 -nt 1 -o hpo_string_cesna`

Create community file with CESNA results

`python Algorithms/number_cesna_results.py --cesna_res snap/examples/cesna/cmtyvv.txt --output cesna_coms.txt --node_names Data/HPO_String_edgelist.nodenames.tsv`

## Greedy, Walk Trap and Belief

`belief.py`, `greedy.py`, and `walktrap.py` are all run in the same manner

`python Algorithms/greedy.py --edgelist Data/HPO_String_edgelist.tsv --out greedy_coms.txt`

## Infomap

`python Algorithms/cluster_infomap.py --edgelist Data/HPO_String_edgelist.numbered.tsv --output infomap_coms.txt --nodenames Data/HPO_String_edgelist.nodenames`

## How to use BOCC to summarize a set of communities

`python BOCC/summarize_clusters.py --coms Data/three_communities.txt --alpha 0.05 --mg2 Data/mygene2_gene_hpo_family.tsv --out results.tsv`

`--coms` should be a tab seporated file where each row represents a community. The first item in each row needs to be the cluster id/name, all subsequent items are members of the community

`--mg2` path to a tab seporated files with the MyGene2 data with the follow columns: gene, HPO, family ID

`--alpha` floating point number 0-1, threshold for signifiance in FDR correction for GO enrichment (default 0.05)

`--out` name of the output file, should end in `.tsv`

## Plot Communities of interest

`python AnalyzeResults/plot_com.py --edgelist Data/HPO_String_edgelist_june_22_2021.tsv --all_coms all_june_22_coms.txt --algo walktrap --com 84 --output walktrap_84.png`

`--edgelist` path to edgelist of the whole network

`--all_coms` path to file with all communities each row formatted like: algorithm \t community_id \t member1 \t member2 ... 

`--algo` name of the algorithm of interst

`--com` community id of interst

`--output` where to save the figure

## Compare to synthetic nulls

### Snowball

`python snowball.py --edgelist edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt
--output snowball.infomap.String_HPO_2015.phenotypic_branch.tsv
--coms Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt
--new_edges Data/new_jenkins_edges.tsv
--reps 100`

`--edgelist` tab separated edge list

`--output` name of file to save the community to

`--coms` file of communities and there members, each row is a com, first item is the com name, all others are members

`--new_edges` file with new edges, tab separated

`--reps` number of repitions, default = 100

## Hierarchial Clustering

Hierarchial clister is done on each community produced by the 4 clustering algorithms. In this case each community is treated as it's own graph and uses a balanced cut with a max size of 200.

`python AnalyzeResults/hierarchical_clustering.py --algo paris --edgelist Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt --coms Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt --output SubComs/infomap.String_HPO_2015.phenotypic_branch.coms`

`--edgelist` tab separated edge list

`--output` prefix used for naming output files (there will be one output file for each community in the `--coms` files

`--coms` file of communities and there members, each row is a com, first item is the com name, all others are members

`--algo` algorithm to be used [`paris`, `ward`, `louvain`]


## Source Files:

`hp_2015_week_46.obo` taken from https://github.com/drseb/HPO-archive/blob/master/2014-2015/2015_week_46/hpo/artefacts/hp.obo.gz

`ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_12_2015.txt`  downloaded from https://github.com/drseb/HPO-archive/blob/master/hpo.annotations.monthly/2015-12-01_00-00-05/archive/annotation/ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes.txt.gz

Old school string

http://string91.embl.de/newstring_cgi/show_download_page.pl?UserId=fcQWm48WiULD&sessionId=PJ_uFqbyyWmH

protein.links.v9.1.txt downloaded from http://string91.embl.de/newstring_download/protein.links.v9.1.txt.gz

Homo Spanien only

`9606.protein.links.v9.1.txt` http://string91.embl.de/newstring_download/protein.links.v9.1/9606.protein.links.v9.1.txt.gz

`9606.protein.links.detailed.v9.1.txt` from http://string91.embl.de/newstring_download/protein.links.detailed.v9.1/9606.protein.links.detailed.v9.1.txt.gz

Modern String

`9606.protein.info.v11.0.txt` from https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz on July 8 2021

`9606.protein.links.v11.0.txt` from https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz on July 8 2021

Modern HPO

`hp_July_8_2021.obo` from http://purl.obolibrary.org/obo/hp.obo on July 8 2021

`genes_to_phenotype.txt` from  http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt on July 8 2021

