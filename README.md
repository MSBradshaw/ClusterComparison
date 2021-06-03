# ClusterComparison

![example workflow](https://github.com/MSBradshaw/ClusterComparison/actions/workflows/ci.yml/badge.svg)

Contained within this reposity are the methods and results of a study on the effectiveness of clustering on biological ontologies

## Installation & set up

Clone this repo

`git clone git@github.com:MSBradshaw/ClusterComparison.git`

Move into the repo

`cd ClusterComparison/`

Set up the conda enviroment (see [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for how to install conda)

`conda config --append channels conda-forge`

`conda create --name bocc --file requirements.txt`

Start the conda enviroment

`conda activate bocc`

Uncompress data

`gunzip BOCC/all_genes_info.json.gz`

## How to use BOCC to summarize a set of communities

`python BOCC/summarize_clusters.py --coms Data/three_communities.txt --pval 0.00005 --out results.tsv`

`--coms` should be a tab seporated file where each row represents a community. The first item in each row needs to be the cluster id/name, all subsequent items are members of the community

`--pval` should be a floating point number used as the threshold for significance in GO enrichment (default 0.00003)

`--out` name of the output file, should end in `.tsv`
