# ClusterComparison

![example workflow](https://github.com/MSBradshaw/ClusterComparison/actions/workflows/ci.yml/badge.svg)

Contained within this reposity are the methods and results of a study on the effectiveness of clustering on biological ontologies

## How to use BOCC to summarize a set of communities

`python BOCC/summarize_clusters.py --coms Data/three_communities.txt --pval 0.00005 --out results.tsv`

`--coms` should be a tab seporated file where each row represents a community. The first item in each row needs to be the cluster id/name, all subsequent items are members of the community

`--pval` should be a floating point number used as the threshold for significance in GO enrichment (default 0.00003)

`--out` name of the output file, should end in `.tsv`
