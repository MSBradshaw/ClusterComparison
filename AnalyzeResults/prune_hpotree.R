library(igraph)
library(ggplot2)
library(purrr)
source('AnalyzeResults/utilities/read_write.R')
source('AnalyzeResults/utilities/plot.R')

G <-
  read_graph(edgelist = 'Data/HPO_String_edgelist_june_22_2021.numbered.tsv',
             nodes = 'Data/HPO_String_edgelist_june_22_2021.nodenames.tsv',
             communities = 'Results/Communities')

H <- igraph::induced_subgraph(G, vids = igraph::V(G)[igraph::V(G)$gene == 0])

# plot HPO tree
# png(file = "AnalyzeResults/figs/hpo-tree.png")
# plot_graph(H, layout = igraph::layout_as_tree, vertex.label = NA, vertex.size = 2, edge.size = 0.5)
# dev.off()

# get distance to root and see how many unique genes are connected at that level
distance_from_root <-
  igraph::distances(H, v = V(H)[V(H)$name == 'HP:0000001'])

neighbors_of_HPOs <-
  neighborhood(G,
               order = 1,
               mindist = 1,
               nodes = V(G)[V(G)$name %in% V(H)$name])

dist_sorted <- sort(unique(distance_from_root[1,]))

nunique_genes <- purrr::map_dbl(dist_sorted, ~ {
  hpos_at_dist <- which(distance_from_root[1,] == .x)
  sum(unique(purrr::reduce(neighbors_of_HPOs[hpos_at_dist], c))$gene == 1)
})

# plot number of unique genes at each distance to root
df <- tibble(dist_from_root = dist_sorted,
             uni_genes = nunique_genes)
ggplot(df, aes(x = factor(dist_from_root), y = uni_genes)) +
  geom_bar(stat = 'identity') +
  labs(x = 'distance from HP:0000001', y = 'number of unique genes', title = 'Number of Unique Genes by Distance from Root') +
  theme_bw()
# ggsave('AnalyzeResults/figs/unigenes-by-distance-from-root-in-hpo-tree.png')

# degree of one hops
cuts <- names(distance_from_root[1,][which(distance_from_root[1,] %in% c(0, 1))])

# remove cuts
Gstar <- igraph::induced_subgraph(G, vids = V(G)[V(G)$name %in% setdiff(V(G)$name, cuts)])
bG <- igraph::betweenness(G)
bGstar <- igraph::betweenness(Gstar)
