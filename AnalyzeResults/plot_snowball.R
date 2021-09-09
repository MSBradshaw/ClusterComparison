library(dplyr)
library(ggplot2)
library(igraph)
source('AnalyzeResults/utilities/read_write.R')
source('AnalyzeResults/utilities/plot.R')

files <-
  list.files('Data/PhenotypicAbnormality',
             pattern = 'snowball',
             full.names = TRUE)

snowball <- files %>%
  purrr::map(., readr::read_tsv) %>%
  purrr::set_names(
    .,
    files %>%
      stringr::str_extract(., 'snowball.*.String') %>%
      stringr::str_remove_all(., 'snowball|String|\\.')
  )

# plot one for reference
snowball_results <- function(snowball_list, method, id) {
  new_edges <- snowball_list %>%
    .[[method]] %>%
    dplyr::filter(., com_id == !!id)
  
  pval <- mean(new_edges$replicate_score >= new_edges$com_score[1])
  
  df <- data.frame(
    method = method,
    com_id = new_edges$com_id[1],
    com_score = new_edges$com_score[1],
    com_size = new_edges$rep_and_com_size[1],
    com_density = new_edges$com_score[1] / choose(new_edges$rep_and_com_size[1], 2),
    pval = pval
  )
  df$com_density[is.nan(df$com_density)] <- 0
  
  return(df)
}

snowball_plot <- function(snowball_list, method, id, save = NULL) {
  new_edges <- snowball_list %>%
    .[[method]] %>%
    dplyr::filter(., com_id == !!id)
  
  h <- hist(new_edges$replicate_score, plot = FALSE)
  p.hist <- ggplot(data = new_edges, aes(x = replicate_score)) +
    geom_histogram(bins = length(h$counts)) +
    geom_vline(xintercept = new_edges$com_score[1],
               linetype = 'dashed') +
    labs(x = 'number of new edges (e)',
         y = 'frequency',
         title = stringr::str_c(stringr::str_to_lower(method), ' : ', prettyNum(id, big.mark = ','))) +
    theme_bw()
  
  p.ecdf <- ggplot(data = new_edges, aes(x = replicate_score)) +
    stat_ecdf(geom = 'step') +
    geom_vline(xintercept = new_edges$com_score[1],
               linetype = 'dashed') +
    labs(x = 'number of new edges (e)',
         y = bquote('P(E' < 'e)'),
         title = stringr::str_c(stringr::str_to_lower(method), ' : ', prettyNum(id, big.mark = ','))) +
    lims(x = c(0, new_edges$com_score[1] + new_edges$com_score[1]/4)) +
    theme_bw()
  
  if(!is.null(save)){
    ggsave(
      filename = stringr::str_c(save$path, '/', method, id, '-hist.png'),
      plot = p.hist,
      height = save$h,
      width = save$w
    )
    ggsave(
      filename = stringr::str_c(save$path, '/', method, id, '-ecdf.png'),
      plot = p.ecdf,
      height = save$h,
      width = save$w
    )
  }
  
  return(list(hist = p.hist, ecdf = p.ecdf))
}

results <- names(snowball) %>%
  purrr::map_dfr(., function(m) {
    ids <- sort(unique(snowball[[m]]$com_id))
    purrr::map_dfr(ids, function(i) snowball_results(snowball, method = m, id = i))
  })

# read in graph and communities
G <-
  read_graph(edgelist = 'Data/PhenotypicAbnormality/String_HPO_2015.phenotypic_branch.numbered.edgelist.txt',
             nodes = 'Data/PhenotypicAbnormality/String_HPO_2015.phenotypic_branch.nodenames.txt',
             communities = 'Data/PhenotypicAbnormality/Communities/2015')

comms <- V(G)$communities %>%
  dplyr::bind_rows() %>%
  dplyr::filter(., !is.na(cluster_id)) %>%
  dplyr::left_join(., y = data.frame(node = V(G)$node, isGene = V(G)$gene),
                   by = 'node') %>%
  dplyr::mutate(., method = stringr::str_extract(method, paste0(unique(results$method), collapse = '|')))

target_results <- results %>%
  dplyr::left_join(., y = comms %>%
                     dplyr::group_by(., method, cluster_id) %>%
                     dplyr::summarise(., mix = any(isGene == 1) & any(isGene == 0)) %>%
                     dplyr::ungroup() %>%
                     dplyr::rename(., com_id = cluster_id),
                   by = c('method', 'com_id')) %>%
  dplyr::filter(., com_size >= 3, mix)

# plots
jitter <- position_jitter(width = 0.2, height = 0)
# full violin
ggplot(data = target_results, aes(x = method, y = pval, fill = method, label = com_id)) +
  geom_violin() +
  geom_point(alpha = 0.6, position = jitter) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  scale_color_brewer('method', palette = 'Set1') +
  scale_fill_brewer('method', palette = 'Set1') +
  labs(x = 'method', y = 'p-value') +
  lims(y = c(0, 1)) +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('AnalyzeResults/figs/violin-full-snowball-pvals.png', width = 6, height = 4)

# significant violin
ggplot(data = dplyr::filter(target_results, pval < 0.1), 
       aes(x = method, y = pval, fill = method, label = com_id)) +
  geom_violin() +
  geom_point(alpha = 0.6, position = jitter) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  scale_color_brewer('method', palette = 'Set1') +
  scale_fill_brewer('method', palette = 'Set1') +
  labs(x = 'method', y = 'p-value') +
  lims(y = c(0, 0.1)) +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('AnalyzeResults/figs/violin-sig-snowball-pvals.png', width = 6, height = 4)

plotly::ggplotly()

# examples of the snowball test
snowball_plot(
  snowball,
  method = 'cesna',
  id = 2,
  save = list(path = 'AnalyzeResults/figs/snowball_examples',
              w = 2.5, h = 2.5)
)

snowball_plot(
  snowball,
  method = 'greedy',
  id = 0,
  save = list(path = 'AnalyzeResults/figs/snowball_examples',
              w = 2.5, h = 2.5)
)

snowball_plot(
  snowball,
  method = 'infomap',
  id = 6,
  save = list(path = 'AnalyzeResults/figs/snowball_examples',
              w = 2.5, h = 2.5)
)

snowball_plot(
  snowball,
  method = 'walktrap',
  id = 53,
  save = list(path = 'AnalyzeResults/figs/snowball_examples',
              w = 2.5, h = 2.5)
)

jitter <- position_jitter(height = 0.01)
ggplot(data = target_results, aes(x = com_size, y = pval)) +
  geom_point(alpha = 0.7, position = jitter) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  scale_x_log10(
    bquote(~log[10]('order of community')),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~method) +
  labs(y = 'p-value') +
  theme_bw()
ggsave('AnalyzeResults/figs/point-snowball-pvalsbyorderofcomm.png', height = 4, width = 4)

jitter <- position_jitter(height = 0.05)
ggplot(data = target_results %>%
         dplyr::mutate(., com_size = dplyr::case_when(
           com_size > 1 & com_size <= 10 ~ '1 - 10',
           com_size > 10 & com_size <= 100 ~ '11 - 100',
           com_size > 100 & com_size <= 1000 ~ '101 - 1000',
           com_size > 1000 & com_size <= 10000 ~ '1001 - 10000',
           com_size > 10000 ~ '10001 - '
         ) %>%
           factor(., levels = c('1 - 10', '11 - 100', '101 - 1000', '1001 - 10000', '10001 - '))),
       aes(x = com_density, y = pval, color = com_size)) +
  geom_point(alpha = 0.7, position = jitter) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  facet_wrap(~method) +
  scale_color_discrete('order') +
  labs(x = 'density of community', y = 'p-value') +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('AnalyzeResults/figs/point-snowball-pvalsbydensitycolororder.png', width = 6, height = 4.5)

ggplot(data = target_results %>%
         dplyr::mutate(., com_size = dplyr::case_when(
           com_size > 1 & com_size <= 10 ~ '1 - 10',
           com_size > 10 & com_size <= 100 ~ '11 - 100',
           com_size > 100 & com_size <= 1000 ~ '101 - 1000',
           com_size > 1000 & com_size <= 10000 ~ '1001 - 10000',
           com_size > 10000 ~ '10001 - '
         ) %>%
           factor(., levels = c('1 - 10', '11 - 100', '101 - 1000', '1001 - 10000', '10001 - '))),
       aes(x = com_density, y = pval, color = com_size)) +
  geom_point(alpha = 0.7, position = jitter) +
  geom_hline(yintercept = 0.1, linetype = 'dashed') +
  scale_x_log10(
    bquote(~log[10]('density of community')),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_color_discrete('order') +
  facet_wrap(~method) +
  labs(y = 'p-value') +
  theme_bw() +
  theme(legend.position = 'bottom')
ggsave('AnalyzeResults/figs/point-snowball-pvalsbydensitycolororder-log.png', width = 6, height = 4.5)

target_results %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(., n = n(), n_sig = sum(pval < 0.1),prop_sig = mean(pval < 0.1))

sig_results <- dplyr::filter(target_results, pval < 0.1)

purrr::map2(.x = sig_results$method, .y = sig_results$com_id, 
            ~ snowball_plot(snowball, method = .x, id = .y, 
                            save =  list(path = 'AnalyzeResults/figs/snowball_significant_examples',
                                         w = 2.5, h = 2.5)))

for(i in 1:nrow(sig_results)){
  cat(i, '\n')
  png(
    filename = stringr::str_c(
      'AnalyzeResults/figs/snowball_significant_examples/comms',
      '/', sig_results$method[i], sig_results$com_id[i],
      '.png'
    ),
    width = 5,
    height = 5,
    units = 'in',
    res = 100
  )
  nodes <- unique(dplyr::filter(comms, method == !!sig_results$method[i], cluster_id == !!sig_results$com_id[i])$node)
  plot_subgraph(G, 
                nodes = nodes, 
                maxn = 400,
                order = 1, 
                vertex.size = 5, 
                vertex.label.cex = 0.8, 
                layout = igraph::layout_with_fr)
  dev.off()
}

# look at characteristics of the target communities with nonsignificant results 
# and the target communities with the target communities with significant results
target_results <- dplyr::mutate(target_results, sig = pval < 0.1)
bocc <- read_BOCCResults('Data/PhenotypicAbnormality/BOCC') %>%
  purrr::set_names(., ~stringr::str_extract(.x, pattern = '^([^.]+)')) %>%
  dplyr::bind_rows(., .id = 'method') %>% 
  dplyr::select(., method, com_id = cluster_id, gene_ratio, HPO_ratio, num_sig_go_enrichment_terms, dplyr::ends_with('specificity'))
target_results <- target_results %>%
  dplyr::left_join(., bocc) %>%
  dplyr::filter(., mix)

# density
ggplot(target_results, aes(x = com_density, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = bquote(log[10]('density of community'))) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$com_density,
        dplyr::filter(target_results, sig == TRUE)$com_density) # < 0.001
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-commdensity.png', width = 5, height = 4)

# disease type specificity
ggplot(target_results, aes(x = max_norm_disease_specificity, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = bquote(log[10]('max norm disease specificity'))) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$max_norm_disease_specificity,
        dplyr::filter(target_results, sig == TRUE)$max_norm_disease_specificity) # 0.009
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-diseasespec.png', width = 5, height = 4)

# number of go terms
ggplot(target_results, aes(x = num_sig_go_enrichment_terms, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = bquote(log[10]('number of significant GO terms'))) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$num_sig_go_enrichment_terms,
        dplyr::filter(target_results, sig == TRUE)$num_sig_go_enrichment_terms) # 0.021
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-siggoterms.png', width = 5, height = 4)

# number of members
ggplot(target_results, aes(x = com_size, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = bquote(log[10]('number of members in community'))) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$com_size,
        dplyr::filter(target_results, sig == TRUE)$com_size) # = 0.031
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-commorder.png', width = 5, height = 4)

# ratio of genes
ggplot(target_results, aes(x = gene_ratio, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = 'gene ratio') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$gene_ratio,
        dplyr::filter(target_results, sig == TRUE)$gene_ratio) # 0.257
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-commgeneratio.png', width = 5, height = 4)

# cell type specificity
ggplot(target_results, aes(x = max_norm_cell_type_specificity, fill = sig)) +
  geom_density(alpha = 0.7) +
  scale_x_log10() +
  scale_fill_brewer('p < 0.1', palette = 'Set1') +
  labs(x = bquote(log[10]('cell type specificity'))) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8))
ks.test(dplyr::filter(target_results, sig == FALSE)$max_norm_cell_type_specificity,
        dplyr::filter(target_results, sig == TRUE)$max_norm_cell_type_specificity) # 0.354
ggsave('AnalyzeResults/figs/snowball_filter_nonsig/density-celltypespec.png', width = 5, height = 4)

