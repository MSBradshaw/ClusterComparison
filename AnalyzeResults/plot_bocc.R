library(purrr)
library(dplyr)
library(ggplot2)
source('AnalyzeResults/utilities/read_write.R')
source('AnalyzeResults/utilities/wrangle.R')
source('AnalyzeResults/utilities/plot.R')

bocc <- read_BOCCResults(path = 'Results/BOCC_Results')

stats2plot <- c('cluster_size', 'gene_ratio', 'HPO_ratio', 'num_sig_go_enrichment_terms', 
                'max_norm_cell_type_specificity',
                'max_norm_disease_specificity')

df <- bocc %>%
  create_stats_df(., stats2plot) %>%
  dplyr::mutate(., dplyr::across(c('gene_ratio', 'HPO_ratio'), ~round(.x, 3)),
                sig_prop = num_sig_go_enrichment_terms/cluster_size)

# raw results
p1 <- ggplot(
  df,
  aes(
    x = max_norm_cell_type_specificity,
    y = max_norm_disease_specificity,
    size = sig_prop,
    color = method
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer('method', palette = 'Set1') +
  scale_size('sig prop') +
  labs(x = 'cell type specificity', y = 'disease type specificity') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.box = 'vertical')

# results without communities of size 1
p2 <- ggplot(
  df %>%
    dplyr::filter(., cluster_size > 1),
  aes(
    x = max_norm_cell_type_specificity,
    y = max_norm_disease_specificity,
    size = sig_prop,
    color = method,
    label = cluster_size,
    label2 = num_sig_go_enrichment_terms,
    label3 = gene_ratio,
    label4 = HPO_ratio
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer('method', palette = 'Set1') +
  scale_size('sig prop') +
  labs(x = 'cell type specificity', y = 'disease type specificity') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.box = 'vertical')
ggsave('AnalyzeResults/figs/bocc_results.png', width = 5, height = 5)

p3 <- p2 + 
  scale_x_log10(
    bquote(~log[10]('cell specificity')),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    bquote(~log[10]('disease specificity')),
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    minor_breaks = log10_minor_break(),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks()
ggsave('AnalyzeResults/figs/bocc_results_loglog.png', width = 5, height = 5)

plotly::ggplotly(p2)

df %>%
  group_by(method) %>%
  summarize(., n = length(cluster_id),
            min_m = min(cluster_size), mean_m = mean(cluster_size), sd_m = sd(cluster_size), max_m = max(cluster_size),
            sum_m = sum(cluster_size)) %>%
  kableExtra::kable(., format = 'latex', digits = 2, booktabs = TRUE) %>%
  kableExtra::kable_styling(position = 'center') %>%
  kableExtra::add_header_above(., list(' ' = 2, 'cluster size' = 5))
