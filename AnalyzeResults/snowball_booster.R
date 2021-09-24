library(dplyr)
library(ggplot2)
library(igraph)
source('AnalyzeResults/utilities/read_write.R')
source('AnalyzeResults/utilities/plot.R')

# read snowball data
files <-
  list.files('Data/PhenotypicAbnormality/Snowball',
             pattern = 'snowball',
             full.names = TRUE)

snowball <- files %>%
  purrr::map(., readr::read_tsv) %>%
  purrr::set_names(
    .,
    files %>%
      stringr::str_extract(., 'snowball.*.String') %>%
      stringr::str_remove_all(., 'snowball|String|\\.')
  ) %>%
  dplyr::bind_rows(., .id = 'method') %>%
  dplyr::group_split(., method, com_id) %>%
  purrr::map_dfr(
    .,
    ~ data.frame(
      method = .x$method[1],
      com_id = .x$com_id[1],
      sb_pval = mean(.x$replicate_score >= .x$com_score[1])
    )
  )
rm(files)

# read bocc data
bocc_wds <- read_BOCCResults('Data/PhenotypicAbnormality/BOCCWithoutDiseaseSpecificity') %>%
  purrr::set_names(., ~stringr::str_extract(.x, pattern = '^([^.]+)')) %>%
  dplyr::bind_rows(., .id = 'method') %>%
  dplyr::select(., method, com_id = cluster_id, dplyr::everything(), 
                -dplyr::matches('go_'), 
                -dplyr::ends_with('_string'), 
                -dplyr::starts_with('mg2_'))

bocc_ds <- read_BOCCResults('Data/PhenotypicAbnormality/BOCC') %>%
  purrr::set_names(., ~stringr::str_extract(.x, pattern = '^([^.]+)')) %>%
  dplyr::bind_rows(., .id = 'method') %>%
  dplyr::select(., method, com_id = cluster_id, num_sig_go_enrichment_terms)
  
bocc <- dplyr::left_join(bocc_wds, bocc_ds, by = c('method', 'com_id'))
rm(bocc_wds, bocc_ds)

# read the data from the communities and get edge counts
G <-
  read_graph(edgelist = 'Data/PhenotypicAbnormality/String_HPO_2015.phenotypic_branch.numbered.edgelist.txt',
             nodes = 'Data/PhenotypicAbnormality/String_HPO_2015.phenotypic_branch.nodenames.txt',
             communities = 'Data/PhenotypicAbnormality/Communities/2015')

communities <- V(G)$communities %>%
  dplyr::bind_rows() %>%
  dplyr::filter(., !is.na(cluster_id)) %>%
  dplyr::left_join(., y = data.frame(node = V(G)$node, isGene = V(G)$gene),
                   by = 'node') %>%
  dplyr::mutate(., method = stringr::str_extract(method, 'greedy|cesna|infomap|walktrap')) %>%
  dplyr::select(., method, com_id = cluster_id, node, isGene)

edge_counts <- communities %>%
  dplyr::group_split(., method, com_id) %>%
  purrr::map_dfr(., ~ {
    sG <- igraph::induced_subgraph(G, vids = V(G)$node %in% .x$node)
    data.frame(method = .x$method[1], com_id = .x$com_id[1], nedges = igraph::ecount(sG))
  })
rm(G, communities)

# create model frame
df <- bocc %>%
  dplyr::left_join(., snowball, by = c('method', 'com_id')) %>%
  dplyr::left_join(., edge_counts, by = c('method', 'com_id')) %>%
  # keep only communties with a mixture of genes and hpos and at least 3 members
  dplyr::filter(., gene_ratio > 0, HPO_ratio > 0, cluster_size > 2, !is.na(sb_pval)) %>%
  dplyr::mutate(., edge_density = nedges/choose(cluster_size, 2),
                method = factor(method, levels = c('cesna', 'greedy', 'infomap', 'walktrap')))

# impute missingness
mdf <- mice::mice(dplyr::select(df, -method, -HPO_ratio, -conductance, -cut_ratio, -newman_girvan_modularity))
cdf <- mdf %>%
  mice::complete() %>%
  dplyr::bind_cols(., dplyr::select(df, method, HPO_ratio, conductance, cut_ratio, newman_girvan_modularity)) %>%
  dplyr::select(., dplyr::any_of(names(df)))
rm(df, mdf)

# fit booster
set.seed(652)

# preprocess data
cdf_preprocessed <- cdf %>%
  dplyr::mutate(., sb_uninteresting_comm = factor(ifelse(sb_pval > 0.999, 'yes', 'no'), 
                                                  levels = c('yes', 'no'))) %>%
  fastDummies::dummy_cols(., select_columns = 'method', remove_first_dummy = TRUE, remove_selected_columns = FALSE) %>%
  dplyr::select(., method, com_id, sb_pval, sb_uninteresting_comm, dplyr::starts_with('method_'), dplyr::everything(), -HPO_ratio)

# train and test split
cdf_split <-
  rsample::initial_split(cdf_preprocessed, prop = 0.6, strata = 'method')
val <- list(
  data = cdf_preprocessed,
  train = rsample::training(cdf_split),
  test = rsample::testing(cdf_split)
)

# run models
fit_control <- caret::trainControl(
  method = 'cv',
  number = 5,
  savePredictions = 'all',
  classProbs = TRUE,
  allowParallel = TRUE,
  verboseIter = TRUE
)

tuners <- list(
  xgbLinear = caretEnsemble::caretModelSpec(method = "xgbLinear", tuneLength = 5),
  xgbTree = caretEnsemble::caretModelSpec(method = "xgbTree", tuneLength = 5),
  xgbDART = caretEnsemble::caretModelSpec(method = "xgbDART", tuneLength = 5)
)

methods <- names(tuners)

model_list <- caretEnsemble::caretList(
  sb_uninteresting_comm ~ .-method-com_id-sb_pval,
  data = val$train,
  trControl = fit_control,
  methodList = methods,
  tuneList = tuners,
  continue_on_fail = FALSE
)

result <- list(data = val,
               method = methods,
               fit_control = fit_control,
               tuners = tuners,
               model_results = model_list)

readr::write_rds('Data/bigp_booster_results.rds')

# assess results
m <- readr::read_rds('Data/bigp_booster_results.rds')

preds <- m$model_results %>%
  purrr::map_dfc(., ~predict(.x, newdata = m$data$test)) %>%
  dplyr::bind_cols(., sb_uninteresting_comm = m$data$test$sb_uninteresting_comm) %>%
  tidyr::pivot_longer(., cols = -sb_uninteresting_comm, names_to = 'model', values_to = 'prediction') %>%
  dplyr::group_split(., model)

# CV best accuracy - xgbDART = 0.81
purrr::map_dbl(m$model_results, ~max(.x$results$Accuracy))

# confusion matrix
model_perf <- preds %>%
  purrr::map(., ~ {
    est <-
      yardstick::conf_mat(.x,
                          truth = sb_uninteresting_comm,
                          estimate = prediction)
    list(
      conf_mat = est,
      acc = sum(diag(est$table)) / sum(est$table),
      acc_given_pred = diag(est$table) / c(sum(est$table[1, ]), sum(est$table[2, ]))
    )
  }) %>%
  purrr::set_names(., names(m$model_results))
model_perf$xgbDART

# model performance of DART
vimp_dt <- m$model_results$xgbDART %>%
  caret::varImp() %>%
  .$importance %>%
  data.frame() %>%
  tibble::rownames_to_column(., 'var') %>%
  dplyr::rename(., var_imp = Overall)

ggplot(vimp_dt, aes(x = factor(var, var), y = var_imp)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = '', y = 'variable importance') +
  theme_bw()

# full results  
m$model_results$xgbDART %>%
  predict(., newdata = m$data$data) %>%
  dplyr::bind_cols(
    sb_uninteresting_comm = m$data$data$sb_uninteresting_comm,
    prediction = .
  ) %>%
  list(xgbDART = .) %>%
  purrr::map(., ~ {
    est <-
      yardstick::conf_mat(.x,
                          truth = sb_uninteresting_comm,
                          estimate = prediction)
    list(
      conf_mat = est,
      acc = sum(diag(est$table)) / sum(est$table),
      acc_given_pred = diag(est$table) / c(sum(est$table[1,]), sum(est$table[2,]))
    )
  })

# full results by type
by_method <- m$model_results$xgbDART %>%
  predict(., newdata = m$data$data) %>%
  dplyr::bind_cols(
    method = m$data$data$method,
    sb_uninteresting_comm = m$data$data$sb_uninteresting_comm,
    prediction = .
  ) %>%
  dplyr::group_split(., method) %>%
  purrr::map(., ~ {
    est <-
      yardstick::conf_mat(.x,
                          truth = sb_uninteresting_comm,
                          estimate = prediction)
    list(
      N = nrow(.x),
      ppl1 = sum(.x$prediction == 'no'),
      conf_mat = est,
      acc = sum(diag(est$table)) / sum(est$table),
      acc_given_pred = diag(est$table) / c(sum(est$table[1,]), sum(est$table[2,]))
    )
  }) %>%
  purrr::set_names(., c('cesna', 'greedy', 'infomap', 'walktrap'))

# contribution of variables
library(SHAPforxgboost)
train <- m$data$train %>%
  dplyr::select(., -dplyr::any_of(c('sb_uninteresting_comm', 'method', 'com_id', 'sb_pval'))) %>%
  as.matrix()

# shap values
shap_values <- SHAPforxgboost::shap.values(m$model_results$xgbDART$finalModel, X_train = train)

# shap plot
plot_shap <- function (data_long, 
                       x_bound = NULL) {
  require(SHAPforxgboost)
  
  label_format = "%.3f"
  N_features <- length(unique(data_long$variable))
  nrow_X <- nrow(data_long)/N_features
  x_bound <- if (is.null(x_bound)) {
    max(abs(data_long$value)) * 1.1
  } else {
    as.numeric(abs(x_bound))
  } 
  plot1 <- ggplot(data = data_long) + 
    coord_flip(ylim = c(-x_bound, x_bound)) + 
    geom_hline(yintercept = 0) + 
    ggforce::geom_sina(
      aes(x = variable, y = value, color = stdfvalue), 
      method = "counts",
      maxwidth = 0.7, 
      alpha = 0.5) + 
    geom_text(
      data = unique(data_long[,
                              c("variable", "mean_value")]),
      aes(
        x = variable,
        y = -Inf,
        label = sprintf(label_format, mean_value)
      ),
      size = 3,
      alpha = 0.7,
      hjust = -0.2,
      fontface = "bold"
    ) + 
    scale_color_gradient(
      low = "#FFCC33",
      high = "#6600CC",
      breaks = c(0, 1),
      labels = c(" Low", "High "),
      guide = guide_colorbar(barwidth = 12, barheight = 0.3)
    ) +
    theme_bw() + theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      axis.title.x = element_text(size = 10)
    ) +
    scale_x_discrete(limits = rev(levels(data_long$variable)),
                     labels = SHAPforxgboost:::label.feature(rev(levels(data_long$variable)))) +
    labs(y = "SHAP value (impact on model output)",
         x = "", color = "feature value  ")
  
  return(plot1)
}
plot_shap(shap_long, x_bound = 1)

shap_long <- SHAPforxgboost::shap.prep(m$model_results$xgbDART$finalModel, X_train = train)
SHAPforxgboost::shap.plot.summary(shap_long, x_bound = 1)

# side by side plots
vars <- stringr::str_subset(names(shap_values$mean_shap_score), 'method_', negate = TRUE)
v1t8 <- list()
v9t16 <- list()
# greedy
v1t8$greedy <- purrr::map(
  vars[1:8],
  ~ SHAPforxgboost::shap.plot.dependence(shap_long, x = .x, color_feature = 'method_greedy', smooth = FALSE) +
    geom_hline(yintercept = 0, linetype = 'dashed')
)
# infomap
v1t8$infomap <- purrr::map(
  vars[1:8],
  ~ SHAPforxgboost::shap.plot.dependence(shap_long, x = .x, color_feature = 'method_infomap', smooth = FALSE) +
    geom_hline(yintercept = 0, linetype = 'dashed')
)
# walktrap
v1t8$walktrap <- purrr::map(
  vars[1:8],
  ~ SHAPforxgboost::shap.plot.dependence(shap_long, x = .x, color_feature = 'method_walktrap', smooth = FALSE) +
    geom_hline(yintercept = 0, linetype = 'dashed')
)

ggpubr::ggarrange(plotlist = v1t8$greedy, nrow = 2, ncol = 4, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(plotlist = v1t8$infomap, nrow = 2, ncol = 4, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(plotlist = v1t8$walktrap, nrow = 2, ncol = 4, common.legend = TRUE, legend = 'bottom')

# variable importance shap
df <- tibble::rownames_to_column(data.frame(mean_shap = shap_values$mean_shap_score), var = 'var')
ggplot(data = df, aes(x = factor(var, levels = var), y = mean_shap)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  labs(x = '', y = 'mean absolute SHAP value') +
  theme_bw()

