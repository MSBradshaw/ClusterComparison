create_stats_df <- function(bocc, stats){
  map2_dfr(.x = bocc, .y = names(bocc), ~{
    df <- data.frame(method = .y, cluster_id = .x$cluster_id)
    vstat <- .x[, stats]
    dplyr::bind_cols(df, vstat)
  }) %>%
    dplyr::select(method, cluster_id, dplyr::everything()) %>%
    dplyr::as_tibble()
}