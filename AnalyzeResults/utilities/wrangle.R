create_stats_df <- function(bocc, stats){
  map2_dfr(.x = bocc, .y = names(bocc), ~{
    df <- data.frame(method = .y, cluster_id = .x$cluster_id)
    vstat <- .x[, stats]
    dplyr::bind_cols(df, vstat)
  }) %>%
    dplyr::select(method, cluster_id, dplyr::everything()) %>%
    dplyr::as_tibble()
}

swap_if <- function(df, to_compare, to_swap) {
  x1 <- df[, to_compare[1]]
  x2 <- df[, to_compare[2]]
  should_swap <- ifelse(x1 <= x2, 0, 1)
  df[should_swap == 1, to_swap] <-
    df[should_swap == 1, rev(to_swap)]
  
  return(df)
}