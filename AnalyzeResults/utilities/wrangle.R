#' Subsets and formats the BOCC results for easy plotting.
#'
#' @param bocc The BOCC results, read using `read_BOCCResults()` in `utilities`.
#' @param stats A vector of column names. These should be columns to include in a plot.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' bocc <- read_BOCCResults(path = 'Results/BOCC_Results')
#' df <- create_stats_df(bocc, c('cluster_size', 'max_norm_cell_type_comma_sep_string'))
create_stats_df <- function(bocc, stats) {
  map2_dfr(.x = bocc, .y = names(bocc), ~ {
    df <- data.frame(method = .y, cluster_id = .x$cluster_id)
    vstat <- .x[, stats]
    dplyr::bind_cols(df, vstat)
  }) %>%
    dplyr::select(method, cluster_id, dplyr::everything()) %>%
    dplyr::as_tibble()
}

#' Swap values in two columns based on the order of values in another two columns of a data frame.
#'
#' @param df Dataframe containing columns to compare and swap.
#' @param to_compare A vector of size two of column names. If the value of the first column is larger than the value of the second column, then the values in the columns `to_swap` are swapped.
#' @param to_swap A vector of size two of column names. If the value of the first column of `to_compare` is larger than the value of the second column of `to_compare`, then the values in the columns `to_swap` are swapped.
#'
#' @return Dataframe with the proper ordering of columns.
#' @export
#'
#' @examples
#' df <- data.frame(a = c(1, 8, 3), b = c(1, 2, 27))
#' swap_if(df, to_compare = c('a', 'b'), to_swap = c('a', 'b'))
swap_if <- function(df, to_compare, to_swap) {
  x1 <- df[, to_compare[1]]
  x2 <- df[, to_compare[2]]
  should_swap <- ifelse(x1 <= x2, 0, 1)
  df[should_swap == 1, to_swap] <-
    df[should_swap == 1, rev(to_swap)]
  
  return(df)
}