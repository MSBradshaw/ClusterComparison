library(dplyr)
library(purrr)

# function to read BOCC results
read_BOCCResults <- function(path){
  path %>%
    list.files(., full.names = TRUE) %>%
    purrr::map(., ~ {
      result <- readr::read_delim(.x, delim = '\t')
      nms <- names(result)[-length(names(result))]
      result <- result[,-1]
      names(result) <- nms
      comma_sep_column <- purrr::map_lgl(result, ~any(stringr::str_detect(as.character(.x), ',')))
      result_nonlist <- result[which(!comma_sep_column)]
      result_list <- result[which(comma_sep_column)]
      convert_list <- list()
      for (i in 1:ncol(result_list)) {
        col <- result_list[[i]]
        convert_list[[i]] <- map(col, ~type.convert(stringr::str_split(.x, ',', simplify = TRUE)))
      }
      names(convert_list) <- names(result_list)
      tbl <- dplyr::bind_cols(result_nonlist, dplyr::as_tibble(convert_list))
      dplyr::select(tbl, tidyselect::any_of(nms))
    }) %>%
    purrr::set_names(., stringr::str_remove_all(list.files(path), '.results'))
}

# function to read discovered communities
read_Communities <- function(path){
  path %>%
    list.files(., full.names = TRUE) %>%
    purrr::map(., ~ {
      result <- readr::read_delim(.x, delim = ' ', col_names = FALSE)
      result$X1 %>%
        map(., ~ {
          v <- stringr::str_split(.x,
                                  pattern = '\t',
                                  n = 2,
                                  simplify = TRUE)
          stringr::str_split(v[-1], pattern = '\t', simplify = TRUE)
        }) %>%
        set_names(., ~ stringr::str_c('C', 0:(length(.) - 1)))
    })
}