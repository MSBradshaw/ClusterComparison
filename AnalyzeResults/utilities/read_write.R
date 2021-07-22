require(dplyr)
require(purrr)
source('AnalyzeResults/utilities/wrangle.R')

# function to read BOCC results
read_BOCCResults <- function(path){
  comma_sep_column <- path %>%
    list.files(., full.names = TRUE) %>%
    purrr::map(., ~ {
      result <- readr::read_delim(.x, delim = '\t', col_names = TRUE)
      nms <- names(result)
      purrr::map_lgl(result, ~any(stringr::str_detect(as.character(.x), ',')))
    }) %>%
    purrr::transpose() %>%
    purrr::map(., ~unlist(.x)) %>%
    purrr::map_lgl(., ~any(.x))
    
  path %>%
    list.files(., full.names = TRUE) %>%
    purrr::map(., ~ {
      result <- readr::read_delim(.x, delim = '\t', col_names = TRUE)
      nms <- names(result)
      result_nonlist <- result[which(!comma_sep_column)]
      result_list <- result[which(comma_sep_column)]
      convert_list <- list()
      for (i in 1:ncol(result_list)) {
        col <- result_list[[i]]
        convert_list[[i]] <- map(col, ~type.convert(stringr::str_split(.x, ',', simplify = TRUE), as.is = TRUE))
      }
      names(convert_list) <- names(result_list)
      tbl <- dplyr::bind_cols(result_nonlist, dplyr::as_tibble(convert_list))
      dplyr::select(tbl, tidyselect::any_of(nms))
    }) %>%
    purrr::set_names(., stringr::str_remove_all(list.files(path), '.results'))
}

# function to read discovered communities
read_Communities <- function(path){
  comms_list <- path %>%
    list.files(., full.names = TRUE) %>%
    purrr::map(., ~ {
      method <- as.vector(stringr::str_match(.x, '.*/(.*?)_coms*'))
      method <- method[length(method)]
      result <- readr::read_delim(.x, delim = ' ', col_names = FALSE)
      comms <-
        purrr::map2_dfr(.x = result$X1, .y = rep(method, length(result$X1)), ~ {
          v <- stringr::str_split(.x,
                                  pattern = '\t',
                                  n = 2,
                                  simplify = TRUE)
          id <- v[1]
          v <-
            as.vector(stringr::str_split(v[-1], pattern = '\t', simplify = TRUE))
          data.frame(method = .y, cluster_id = id, name = v)
        })
    })
  methods <- path %>% 
    list.files(., full.names = TRUE) %>% 
    stringr::str_match(., '.*/(.*?)_coms*') %>% 
    .[,2]
  comms_list %>%
    purrr::set_names(., methods)
}

read_graph <- function(edgelist, nodes, communities){
  el <- suppressWarnings(
    edgelist %>%
      readr::read_delim(., delim = '\t', col_names = FALSE) %>%
      purrr::set_names(., c('from', 'to')) %>%
      # filter out multiedges
      swap_if(., to_compare = c('to', 'from'), to_swap = c('to', 'from')) %>%
      dplyr::distinct(., .keep_all = TRUE)
  )
  nl <- nodes %>%
    readr::read_delim(., delim = '\t', col_names = FALSE) %>%
    purrr::set_names(., c('nodeID', 'name')) %>%
    dplyr::mutate(., node = nodeID,
                  gene = as.numeric(stringr::str_sub(name, end = 3) != "HP:"))
  comms <- read_Communities(communities)
  
  for(method in names(comms)){
    comms2bind <- comms[[method]] %>%
      dplyr::select(., name, cluster_id) %>%
      purrr::set_names(., c('name', method))
    nl <- dplyr::left_join(nl, comms2bind, by = 'name')
  }
  
  nl <- nl %>%
    readr::type_convert() %>%
    tidyr::pivot_longer(., cols = tidyselect::any_of(names(comms)), names_to = 'method', values_to = 'cluster_id') %>%
    tidyr::nest(., communities = c(node, method, cluster_id)) %>%
    dplyr::mutate(., node = nodeID) %>%
    dplyr::select(., nodeID, node, name, dplyr::everything())
  
  G <- igraph::graph_from_data_frame(d = el, vertices = nl, directed = FALSE)
  
  return(G)
}

file <- 'Data/raw_data/genes_to_phenotype.txt'
read_genehpo_links <- function(file){
  `%>%` <- dplyr::`%>%`
  
  con <- file(file, 'r')
  first_line <- readLines(con, n = 1)
  close(con)
  
  first_line <- first_line %>%
    stringr::str_remove(., '\\#Format:') %>%
    stringr::str_squish() %>%
    stringr::str_split(., pattern = '<tab>', simplify = TRUE) %>%
    as.vector() %>%
    janitor::make_clean_names(.)
  
  tsv_file <- file %>%
    readr::read_tsv(., col_names = FALSE, skip = 1) %>%
    setNames(., first_line)
  
  gene_first <- stringr::str_detect(names(tsv_file[,1]), 'gene')
  if(gene_first){
    dplyr::select(tsv_file, dplyr::starts_with('entrez'), dplyr::starts_with('hpo'), dplyr::starts_with('disease'))
  }else{
    dplyr::select(tsv_file, dplyr::starts_with('hpo'), dplyr::starts_with('entrez'), dplyr::starts_with('disease'))
  }
}
