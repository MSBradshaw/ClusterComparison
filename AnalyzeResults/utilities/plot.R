#' Add log10 minor breaks to a ggplot with ggplot2::scale_x/y_log10() activated.
#'
#' @param ... Added via a `+` in a `ggplot2` chain.
#'
#' @return A log10 plot with minor breaks.
#' @export
#'
#' @examples
#' ggplot2::ggplot() + ggplot2::scale_x_log10('x', minor_breaks = log10_minor_breaks())
log10_minor_break <- function (...) {
  function(x) {
    minx <- floor(min(log10(x), na.rm = T)) - 1
    maxx <- ceiling(max(log10(x), na.rm = T)) + 1
    n_major <- maxx - minx + 1
    major_breaks <- seq(minx, maxx, by = 1)
    minor_breaks <-
      rep(log10(seq(1, 9, by = 1)), times = n_major) +
      rep(major_breaks, each = 9)
    
    return(10 ^ (minor_breaks))
  }
}

#' A better implementation of `plot.igraph`.
#'
#' @param G An igraph object.
#' @param legend `TRUE` or `FALSE` depending on whether `G` is heterogeneous and one would like a legend indicating node type.
#' @param ... Any other `plot.igraph` options. See the [help guide](https://igraph.org/r/doc/plot.common.html) for more details.
#'
#' @return A plot.
#' @export
plot_graph <- function(G, legend = FALSE, ...) {
  V(G)$color <- purrr::map_chr(V(G)$gene + 1, ~RColorBrewer::brewer.pal(9, 'Set1')[.x])
  par(bg = NA, mar = c(0, 0, 0, 0) + 0.1)
  plot(G, ...)
  if(legend){
    types <- sort(unique(V(G)$gene))
    legend("bottomleft",
           bty = "n",
           legend = types,
           fill = RColorBrewer::brewer.pal(9, 'Set1'), 
           ncol = min(length(types), 4),
           border = NA)
  }
}

plot_subgraph <- function(G, nodes, maxn = 100, order = 1, ...) {
  if(length(nodes) <= maxn){
    lv <- V(G)$node %in% nodes
    sG <- igraph::induced_subgraph(G, vids = lv)
    neighbor_nodes <- igraph::neighborhood(G, nodes = lv, order = order) %>%
      purrr::reduce(., c) %>%
      unique() %>%
      .$node
    if(length(neighbor_nodes) >= maxn){
      neighbor_nodes <- nodes
      cat('Too big for neighbors', '\n')
    }
    sG_neigh <- igraph::induced_subgraph(G, vids = V(G)$node %in% neighbor_nodes)
    # set all colors to dark red and blue
    V(sG_neigh)$color <- RColorBrewer::brewer.pal(9, 'Set1')[V(sG_neigh)$gene + 1]
    # set to pale if a neighboor
    V(sG_neigh)[!V(sG_neigh)$node %in% V(sG)$node & V(sG_neigh)$gene == 1]$color <- '#c7e6ff' # grey-blue
    V(sG_neigh)[!V(sG_neigh)$node %in% V(sG)$node & V(sG_neigh)$gene == 0]$color <- '#fca9aa' # grey-red
    
    # set vertex label
    V(sG_neigh)$vertex_label <- NA
    if (igraph::vcount(sG) < 8){
      V(sG_neigh)[V(sG_neigh)$node %in% V(sG)$node]$vertex_label <-
        V(sG_neigh)[V(sG_neigh)$node %in% V(sG)$node]$name
    }
    
    # change color and width of within community ties
    E(sG_neigh)$width <- 1
    E(sG_neigh)$color <- 'lightgray'
    oG_edges <- igraph::as_data_frame(sG)
    neigh_edges <- sG_neigh %>%
      igraph::as_data_frame() %>%
      tibble::rowid_to_column()
    eids <- setdiff(neigh_edges$rowid, dplyr::anti_join(neigh_edges, oG_edges, by = c('from', 'to'))$rowid)
    E(sG_neigh)[eids]$width <- 2
    E(sG_neigh)[eids]$color <- 'black'
    
    # create plot
    par(bg = NA, mar = c(0, 0, 0, 0) + 0.1)
    plot(
      sG_neigh,
      vertex.size = 6,
      vertex.label = V(sG_neigh)$vertex_label,
      vertex.color = V(sG_neigh)$color,
      vertex.label.color = 'black'
    )
  }else{
    NULL
  }
}
