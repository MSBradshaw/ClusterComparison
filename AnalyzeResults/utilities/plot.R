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
