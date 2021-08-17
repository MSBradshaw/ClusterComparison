#' Calculate the average jaccard similarity coefficient across two sets of communities. This measure gauges the similarity between two community strucutres with possible overlap.
#'
#' @param A A list of vectors.
#' @param B A list of vectors.
#'
#' @return A double between 0 (no correspondence) and 1 (perfect correspondence).
#' @export
#'
#' @examples
#' A <- list(1:5, 6:10)
#' B <-  purrr::map(1:3, ~ sort(sample(1:10, size = 3, replace = FALSE)))
#' avg_jaccard(A, B)
avg_jaccard <- function(A, B){
  if(!is.list(A) | !is.list(B)) stop('A and B must be lists of vectors.')
  # function to calculate jaccard for two sets
  jaccard <- function(Ai, Bi) {
    return(length(base::intersect(Ai, Bi)) / length(base::union(Ai, Bi)))
  }
  # function to calculate expand.grid fast
  expand_grid <- function(v1, v2) {
    data.frame(A = rep.int(v1, length(v2)),
               B = rep.int(v2, rep.int(length(v1), length(v2))))
  }
  
  # extend to list of sets, averaging
  ## find grid of all possible pairs and calculate jaccard
  grid <- expand_grid(1:length(A), 1:length(B))
  obs_jaccard <-
    purrr::map2_dbl(.x = grid$A, .y = grid$B, ~ jaccard(Ai = A[[.x]], Bi = B[[.y]]))
  grid$jaccard <- obs_jaccard
  ## take the average and return
  cB <-
    base::sum(base::tapply(
      X = grid$jaccard,
      INDEX = grid$B,
      FUN = function(x)
        max(x, na.rm = TRUE)
    )) / (2 * length(B))
  cA <- 
    base::sum(base::tapply(
      X = grid$jaccard,
      INDEX = grid$A,
      FUN = function(x)
        max(x, na.rm = TRUE)
    )) / (2 * length(A))
  
  return(list(avg = cA + cB, grid = grid))
}

