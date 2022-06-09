#' Graph-based multiple change-point detection
#' @description This function performs gMulti: multiple change-point detection by graph-based method. The whole process involves 1. searching and 2. pruning step. A change-point dedrogram will be plotted for understanding structure of detected change-points.
#' @param y A n x d matrix to be scanned with n observations and d dimensions.
#' @param search.type Greedy method to be performed in the first step.
#'
#'      "WBS" specifies graph-based wild binary segmentation
#'
#'      "SBS" specifies graph-based seeded binary segmentation

#' @param model.selection Perform the pruning (model selection) step.
#' @param dendrogram Plot change-point dendrogram
#' @param ... Arguments to be passed to \code{\link[gMulti]{gWBS}}, \code{\link[gMulti]{gSBS}}, \code{\link[gMulti]{gBE}} and \code{\link[gMulti]{cpdendrogram}}
#'
#' @return
#' \item{tauhat}{Final selected change-points after pruning.}
#' \item{tautilde}{Searched change-points in the first step.}
#'
#' @examples
#' set.seed(1)
#' y <- rbind(
#'   matrix(rnorm(200), ncol = 10), matrix(rexp(200), ncol = 10),
#'   matrix(rnorm(200), ncol = 10), matrix(rexp(200), ncol = 10)
#' )
#' gMulti(y)
#' @export
gMulti <- function(y, search.type = "WBS", model.selection = TRUE, dendrogram = TRUE, ...) {
  if (search.type == "WBS") {
    searchr <- gWBS(y = y, ...)
  } else if (search.type == "SBS") {
    searchr <- gSBS(y = y, ...)
  } else {
    stop(paste('No method called "', search.type, '". Use "WBS" or "SBS" instead.', sep = ""))
  }
  result <- list(tautilde = searchr)
  if (model.selection) {
    if (min(diff(sort(searchr))) <= 20) {
      prune <- gBE(y, searchr, kMax = 3, detail = TRUE, ...)
    } else {
      prune <- gBE(y, searchr, detail = TRUE, ...)
    }
    if (dendrogram) {
      if (length(prune$tauhat) > 1) {
        cpdendrogram(y, prune)
      } else {
        warning("only one change-point detected, cannot plot a dendrogram")
      }
    }
    result$tauhat <- prune$tauhat
  }
  return(result)
}
