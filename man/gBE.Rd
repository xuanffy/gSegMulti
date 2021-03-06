% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gBE.R
\name{gBE}
\alias{gBE}
\title{Change-Point selection by graph-based backward elimination}
\usage{
gBE(
  y,
  tauhat,
  c = 2,
  graph = "mst",
  kMax = 5,
  distType = "euclidean",
  minCpNum = 1,
  maxSeq = 5000,
  detail = FALSE
)
}
\arguments{
\item{y}{A n x d matrix to be scanned with n observations and d dimensions.}

\item{tauhat}{A vector of candidate change-points to be selected.}

\item{c}{Penalty multiplier.}

\item{graph}{Type of similarity graphs.\preformatted{"mst" specifies the minimum spanning tree;

"knn" specifies the nearest neighbor graph.
}}

\item{kMax}{Max k of the similarity graph (k-mst or knn).}

\item{distType}{The distance measure to be used in the \code{\link[stats]{dist}} function.}

\item{minCpNum}{Minimum number of change-points.}

\item{maxSeq}{Maximum length of sequence to be tested.}

\item{detail}{Return detailed information.}
}
\value{
If \code{detail = FALSE}, the function returns selected change-points. If \code{detail = TRUE}, the function returns a list.
\item{tauhat}{Selected change-points.}
\item{gofSeq}{Goodness-of-fit statistic values in each step of backward elimination.}
\item{mergeSeq}{Order of candidate change-points removed in backward elimination.}
}
\description{
This function prunes candidate change-points and choose the model with the largest exteded pseudo-BIC values.
}
\examples{
set.seed(1)
y <- rbind(
  matrix(rnorm(200), ncol = 10), matrix(rexp(200), ncol = 10),
  matrix(rnorm(200), ncol = 10), matrix(rexp(200), ncol = 10)
)
gBE(y, c(20, 30, 40, 60), c = 2, kMax = 2, detail = TRUE)
}
