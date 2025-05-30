% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downsample.R
\name{readSubset}
\alias{readSubset}
\title{\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#superseded}{\figure{lifecycle-superseded.svg}{options: alt='[Superseded]'}}}{\strong{[Superseded]}} See \code{\link{downsample}}}
\usage{
readSubset(
  object,
  slot.use = "normData",
  balance = NULL,
  max.cells = 1000,
  chunk = 1000,
  datasets.use = NULL,
  genes.use = NULL,
  rand.seed = 1,
  verbose = getOption("ligerVerbose", TRUE)
)
}
\arguments{
\item{object}{\linkS4class{liger} object}

\item{slot.use}{Only create subset from one or more of \code{"rawData"},
\code{"normData"} and \code{"scaleData"}. Default \code{NULL} subsets the
whole \code{object} including downstream results.}

\item{balance}{\code{"all"} for sampling \code{maxCells} cells from all
datasets specified by \code{useDatasets}. \code{"cluster"} for sampling
\code{maxCells} cells per cluster per dataset. \code{"dataset"} for
\code{maxCells} cells per dataset.}

\item{max.cells}{Max number of cells to sample from the grouping based on
\code{balance}.}

\item{chunk}{Integer. Number of maximum number of cells in each chunk,
Default \code{1000}.}

\item{datasets.use}{Index selection of datasets to consider. Default
\code{NULL} for using all datasets.}

\item{genes.use}{Character vector. Subset features to this specified range.
Default \code{NULL} does not subset features.}

\item{rand.seed}{Random seed for reproducibility. Default \code{1}.}

\item{verbose}{Logical. Whether to show information of the progress. Default
\code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.}
}
\value{
Subset of \linkS4class{liger} \code{object}.
}
\description{
This function mainly aims at downsampling datasets to a size
suitable for plotting.
}
\seealso{
\code{\link{downsample}}, \code{\link{subsetLiger}},
\code{\link{subsetLigerDataset}}
}
