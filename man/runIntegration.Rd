% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integration.R
\name{runIntegration}
\alias{runIntegration}
\alias{runIntegration.liger}
\alias{runIntegration.Seurat}
\title{Integrate scaled datasets with iNMF or variant methods}
\usage{
runIntegration(
  object,
  k = 20,
  lambda = 5,
  method = c("iNMF", "onlineINMF", "UINMF"),
  ...
)

\method{runIntegration}{liger}(
  object,
  k = 20,
  lambda = 5,
  method = c("iNMF", "onlineINMF", "UINMF"),
  seed = 1,
  verbose = getOption("ligerVerbose", TRUE),
  ...
)

\method{runIntegration}{Seurat}(
  object,
  k = 20,
  lambda = 5,
  method = c("iNMF", "onlineINMF"),
  datasetVar = "orig.ident",
  useLayer = "ligerScaleData",
  assay = NULL,
  seed = 1,
  verbose = getOption("ligerVerbose", TRUE),
  ...
)
}
\arguments{
\item{object}{A \linkS4class{liger} object or a Seurat object with
non-negative scaled data of variable features (Done with
\code{\link{scaleNotCenter}}).}

\item{k}{Inner dimension of factorization (number of factors). Generally, a
higher \code{k} will be needed for datasets with more sub-structure. Default
\code{20}.}

\item{lambda}{Regularization parameter. Larger values penalize
dataset-specific effects more strongly (i.e. alignment should increase as
\code{lambda} increases). Default \code{5}.}

\item{method}{iNMF variant algorithm to use for integration. Choose from
\code{"iNMF"}, \code{"onlineINMF"}, \code{"UINMF"}. Default \code{"iNMF"}.}

\item{...}{Arguments passed to other methods and wrapped functions.}

\item{seed}{Random seed to allow reproducible results. Default \code{1}.}

\item{verbose}{Logical. Whether to show information of the progress. Default
\code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.}

\item{datasetVar}{Metadata variable name that stores the dataset source
annotation. Default \code{"orig.ident"}.}

\item{useLayer}{For Seurat>=4.9.9, the name of layer to retrieve input
non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
always retrieve from \code{scale.data} slot.}

\item{assay}{Name of assay to use. Default \code{NULL} uses current active
assay.}
}
\value{
Updated input object. For detail, please refer to the refered method
linked in Description.
}
\description{
LIGER provides dataset integration methods based on iNMF (integrative
Non-negative Matrix Factorization [1]) and its variants (online iNMF [2]
and UINMF [3]). This function wraps \code{\link{runINMF}},
\code{\link{runOnlineINMF}} and \code{\link{runUINMF}}, of which the help
pages have more detailed description.
}
\examples{
pbmc <- normalize(pbmc)
pbmc <- selectGenes(pbmc)
pbmc <- scaleNotCenter(pbmc)
if (requireNamespace("RcppPlanc", quietly = TRUE)) {
    pbmc <- runIntegration(pbmc)
}
}
\references{
\enumerate{
\item{Joshua D. Welch and et al., Single-Cell Multi-omic Integration Compares
and Contrasts Features of Brain Cell Identity, Cell, 2019}
\item{Chao Gao and et al., Iterative single-cell multi-omic integration using
online learning, Nat Biotechnol., 2021}
\item{April R. Kriebel and Joshua D. Welch, UINMF performs mosaic integration
of single-cell multi-omic datasets using nonnegative matrix factorization,
Nat. Comm., 2022}
}
}
