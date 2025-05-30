% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{scaleNotCenter}
\alias{scaleNotCenter}
\alias{scaleNotCenter.dgCMatrix}
\alias{scaleNotCenter.DelayedArray}
\alias{scaleNotCenter.ligerDataset}
\alias{scaleNotCenter.ligerMethDataset}
\alias{scaleNotCenter.liger}
\alias{scaleNotCenter.Seurat}
\title{Scale genes by root-mean-square across cells}
\usage{
scaleNotCenter(object, ...)

\method{scaleNotCenter}{dgCMatrix}(object, features, scaleFactor = NULL, ...)

\method{scaleNotCenter}{DelayedArray}(
  object,
  features,
  scaleFactor = NULL,
  geneRootMeanSq = NULL,
  overwrite = FALSE,
  chunk = getOption("ligerChunkSize", 20000),
  verbose = getOption("ligerVerbose", TRUE),
  ...
)

\method{scaleNotCenter}{ligerDataset}(
  object,
  features = NULL,
  scaleFactor = NULL,
  chunk = getOption("ligerChunkSize", 20000),
  verbose = getOption("ligerVerbose", TRUE),
  ...
)

\method{scaleNotCenter}{ligerMethDataset}(
  object,
  features = NULL,
  verbose = getOption("ligerVerbose", TRUE),
  ...
)

\method{scaleNotCenter}{liger}(
  object,
  useDatasets = NULL,
  features = varFeatures(object),
  verbose = getOption("ligerVerbose", TRUE),
  remove.missing = NULL,
  ...
)

\method{scaleNotCenter}{Seurat}(
  object,
  assay = NULL,
  layer = "ligerNormData",
  save = "ligerScaleData",
  datasetVar = "orig.ident",
  features = NULL,
  ...
)
}
\arguments{
\item{object}{\linkS4class{liger} object, \linkS4class{ligerDataset} object,
\link[Matrix]{dgCMatrix-class} object, or a Seurat object.}

\item{...}{Arguments passed to other methods. The order goes by: "liger"
method calls "ligerDataset" method", which then calls "dgCMatrix" method.
"Seurat" method directly calls "dgCMatrix" method.}

\item{features}{Character, numeric or logical index that choose the variable
feature to be scaled. "liger" method by default uses
\code{\link{varFeatures}(object)}. "ligerDataset" method by default uses all
features. "Seurat" method by default uses
\code{Seurat::VariableFeatures(object)}.}

\item{scaleFactor}{Numeric vector of scaling factor to normalize the raw
counts to unit sum. This pre-calculated at liger object creation (stored as
\code{object$nUMI} and internally specified in S3 method chains, thus is
generally not needed to be specified by users.}

\item{geneRootMeanSq}{Numeric vector of root-mean-square of unit-normalized
expression for each gene. This is pre-calculated at the call of
\code{\link{selectBatchHVG}} (stored at
\code{featureMeta(dataset(object, "datasetName"))$rootMeanSq} and internally
specified in S3 method chains, thus is generally not needed to be specified
by users.}

\item{overwrite}{Logical. When writing newly computed HDF5 array to a
separate HDF5 file, whether to overwrite the existing file. Default
\code{FALSE} raises an error when the file already exists.}

\item{chunk}{Integer. Number of maximum number of cells in each chunk, when
scaling is applied to any HDF5 based dataset. Default \code{20000}.}

\item{verbose}{Logical. Whether to show information of the progress. Default
\code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.}

\item{useDatasets}{A character vector of the names, a numeric or logical
vector of the index of the datasets to be scaled but not centered. Default
\code{NULL} applies to all datasets.}

\item{remove.missing}{\bold{Deprecated}. The functionality of this is covered
through other parts of the whole workflow and is no long needed. Will be
ignored if specified.}

\item{assay}{Name of assay to use. Default \code{NULL} uses current active
assay.}

\item{layer}{For Seurat>=4.9.9, the name of layer to retrieve normalized
data. Default \code{"ligerNormData"}. For older Seurat, always retrieve from
\code{data} slot.}

\item{save}{For Seurat>=4.9.9, the name of layer to store normalized data.
Default \code{"ligerScaleData"}. For older Seurat, stored to
\code{scale.data} slot.}

\item{datasetVar}{Metadata variable name that stores the dataset source
annotation. Default \code{"orig.ident"}.}
}
\value{
Updated \code{object}
\itemize{
\item{dgCMatrix method - Returns scaled dgCMatrix object}
\item{ligerDataset method - Updates the \code{scaleData} and
\code{scaledUnsharedData} (if unshared variable feature available) slot
of the object}
\item{liger method - Updates the \code{scaleData} and
\code{scaledUnsharedData} (if unshared variable feature available) slot
of chosen datasets}
\item{Seurat method - Adds a named layer in chosen assay (V5), or update the
\code{scale.data} slot of the chosen assay (<=V4)}
}
}
\description{
This function scales normalized gene expression data after
variable genes have been selected. We do not mean-center the data before
scaling in order to address the non-negativity constraint of NMF.
Computation applied to each normalized dataset matrix can form the following
equation:

\deqn{S_{i,j}=\frac{N_{i,j}}{\sqrt{\sum_{p}^{n}\frac{N_{i,p}^2}{n-1}}}}

Where \eqn{N} denotes the normalized matrix for an individual dataset,
\eqn{S} is the output scaled matrix for this dataset, and \eqn{n} is the
number of cells in this dataset. \eqn{i, j} denotes the specific gene and
cell index, and \eqn{p} is the cell iterator.

Please see detailed section below for explanation on methylation dataset.
}
\note{
Since the scaling on genes is applied on a per dataset base, other scaling
methods that apply to a whole concatenated matrix of multiple datasets might
not be considered as equivalent alternatives, even if options like
\code{center} are set to \code{FALSE}. Hence we implemented an efficient
solution that works under such circumstance, provided with the Seurat S3
method.
}
\section{Methylation dataset}{

Because gene body mCH proportions are negatively correlated with gene
expression level in neurons, we need to reverse the direction of the
methylation data before performing the integration. We do this by simply
subtracting all values from the maximum methylation value. The resulting
values are positively correlated with gene expression. This will only be
applied to variable genes detected in prior. Please make sure that argument
\code{modal} is set accordingly when running \code{\link{createLiger}}. In
this way, this function can automatically detect it and take proper action.
If it is not set, users can still manually have the equivalent processing
done by doing \code{scaleNotCenter(lig, useDataset = c("other", "datasets"))},
and then \code{\link{reverseMethData}(lig, useDataset = c("meth", "datasets"))}.
}

\examples{
pbmc <- selectBatchHVG(pbmc, n = 10)
pbmc <- scaleNotCenter(pbmc)
}
