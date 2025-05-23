# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename) {
    .Call(`_rliger_RunModularityClusteringCpp`, SNN, modularityFunction, resolution, algorithm, nRandomStarts, nIterations, randomSeed, printOutput, edgefilename)
}

moe_correct_ridge_cpp <- function(Z_orig, R, lambda, Phi, B, N) {
    .Call(`_rliger_moe_correct_ridge_cpp`, Z_orig, R, lambda, Phi, B, N)
}

normalize_byCol_dense_rcpp <- function(x) {
    .Call(`_rliger_normalize_byCol_dense_rcpp`, x)
}

colNormalize_dense_cpp <- function(x, L) {
    .Call(`_rliger_colNormalize_dense_cpp`, x, L)
}

scaleNotCenter_byRow_rcpp <- function(x) {
    .Call(`_rliger_scaleNotCenter_byRow_rcpp`, x)
}

safe_scale <- function(x, center, scale) {
    .Call(`_rliger_safe_scale`, x, center, scale)
}

scaleNotCenter_byCol_dense_rcpp <- function(x) {
    .Call(`_rliger_scaleNotCenter_byCol_dense_rcpp`, x)
}

scaleNotCenter_byRow_perDataset_rcpp <- function(x, ann, n) {
    .Call(`_rliger_scaleNotCenter_byRow_perDataset_rcpp`, x, ann, n)
}

rowVars_sparse_rcpp <- function(x, means, ncol) {
    .Call(`_rliger_rowVars_sparse_rcpp`, x, means, ncol)
}

rowDivide_rcpp <- function(x, v) {
    .Call(`_rliger_rowDivide_rcpp`, x, v)
}

sumSquaredDeviations <- function(x, means) {
    .Call(`_rliger_sumSquaredDeviations`, x, means)
}

denseZScore <- function(x, m) {
    .Call(`_rliger_denseZScore`, x, m)
}

rowVarsDense <- function(x, m) {
    .Call(`_rliger_rowVarsDense`, x, m)
}

SparseRowVarStd <- function(x, mu, sd, vmax) {
    .Call(`_rliger_SparseRowVarStd`, x, mu, sd, vmax)
}

colAggregateSums_sparse <- function(x, group, ngroups) {
    .Call(`_rliger_colAggregateSums_sparse`, x, group, ngroups)
}

colAggregateMedian_dense_cpp <- function(x, group, n) {
    .Call(`_rliger_colAggregateMedian_dense_cpp`, x, group, n)
}

sample_cpp <- function(x, size) {
    .Call(`_rliger_sample_cpp`, x, size)
}

updatePseudoBulkRcpp <- function(psdBulk, sparseRaw, featureIdx, repIdx) {
    invisible(.Call(`_rliger_updatePseudoBulkRcpp`, psdBulk, sparseRaw, featureIdx, repIdx))
}

updateNCellExprRcpp <- function(out, sparseRaw, featureIdx, groupVar) {
    invisible(.Call(`_rliger_updateNCellExprRcpp`, out, sparseRaw, featureIdx, groupVar))
}

#' Fast calculation of feature count matrix
#'
#' @param bedmat A feature count list generated by bedmap
#' @param barcodes A list of barcodes
#'
#' @return A feature count matrix with features as rows and barcodes as
#' columns
#' @export
#' @examples
#' \dontrun{
#' gene.counts <- makeFeatureMatrix(genes.bc, barcodes)
#' promoter.counts <- makeFeatureMatrix(promoters.bc, barcodes)
#' samnple <- gene.counts + promoter.counts
#' }
makeFeatureMatrix <- function(bedmat, barcodes) {
    .Call(`_rliger_makeFeatureMatrix`, bedmat, barcodes)
}

objErr_i <- function(H, W, V, E, lambda) {
    .Call(`_rliger_objErr_i`, H, W, V, E, lambda)
}

cluster_vote_rcpp <- function(nn_ranked, clusts) {
    .Call(`_rliger_cluster_vote_rcpp`, nn_ranked, clusts)
}

max_factor_rcpp <- function(H, dims_use, center = FALSE) {
    .Call(`_rliger_max_factor_rcpp`, H, dims_use, center)
}

ComputeSNN <- function(nn_idx, prune) {
    .Call(`_rliger_ComputeSNN`, nn_idx, prune)
}

WriteEdgeFile <- function(snn, filename, display_progress) {
    invisible(.Call(`_rliger_WriteEdgeFile`, snn, filename, display_progress))
}

DirectSNNToFile <- function(nn_ranked, prune, display_progress, filename) {
    .Call(`_rliger_DirectSNNToFile`, nn_ranked, prune, display_progress, filename)
}

cpp_rank_matrix_dgc <- function(x, p, nrow, ncol, showProgress = FALSE) {
    .Call(`_rliger_cpp_rank_matrix_dgc`, x, p, nrow, ncol, showProgress)
}

rowAggregateSum_sparse <- function(X, groups, ngroups) {
    .Call(`_rliger_rowAggregateSum_sparse`, X, groups, ngroups)
}

colAggregateSum_sparse <- function(X, groups, ngroups) {
    .Call(`_rliger_colAggregateSum_sparse`, X, groups, ngroups)
}

colNNZAggr_sparse <- function(X, groups, ngroups) {
    .Call(`_rliger_colNNZAggr_sparse`, X, groups, ngroups)
}

rowNNZAggr_sparse <- function(X, groups, ngroups) {
    .Call(`_rliger_rowNNZAggr_sparse`, X, groups, ngroups)
}

