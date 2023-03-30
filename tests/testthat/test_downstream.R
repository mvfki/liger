data("pbmc", package = "rliger2")

process <- function(object, f = TRUE, q = TRUE) {
    object <- normalize(object)
    object <- selectGenes(object)
    object <- scaleNotCenter(object)
    if (f) object <- online_iNMF(object, k = 20, miniBatch_size = 100)
    if (q) object <- quantileNorm(object)
    return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Clustering")
test_that("clustering", {
    pbmc <- process(pbmc, f = FALSE, q = FALSE)
    expect_error(runLeidenCluster(pbmc), "No factor loading ")
    expect_error(runLouvainCluster(pbmc), "No factor loading ")

    pbmc <- online_iNMF(pbmc, k = 20, miniBatch_size = 100)
    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1),
                   "Leiden clustering on unnormalized")
    expect_message(runLouvainCluster(pbmc, nRandomStarts = 1),
                   "Louvain clustering on unnormalized")

    pbmc <- quantileNorm(pbmc)
    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1),
                   "Leiden clustering on quantile normalized")
    expect_message(runLouvainCluster(pbmc, nRandomStarts = 1),
                   "Louvain clustering on quantile normalized")

    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1,
                             partitionType = "CPMVertexPartition"),
                   "63 singletons identified. 112 final clusters")
    pbmc <- runLeidenCluster(pbmc, nRandomStarts = 1,
                             partitionType = "CPMVertexPartition",
                             groupSingletons = FALSE)
    expect_equal(nlevels(pbmc$leiden_cluster), 113)
    expect_true("singleton" %in% pbmc$leiden_cluster)

    expect_warning(louvainCluster(pbmc), "was deprecated in")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensionality reduction
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("dimensionality reduction")
test_that("dimensionality reduction", {
    pbmc <- process(pbmc)
    expect_message(runUMAP(pbmc, useRaw = TRUE),
                   "Generating UMAP on unnormalized")
    expect_message(pbmc <- runUMAP(pbmc, useRaw = FALSE),
                   "Generating UMAP on quantile normalized")
    expect_equal(dim(pbmc$UMAP), c(ncol(pbmc), 2))

    expect_message(runTSNE(pbmc, useRaw = TRUE),
                   "Generating TSNE \\(Rtsne\\) on unnormalized")
    expect_message(pbmc <- runTSNE(pbmc, useRaw = FALSE),
                   "Generating TSNE \\(Rtsne\\) on quantile normalized")
    expect_equal(dim(pbmc$TSNE), c(ncol(pbmc), 2))

    expect_error(runTSNE(pbmc, method = "fft"),
                 "Please pass in path to FIt-SNE directory as fitsne.path.")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Differential Expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Differential Expression")
test_that("wilcoxon", {
    expect_error(runWilcoxon(pbmc),
                 "All datasets being involved has to be normalized")
    pbmc <- process(pbmc)
    pbmc <- runLeidenCluster(pbmc, nRandomStarts = 1)
    expect_error(runWilcoxon(pbmc, method = "dataset", useDatasets = 1),
                 "Should have at least 2 datasets as input ")
    res1 <- runWilcoxon(pbmc)
    expect_equal(dim(res1), c(1992, 10))
    expect_equal(res1[1,4], -3.6828172, tolerance = 1e-6)
    res2 <- runWilcoxon(pbmc, method = "dataset")
    expect_equal(dim(res2), c(3984, 10))
    expect_equal(res2[1,7], 2.936397e-24, tolerance = 1e-6)
    hm1 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "l")
    hm2 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "p")
    expect_is(hm1, "HeatmapList")
    expect_is(hm2, "HeatmapList")
    expect_is(plotVolcano(res1, 0), "ggplot")
    expect_is(plotEnhancedVolcano(res1, 0), "ggplot")

    expect_error(getFactorMarkers(pbmc, "ctrl", "stim", factorShareThresh = 0),
                 "No factor passed the dataset specificity threshold")
    expect_warning(
        getFactorMarkers(pbmc, "ctrl", "stim", factorShareThresh = 0.01),
        "Only 1 factor passed the dataset specificity threshold"
    )
    expect_warning(
        expect_message(
            getFactorMarkers(pbmc, "ctrl", "stim", printGenes = TRUE),
            "GAPDH, LGALS1, CXCR4, ACTB, FTL, ISG15, GBP1, SELL, RSAD2, TEX264"
        ),
        "Factor 16 did not appear as max in any cell in either dataset"
    )
    expect_warning(res3 <- getFactorMarkers(pbmc, "ctrl", "stim"),
                   "Factor 16 did not appear as max in any cell in either")
    expect_is(res3, "list")
    expect_identical(names(res3), c("ctrl", "shared", "stim", "num_factors_V1",
                                   "num_factors_V2"))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GSEA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("GSEA")
custom <- list(
    `Immune System` = c("9636", "2633", "6282", "6280", "6279", "2207", "2214",
                        "6402", "91543", "6233", "10578", "3553", "5473",
                        "3627", "51316", "929", "972")
)
test_that("gsea", {
    expect_warning({
        expect_is(runGSEA(pbmcPlot, genesets = "Immune System"), "list")
        expect_is(runGSEA(pbmcPlot, customGenesets = custom), "list")
    })

    res1 <- runFactorGeneGO(pbmcPlot)
    expect_identical(unique(res1$result$query), c("Factor_13", "Factor_19"))
    res2 <- runFactorGeneGO(pbmcPlot, sumLoading = FALSE)
    expect_identical(unique(res2$result$query),
                     c("ctrl_Factor_13", "ctrl_Factor_19",
                       "stim_Factor_13", "stim_Factor_19"))
})