## rliger Next

- Standardized H5 IO specification that can be shared with other platforms.
  - ~~Will move to use HDF5Array (TENxMatrix, H5ADMatrix)/ or BPCells for backed data representation, depending on easiness, cleanliness and future-sustainability of the implementation and cross-platform interoperability.~~
  - Read feature metadata (e.g. id, name, ...) if available; Allow setting "id" as rownames, "name" for visualization.
  - rawData - coming from the original input, read only (qc filtering should be just stored in the object, no IO)
  - ~~preprocessing metrics - nUMI, nGene and etc, still go "chunkApply" so the file is read only once~~
  - ~~normData - delayed computed data from rawData, no on disk representation~~
  - ~~scaleData - new on-disk file and then create object back, because RcppPlanc won't be able to handle delayed computation~~
- Ability to reorganize datasets
  - Allow doing something like `reorganize(ligerObj, variable = "somethingNotDataset")` and resulting in a new liger object with different ligerDataset grouping.
- Ability to do downstream analysis on H5 data
  - Pseudo-bulk should be easy because we are just aggregating cells.
  - Wilcoxon might be a bit harder because ranks are calculated per gene but the 
  H5 sparse data is column majored. Might need to find a fast on-disk 
  transposition method, ~~which would also enhance RcppPlanc performance when 
  running ANLS on H5 data~~.

## rliger 2.2.0

- Implemented highly efficient on-disk iNMF that scales to a million cells using
slightly more time than in-memory version, requiring only laptop-level memory.
- Added 10X H5 data and H5AD loading function that loads the data into regular dgCMatrix 
in memory or the DelayedArray representation backed on disk, the latter is used
for on-disk iNMF implementation.
- Added `selectBatchHVG()` which implements another HVG selection strategy, credit to SCIB
- Adding `suggestK()` back with new methodology
- Clarified optimal `runGOEnrich()` workflow and added fold enrichment metric in 
the returned result
- Fixed important bug in online iNMF scenario 2
- Fixed multiple problems related to ATAC analysis
  - Fixed Wilcoxon rank-sum test bug when using ATAC peak counts
  - Fixed gene coordinate parsing bug from BED file
  - Optimized peak parsing speed

## rliger 2.1.0

- Added `centroidAlign()` for new cell factor loading alignment method
- Added `plotProportionBox()` for visualizing compositional analysis
- Added `plotClusterGeneViolin()` for visualizing gene expression in clusters
- Added `plotBarcodeRank()` for basic QC visualization
- Added `plotPairwiseDEGHeatmap()` for visualizing pairwise DEG results
- Added `plotGODot()` for visualizing GO enrichment results
- Added `calcNMI()` for evaluating clustering results against ground truth
- Added `ligerToH5AD()` allowing reticulate/Python free export of liger object to H5AD format. This is presented in extension source code (i.e. not loaded with `library(rliger)`).
- Added organism support in `runGeneralQC()` and refined hemoglobin gene matching regex pattern.
- Optimized DE test memory usage scalability for both pseudo-bulk method and wilcoxon test
- Optimized `plotProportionPie()` by adding argument `circleColors`
- Optimized `plotVolcano()` text annotation positioning and gene highlighting logic.
- Optimized visualization function additional argument documentation
- Changed `runMarkerDEG()` and `runPairwiseDEG()` default method from `"wilcoxon"` to `"pseudoBulk"`
- Fixed `runMarkerDEG(method = "pseudobulk")` bug in assigning pseudo-replicates, and optimized error/warning signaling.
- Fixed bug in `calcAlignment()`, `subsetMemLigerDataset()`, `cellMeta()`
- Fixed bug in old version updating functions

## rliger 2.0.1

- Fixed wrong UINMF aborting criteria
- Fixed example/test skipping criteria for non-existing dependencies
- Fixed file access issue when checking on CRAN
- Updated installed data file `system.file("extdata/ctrl.h5", "extdata/stim.h5")` to be of standard 10X H5 format
- Updated `quantileNorm()` automatic reference selection according to #297
- Other minor fixes (including #308)

## rliger 2.0.0

- Added `ligerDataset` class for per-dataset information storage, with inheritance for specific modalities
- Added a number of plotting functions with clear function names and useful functionality
- Added Leiden clustering method, now as default rather than Louvain
- Added pseudo-bulk DEG method
- Added DEG analysis with one-vs-rest marker detection in `runMarkerDEG()` and pairwise comparison in `runPairwiseDEG()`
- Added gene name pattern for expression percentage QC metric
- Added native Seurat object support for the core integration workflow
- Added a documentation website built with pkgdown
- Added new iNMF variant method, consensus iNMF (c-iNMF), in `runCINMF()`. Not stable.
- Added GO enrichment dowsntream analysis in `runGOEnrich()`
- Changed `liger` object class structure
- Moved iNMF (previously `optimizeALS()`), UINMF (previously `optimizeALS(unshared = TRUE)`) and online iNMF (previously `online_iNMF()`) implementation to new package *RcppPlanc* with vastly improved performance. Now wrapped in `runINMF()`, `runUINMF()` and `runOnlineINMF()` respectively, and all can be invoked with `runIntegration()`.
- Updated H5AD support to match up with Python anndata package 0.8.0 specs
- Renamed many function/argument names to follow camelCase style, original names are still available while deprecation warnings are issued

## rliger 1.0.1

- Allow setting mito pattern in `getMitoProportion()` #271
- Fix efficiency issue when taking the log of norm.data (e.g. `runWilcoxon`)
- Add runable examples to all exported functions when possible
- Fix typo in online_iNMF matrix initialization
- Adapt to Seurat5
- Other minor fixes

