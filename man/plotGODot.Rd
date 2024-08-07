% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEA.R
\name{plotGODot}
\alias{plotGODot}
\title{Visualize GO enrichment test result in dot plot}
\usage{
plotGODot(
  result,
  group = NULL,
  pvalThresh = 0.05,
  n = 20,
  termIDMatch = "^GO",
  colorPalette = "E",
  colorDirection = 1,
  xlab = "-log10(P-value)",
  ylab = "Term name",
  ...
)
}
\arguments{
\item{result}{Returned list object from \code{\link{runGOEnrich}}.}

\item{group}{Character vector of group names, must be available in
\code{names(result)}. Default \code{NULL} make plots for all groups.}

\item{pvalThresh}{Numeric scalar, cutoff for p-value where smaller values are
considered as significant. Default \code{0.05}.}

\item{n}{Number of top terms to be shown, ranked by p-value. Default
\code{20}.}

\item{termIDMatch}{Regular expression pattern to match the term ID. Default
\code{"^GO"} for only using GO terms from returned results.}

\item{colorPalette, colorDirection}{Viridis palette options. Default
\code{"E"} and \code{1}.}

\item{xlab, ylab}{Axis title for x and y axis. Default
\code{"-log10(P-value)"} and \code{"Term name"}, respectively.}

\item{...}{
  Arguments passed on to \code{\link[=.ggplotLigerTheme]{.ggplotLigerTheme}}
  \describe{
    \item{\code{legendColorTitle}}{Legend title text for color aesthetics, often used
for categorical or continuous coloring of dots. Default \code{NULL} shows the
original variable name.}
    \item{\code{legendSizeTitle}}{Legend title text for size aesthetics, often used for
sizing dots by continuous variable. Default \code{NULL} shows the original
variable name.}
    \item{\code{showLegend}}{Whether to show the legend. Default \code{TRUE}.}
    \item{\code{legendPosition}}{Text indicating where to place the legend. Choose from
\code{"top"}, \code{"bottom"}, \code{"left"} or \code{"right"}. Default
\code{"right"}.}
    \item{\code{baseSize}}{One-parameter control of all text sizes. Individual text
element sizes can be controlled by other size arguments. "Title" sizes are
2 points larger than "text" sizes when being controlled by this.}
    \item{\code{titleSize,xTitleSize,yTitleSize,legendTitleSize}}{Size of main title,
axis titles and legend title. Default \code{NULL} controls by
\code{baseSize + 2}.}
    \item{\code{subtitleSize,xTextSize,yTextSize,legendTextSize}}{Size of subtitle text,
axis texts and legend text. Default \code{NULL} controls by \code{baseSize}.}
    \item{\code{plotly}}{Whether to use plotly to enable web based interactive browsing
for the plot. Requires installation of package "plotly". Default
\code{FALSE}.}
  }}
}
\value{
A ggplot object if only one group or a list of ggplot objects.
}
\description{
Visualize GO enrichment test result in dot plot
}
\examples{
defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
# Test the DEG between "stim" and "ctrl", within each cluster
result <- runPairwiseDEG(
    pbmc,
    groupTest = "stim",
    groupCtrl = "ctrl",
    variable1 = "dataset",
    splitBy = "defaultCluster"
)
# Setting `significant = FALSE` because it's hard for a gene list obtained
# from small test dataset to represent real-life biology.
\donttest{
if (requireNamespace("gprofiler2", quietly = TRUE)) {
    go <- runGOEnrich(result, group = "0.stim", significant = FALSE)
    # The toy example won't have significant result.
    plotGODot(go)
}
}
}
