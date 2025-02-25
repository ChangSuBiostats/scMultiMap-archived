% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{make_gene_to_peak_link}
\alias{make_gene_to_peak_link}
\title{make_gene_to_peak_link}
\usage{
make_gene_to_peak_link(df, gene, obj, peak_assay = "peak")
}
\arguments{
\item{df}{A data frame of peak-gene associations results, containing columns of \code{gene}, \code{peak}, and pvar (see below).}

\item{gene}{Character, gene of interest.}

\item{obj}{Seurat object that contains the data on the gene and the peaks.}

\item{peak_assay}{Name for peak assay. Default to 'peak'.}
}
\value{
A grange object with genomic ranges corresponding to the gene and the peak, and with the following metadata:
\describe{
\item{score}{correlation from the peak-gene association results}
\item{gene}{gene name}
\item{peak}{peak name}
\item{test_stat}{test statistics from the peak-gene association results`}
\item{pval}{p-value from the peak-gene association results}
}
}
\description{
make a data frame that links peak to a given gene as input to Signac::LinkPlot()
}
