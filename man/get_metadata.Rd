% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_circlize.R
\name{get_metadata}
\alias{get_metadata}
\title{Get metadata from a Seurat object}
\usage{
get_metadata(seu_obj, reductions = "umap", coord_scale = 0.8, color)
}
\arguments{
\item{seu_obj}{SeuratObject}

\item{reductions}{reductions methods, e.g."umap" or "tsne".}

\item{coord_scale}{value from c(0,1) to adjust the UMAP/tSNE coordinates.}

\item{color}{Colors assigned to the cell clusters}
}
\value{
A metadata dataframe
}
\description{
This function extracts the metadata from a Seurat object and transforms the
UMAP/tSNE coordinates.
}
