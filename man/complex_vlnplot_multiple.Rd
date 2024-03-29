% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_violin.R
\name{complex_vlnplot_multiple}
\alias{complex_vlnplot_multiple}
\title{Violin plot for multiple genes across groups}
\usage{
complex_vlnplot_multiple(
  seu_obj,
  features,
  celltypes = NULL,
  group,
  add.dot = T,
  font.size = 12,
  pt.size = 0.1,
  alpha = 0.01,
  strip.color = NULL
)
}
\arguments{
\item{seu_obj}{A complete Seurat object}

\item{features}{Gene name. Only one gene is allowed.}

\item{celltypes}{Cell types of interest. By default, all cell types are included.}

\item{group}{Only one groupID is allowed.}

\item{add.dot}{Whether or not to add points on the violins.}

\item{font.size}{Font size for the labels.}

\item{pt.size}{Point size for the data points on the violin}

\item{alpha}{Point transparency. value from 0 to 1.}

\item{strip.color}{Colors for the strip background}
}
\value{
A ggplot object
}
\description{
This function generates violin plot(s) to compare the expression of multiple genes across
different groups or cell types. It is designed for visualizing a complicated scenario:
Gene expression of multiple genes on multiple cell types across groups.
}
