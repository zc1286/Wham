% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wham.vis.R
\name{TaxaBarPlot}
\alias{TaxaBarPlot}
\title{TaxaBarPlot}
\usage{
TaxaBarPlot(
  object,
  taxa.level = c("k", "p", "o", "c", "f", "g", "s", "otu"),
  filter = c("all", "custom", "DE.filter"),
  display.number = 30,
  custom,
  p.cutoff = 0.05,
  effect.size.range,
  annotation,
  merge_group = F,
  abundance_cutoff = 0.001,
  relative_abundance = c("to_selection", "to_total")
)
}
\arguments{
\item{object}{required. The output of WhamBiobakery, WhamEBI or Wham16s}

\item{taxa.level}{Required, c("k","p","o","c","f","g","s").Collapse taxa at the provided level.}

\item{filter}{c("all","custom","DE.filter"). Default is "all". "all" shows all samples and top # taxa(indicated in display.number value (default:30)). "custom" shows taxa provided by user and all samples. "DE.filter" shows taxa which pass DE test cutoff(p.cutoff and effect.size.range) and samples which are in the DE.result contrast.}

\item{display.number}{applied when filter is set as "all". Default is 30.}

\item{custom}{required when filter is set as "custom". A vector which contain taxa to be shown}

\item{p.cutoff}{applied when filter is set as "DE.filter". Default is 0.05.}

\item{effect.size.range}{applied when filter is set as "DE.filter" and aldex.module is set as t.test.Default range is c(0,0)}

\item{annotation}{optional. Default is the last variable in design. A character vector indicates what column of metadata will be used for annotation bar}

\item{merge_group}{Default F. Merging samples under the condition.(mean of relative abundance)}

\item{abundance_cutoff}{When setting filter argument as "all", the bacterium whose the proportion is less than abundance_cutoff will be removed.}

\item{relative_abundance}{optional. Default is "to_selection". "to_selection" means counts are divided by the total sum of selected taxonomic component in each library. "to_total" means counts are divided by the total sum of all components in each library}
}
\value{
a barplot of taxonomy
}
\description{
Visualization of Wham upstream analysis
}
\examples{
# see the vignette for more details
}
