% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Wham.vis.R
\name{FeatureCorr}
\alias{FeatureCorr}
\title{feature correlation}
\usage{
FeatureCorr(
  object,
  feature,
  group,
  UseDE.result = F,
  p.cutoff = 0.05,
  effect.size.range = c(0, 0),
  method = c("pearson", "spearman", "kendall"),
  adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  count = c("rawCount", "WhamTransformation")
)
}
\arguments{
\item{object}{required. The output of WhamBiobakery, WhamEBI or Wham16s}

\item{feature}{default is NULL. if not provided and UseDE.result is FALSE, all features will be used.}

\item{group}{default is NULL. if not provided and UseDE.result is FALSE, all samples will be used}

\item{UseDE.result}{default value is FALSE. FALSE means using the argument feature and argument group. TRUE means using the features and samples from DE.result (after selection of parameters below).}

\item{p.cutoff}{Keep genes whose we.eHB values are less than cutoff. The default value is 0.05. Use all gene from DE.result}

\item{effect.size.range}{Keep genes whose effect size values are out of the range.The default value is c(0,0), Use all gene from DE.result. This only work when the aldex.module is set as t.test}

\item{method}{the method of correlation. "spearman","pearson", and "kendall"}

\item{adjust}{the method of multiple test correction. "holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"}

\item{count}{c("rawCount","WhamTransformation"). Default value is "rawCount". Setting as "rawCount" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)}
}
\value{
a correlation heatmap of gene family
}
\description{
Visualization of Wham upstream analysis
}
\examples{

# see the vignette for more details
}
