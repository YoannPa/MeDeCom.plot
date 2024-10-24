% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/medecom_plots.R
\name{mds_plot}
\alias{mds_plot}
\title{Plots a MDS of LMCs and samples.}
\usage{
mds_plot(MDCset, D, k, lambda, cor.method = "pearson", plot.title = NULL)
}
\arguments{
\item{MDCset}{A \code{MeDeComSet} obtained with \link[MeDeCom]{runMeDeCom}.}

\item{D}{A methylation \code{matrix} you ran MeDeCom on.}

\item{k}{Am \code{integer} to specify the value of the kappa
parameter.}

\item{lambda}{A \code{numeric} to specify the value of the lambda
parameter.}

\item{cor.method}{A \code{character} to specify the method to use for
computing correlation between LMCs and samples
(Default: cor.method = 'pearson'; For a list of all
supported correlation methods see \link[stats]{cor}).}

\item{plot.title}{A \code{character} to specify another title than the
default one for your plot.}
}
\value{
A \code{gg} plot of the MDS on LMCs and samples.
}
\description{
Plots a MDS of LMCs and samples.
}
\examples{
# Load the example data sets
data(example.dataset, package="MeDeCom")
# Run MeDeCom (WARNING: takes 1 or 2 hours to complete)
medecom.result <- runMeDeCom(
  D, 2:10, c(0,10^(-5:-2)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=9)
# Plot Pearson-based MDS
mds_plot(MDCset = medecom.result, D = D, k = 5, lambda = 10^-2)
# Plot Spearman-based MDS
mds_plot(
  MDCset = medecom.result, D = D, k = 5, lambda = 10^-2,
  cor.method = "spearman")
}
\references{
Lutsik P. et al., MeDeCom: discovery and quantification of latent
components of heterogeneous methylomes. Genome Biol. 2017 Mar
24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624;
PMCID: PMC5366155.
}
\author{
Yoann Pageaud.
}
