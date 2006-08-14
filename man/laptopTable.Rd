\name{laptopTable}

\alias{laptopTable}

\title{Table of Top Genes from Laplace Linear Model Fit}

\description{Extract a table of the top-ranked genes from a Laplace mixture model fit.}

\usage{laptopTable(res, number=res$med.number, sort.by='L')}

\arguments{
  \item{res}{the output from the \code{lapmix.Fit} routine}
  \item{number}{how many genes to pick out; if missing: number is determined by the median rule}
  \item{sort.by}{character string specifying statistic to sort the selected genes by in the output data.frame.}
}

\value{
  A dataframe with a row for the \code{number} top genes and the following columns:
  \item{M}{average log fold change}
  \item{log.odds}{log posterior odds that the gene is differentially expressed}
}

\details{
This function summarizes a Laplace mixture model fit object produced by \code{lapmix.Fit} by selecting the top-ranked genes according to the posterior log-odds or M-values. 

The \code{sort.by} argument specifies the criterion used to select the top genes. Only two choices at the moment: \code{"M"} to sort by the (absolute) coefficient representing the log-fold-change, and \code{"L"} to sort by the posterior odds of differential expression under the Laplace mixture model.
}

\author{Yann Ruffieux}

\examples{
#  See lapmix.Fit example
}

\keyword{htest}
