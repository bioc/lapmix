\name{lap.volcanoplot}

\alias{lap.volcanoplot}

\title{Laplace Volcano Plot}

\description{Creates a volcano plot of log-fold changes versus log-odds of
     differential expression under the Laplace mixture model.
}

\usage{lap.volcanoplot(res, highlight=0, ...)}

\arguments{
  \item{res}{the output from the \code{lapmix.Fit} routine}
  \item{highlight}{number of genes to be highlighted; genes are highlighted in descending order of their posterior odds}
  \item{...}{additional arguments given to the plot function}
}

\value{A plot is created on the current graphics device.}

\details{A volcano plot is any plot which displays fold changes versus a measure of statistical significance of the change.}

\author{Yann Ruffieux}

\examples{
#  See lapmix.Fit example
}

\keyword{hplot}
