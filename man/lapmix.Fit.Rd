\name{lapmix.Fit}

\alias{lapmix.Fit}

\title{Empirical Bayes Statistics for Differential Expression under Laplace Model}

\description{Computes posterior odds of differential expression under the Laplace mixture model, with parameters estimated using an empirical Bayes approach.}

\usage{
lapmix.Fit(Y, asym=FALSE, fast=TRUE, two.step=TRUE, 
           w=0.1, V=10, beta=0, gamma=1, alpha=0.1)
}

\arguments{
  \item{Y}{data.frame or matrix containing the log relative expression levels, where each row represents a gene. Alternatively, a list of arrays of possibly different sizes, or an object of class \code{eSet}, \code{exprSet}, or \code{ExpressionSet}.}
  \item{asym}{indicates whether the asymmetric Laplace model is used rather than the symmetric one.}
  \item{fast}{indicates whether the 'fast' estimation method is used.}
  \item{two.step}{indicates whether the two-step estimation method is used; otherwise the marginal likelihood is maximised in one step.}
  \item{w,V,beta,gamma,alpha}{initial values given to the optimization algorithms for the estimation of the hyperparameters.}
}

\value{
  \item{lap.lods}{numeric vector containing the posterior log-odds of differential expression}
  \item{prob}{numeric vector containing the posterior probabilities of differential expression}
  \item{med.number}{number of differentially expressed genes according to the median rule}
  \item{M}{numeric vector with average log fold changes within genes}
  \item{s_sq}{numeric vector with sample variances within genes}
  \item{nb.rep}{numeric vector with number of replicates within genes}
  \item{estimates}{list containing the empirical Bayes estimates of the hyperparameters}
  \item{code}{integer indicating why the likelihood optmization terminated, cf. \code{nlm} routine. There are two such indicators in the case of the two-step estimation.}
}

\details{This method fits the results of a microarray experiment to a Laplace mixture model. These results are assumed to take the form of normalized base 2 logarithm of the expression ratios. An empirical Bayes approach is used to estimate the hyperparameters of the model. The \code{lap.lodds} is sometimes known as the L-statistic (if the symmetric model is used) or the AL-statistic (if the asymmetric model is used). These statistics can be used to rank the genes according to the posterior odds of differential expression, via the routine \code{laptopTable}. They can be visualized using the \code{lap.volcanoplot} function.

If there are different numbers of replicates between genes, one may wish to write the data in a list of arrays. If a matrix representation is desired, one can stick in NaN's where appropriate.

The `fast' estimation method ignores the integrals which cannot be computed with the t-distribution function. This method is suggested, since these problematic integrals are few and far between. The estimates are practically not affected, and we avoid the potential problems that arise when integrating numerically with the \code{integrate} function.}

\author{Yann Ruffieux}

\references{Bhowmick, D., Davison, A.C., and Goldstein, D.R. (2006). A Laplace mixture model for identification of differential expression in microarray experiments.}

\examples{

#  Simulate gene expression data under Laplace mixture model: 3000 genes with
#  4 duplicates each; one gene in ten is differentially expressed.

G <- 3000
Y <- NULL
sigma_sq <- 1/rgamma(G, shape=2.8, scale=0.04)
mu <- rexp(G, rate=1/(sigma_sq*1.2))-rexp(G, rate=1/(sigma_sq*1.2))
is.diff <- sample(c(0,1), replace=TRUE, prob=c(0.9,0.1), size=G)
mu <- mu*is.diff
for(g in 1:G)
	Y <- rbind(Y, rnorm(4,mu[g], sd=sqrt(sigma_sq[g])))

# with symmetric model
res <- lapmix.Fit(Y)
res$estimates
laptopTable(res, 20)
lap.volcanoplot(res, highlight=res$med.number)

# with asymmetric model
res2 <- lapmix.Fit(Y, asym=TRUE)
res2$estimates
laptopTable(res2, 20)
lap.volcanoplot(res2, highlight=res2$med.number)
}

\keyword{htest}
