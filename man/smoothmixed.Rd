\name{smoothmixed}
\alias{smoothmixed}
\alias{heterogeneity model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixed Model with a Smooth Random Effects Distribution
}
\description{
This function fits a random effects model with a smooth random effects distribution to data. Usually, the random effects are assumed to be normally distributed in a mixed model, the function in this packages replaces the normal distribution by a constrained mixture distribution.  A mixture of
normal distributions provides a flexible tool for modelling
non-normal random effects distributions, such as a skewed or
bimodal distribution. This function fits mixed models with right or left skewed distributions for the random effects. This leads to an asymmetric distribution around the mean. For instance, a right skewed distribution indicates that the tail on the right side is longer than the left side. The precise form of the random effects distribution is described in the package vignette.
}
\usage{
smoothmixed(formula, data, maxiter = 300, n1 = 3, n2 = 3, verbose = FALSE,
rightskewed = c(T,T), eps = .0001)
}
%- maybe also 'usage' for other objects documented here. fixed maybe mixed
\arguments{
  \item{formula}{
a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars (|) separating expressions for design matrices from grouping factors.
}
  \item{data}{
a data frame containing the variables named in formula. By default the variables are taken from the environment from which smoothmixed is called.
}
  \item{maxiter}{
optional maximum number of iterations for the iterative algorithm. By default, maxiter=300.
}
  \item{n1}{
number of latent classes of the random effect for the intercept.
}
  \item{n2}{
number of latent classes for the slope, this is usually the time variable.
}
  \item{verbose}{
if TRUE, prints the history of iterations
}
  \item{eps}{
convergence criterion
}

  \item{rightskewed}{ A list of booleans indicating the shape of the skewed random effects distribution.  A right-skewed distribution has a long right tail. The mean of a right skewed distribution is to the right of the mode and the median, whereas the opposite is true for a left skewed distribution.
  }

}
\value{
The list returned is:
\item{coef}{the coefficients of the fixed effects in the mixture model. The first two rows of this table refer to the distances between the class means of the intercept and the slope}
\item{crosstable}{table of mixture probabilities. The rows of this table refer to the the class means of the intercept, and the columns refer to the class means of the slopes}
\item{hlme}{an object of class hlme}
\item{niter}{number of iterations of the fitting algorithm}
\item{p1}{an object containing the mixture probabilities and theta parameter that generates these              probabilities. p1 contains these parameters for the distribution of the intercept}
\item{p2}{an object containing the mixture probabilities and theta parameter that generates these              probabilities. p2 contains these parameters for the distribution of the slope}
\item{ll_null}{loglikelihood of the model with a normal distribution for the random effects}
\item{ll_model}{loglikelihood of the model with a smooth distribution for the random effects}

}
\examples{
\dontrun{
data("pbc2", package = "JM")
pbc2$lsgot <- log(pbc2$SGOT)
lme <- smoothmixed(lsgot ~ year + factor(sex) + (1 + year | id), maxiter = 150,
                   data = pbc2, n1 = 3, n2 = 3)
lme$crosstable
plot(lme$hlme)
}
}
\details{
This function is an extension of the hlme function in the R package lcmm. The plot method and also other accessor methods such as coef also work for hlme component in this object. For instance,
the line \code{plot(lme$hlme)}
produces plots of residuals and fitted values using the plot method of the package lcmm.
The iterative maximization procedure cycles thru blocks op parameters. The parameters in a block are updated, while keep the parameter values in the other blocks fixed at their current values. We distinguish the following blocks of parameters: the variance components, the probabilities and the fixed effects in the model. In the first step we estimate the fixed effects using a weighted least squares estimator while keeping the other parameters fixed at their current values. In the second step we update the variance components, that is the residual error and the covariance of the random effects. The parameters for the variance components are extracted from a constrained fit of the lcmm model.
From the fitted lcmm model the posterior probabilities of the components are obtained. To these probabilities an independence model for the two random effects is fitted. Using a line search the power law constraints are imposed on these two sets of probabilities. These steps are an implementation of a constrained EM algorithm, where first the expected loglikelihood is calculated, and in the second step the parameters are estimated. This algorithm is a variant of a block relaxation algorithm. The advantage of this approach is that we can use the lcmm package as a building block in the algorithm.
}
\references{

Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed
Models Using Latent Classes and Latent Processes: the R package lcmm.
Journal of Statistical Software, 78(2), 1-56

Verbeke G and Lesaffre E (1996). A linear mixed-effects model with heterogeneity in
the random-effects population. Journal of the American Statistical Association 91, 217-21

}
\keyword{Mixed Model }% __ONLY ONE__ keyword per line

