\name{smoothmixedSE}
\alias{smoothmixedSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates Standard Errors using the Bootstrap
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function takes a model fit from a smoothmixed model, and uses the bootstrap to       calculate standard errors, with optional confidence intervals, for the
 coefficients in the model. To each bootstrap dataset the model is fitted and the results are saved in a new dataset.
}
\usage{
smoothmixedSE(x, nboot = 100, maxiter = 100, eps = .001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class smoothmixed
}
  \item{nboot}{
argument specifying the number of bootstrap samples. Note that at least nboot = 100 is required in order for the function to return non-zero confidence intervals.
}
  \item{maxiter}{the maximum number of iterations for each bootstrap sample}
  \item{eps}{the convergence criterion}

}
\details{
Standard errors and confidence intervals are obtained by repeated fitting of the  model to bootstrap samples of the original longitudinal and survival data. As a rule of thumb, 100 bootstrap samples are usually sufficient for estimating a standard error, whereas the number of bootstrap samples needed for confidence intervals can be as large as 1000.
}
\value{
A matrix where the rows refer to the parameters and  the columns estimator, the mean, standard deviation, lower limit of the confidence interval, and upper limit of the confidence interval.
}
\keyword{Mixed Model }% __ONLY ONE__ keyword per line
