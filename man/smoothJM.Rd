\name{smoothJM}
\alias{smoothJM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Joint Model for the Analysis of Survival Data and Biomarkers
}
\description{
This function fits a joint model based on the heterogeneity mixed model described in Verbeke and Lesaffre (1996). The means of the mixture distribution specify the distribution of the random effects in the longitudinal part of the model and are also used as auxiliary variables in the survival part of the model. These means lie on a rectangular lattice, where the horizontal dimension represents the intercept, and the second dimension represents the slope. The means of the mixture distribution for the intercept and slope are included in the survival submodel to study the association between the biomarker trajectory and the time to event.
}
\usage{
smoothJM(formula, data, survival, maxiter = 200, rightskewed = c(T,T), n1 = 3, n2 = 3,
          eps = 0.0001, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a two-sided linear formula object describing both the fixed-effects and random-effects part of the model, with the response on the left of a ~ operator and the terms, separated by + operators, on the right. Random-effects terms are distinguished by vertical bars (|) separating expressions for design matrices from the grouping factor.
}
  \item{data}{
a data frame containing the variables named in the formula for the mixed model and the formula in the survival model. By default the variables are taken from the environment from which smoothLME is called.
}
\item{survival}{
a formula expression as for other survival models. The response is a survival object as returned by the Surv function, where the first argument is the time, and the second argument is the event. Only right censored survival times are supported in this model. See the documentation for Surv, survreg and formula for details.
}
  \item{maxiter}{
optional maximum number of iterations for the iterative algorithm. By default, maxiter=30.
}
  \item{n1}{
number of latent classes of the random effect for the intercept.
}
  \item{n2}{
number of latent classes for the variable, usually the time variable.
}
   \item{eps}{convergence criterion, the default value is 0.001}
   \item{verbose}{print the iteration history, the default is FALSE}
     \item{rightskewed}{ Boolean indicating whether the skewed random effects distribution is skewed to the left or skewed to the right. If \code{TRUE} a right-skewed distribution is fitted.  A right-skewed distribution has a long tail to the right of the mode, and the mean in this distributions is located to the right of the mode and the median.
  }

}
\details{
See the package vignette
}
\value{
The list returned is:
\item{jlcmm}{an object of class "jointlcmm"}
}
\references{
Verbeke G and Lesaffre E (1996). A linear mixed-effects model with heterogeneity in
the random-effects population. Journal of the American Statistical Association 91, 217-21
}
\author{
Berrie Zielman
}
\examples{
\dontrun{
library(lcmm)
data("pbc2", package = "JM")
pbc2$lsgot <- log(pbc2$SGOT)
model1 <- smoothJM(lsgot~ year + (1 + year | id), maxiter = 100,
      survival = Surv(years, status2) ~ factor(drug) + factor(sex), data = pbc2,
      n1 = 3, n2 = 3, rightskewed = c(T,T), verbose = TRUE)
plotrandomeffects(model1, which = "slope" , xlab = "slope")
plotrandomeffects(model1, which = "intercept", xlab = "intercept")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Survival Analysis}% use one of  RShowDoc("KEYWORDS")
\keyword{Mixed Model}% __ONLY ONE__ keyword per line
