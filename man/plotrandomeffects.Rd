\name{plotrandomeffects}
\alias{plotrandomeffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of the random effects distribution
}
\description{
Plots of the random effects distribution in the mixed effects model. The marginal distributions of the random intercepts or slopes can be plotted.
}
\usage{
plotrandomeffects(x, which = 'intercept', overlaid = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class smoothmixed
}
  \item{which}{indicates which plot to draw. By default a plot of the distribution of the intercept is produced. Other possibilities are "slope" for the distribution of the slope.}
  \item{overlaid}{The marginal distribution is shown with an overlaid normal distribution with the same mean and variance as the mixture distribution. This option is usefull for detecting deviations from the normal distribution of the fitted distribution.
  }
  \item{...}{further plot parameters}
}
\keyword{Model checking}% __ONLY ONE__ keyword per line
\keyword{Random effects distribution}
