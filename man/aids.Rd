\name{aids}
\alias{aids}
\docType{data}
\title{
Evolution of CD4 counts in aids patients
}
\description{
A randomized clinical trial in which both longitudinal and survival data were collected to compare the efficacy and safety of two antiretroviral drugs in treating patients who had failed or were intolerant of zidovudine (AZT) therapy.
}
\usage{data("aids")}
\format{
  A data frame with 1405 observations on the following 9 variables.
  \describe{
    \item{\code{patient}}{a factor with 467 levels}
    \item{\code{Time}}{a numeric vector}
    \item{\code{death}}{a numeric vector}
    \item{\code{CD4}}{a numeric vector}
    \item{\code{obstime}}{a numeric vector}
    \item{\code{drug}}{a factor with levels \code{ddC} \code{ddI}}
    \item{\code{gender}}{a factor with levels \code{female} \code{male}}
    \item{\code{prevOI}}{a factor with levels \code{noAIDS} \code{AIDS}}
    \item{\code{AZT}}{a factor with levels \code{intolerance} \code{failure}}
  }
}

\references{
Goldman, A., Carlin, B., Crane, L., Launer, C., Korvick, J., Deyton, L. and Abrams, D. (1996) Response of CD4+ and
clinical consequences to treatment using ddI or ddC in patients with advanced HIV infection. \emph{Journal of Acquired
Immune Deficiency Syndromes and Human Retrovirology} \bold{11}, 161--169.

Guo, X. and Carlin, B. (2004) Separate and joint modeling of longitudinal and event time data using standard
computer packages. \emph{The American Statistician} \bold{58}, 16--24.
}
\keyword{Datasets}
