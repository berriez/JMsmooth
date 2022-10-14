\docType{package}
\name{smoothJM-package}
\alias{smoothJM-package}
\title{Estimation of a joint model for the analysis of repeated measures and time to even}
\description{
We restrict to data in which the event times are either known exactly, or are right censored. To estimate this model, the data needed are the minimum of the event time and censoring time, denoted by \eqn{min(t_i, c_i).}

In addition we have a variable \eqn{\delta_i} indicating whether \eqn{t_i} was observed \eqn{\delta_i = 1} or not \eqn{\delta_i= 0}

The function for the class means in the linear mixed model is
\deqn{f_g(t) = \mu_{g1} + t \mu_{g2}}, where \eqn{\mu_{g1}} is the class mean of the intercept and \eqn{\mu_{g2}} is the class mean for the slope. The repeated measures are connected to the survival times by including the mean of the intercept and the rate of change  from the longitudinal model in the design matrix of the survival model.
}
\details{
provide some details
}
\references{
Proust-Lima C, Philipps V, Liquet B (2017). Estimation of Extended Mixed
Models Using Latent Classes and Latent Processes: The R Package lcmm.
Journal of Statistical Software, 78(2), 1-56. doi:10.18637/jss.v078.i02

Lin, Turnbull, McCulloch and Slate (2002). Latent class models for joint
analysis of longitudinal biomarker and event process data: application to
longitudinal prostate-specific antigen readings and prostate cancer. Journal
of the American Statistical Association 97, 53-65.

Muthen and Shedden (1999). Finite mixture modeling with mixture outcomes
using the EM algorithm. Biometrics 55, 463-9

Proust and Jacqmin-Gadda (2005). Estimation of linear mixed models with a
mixture of distribution for the random-effects. Comput Methods Programs
Biomed 78:165-73

Proust, Jacqmin-Gadda, Taylor, Ganiayre, and Commenges (2006). A nonlinear
model with latent process for cognitive evolution using multivariate
longitudinal data. Biometrics 62, 1014-24.

Proust-Lima, Dartigues and Jacqmin-Gadda (2011). Misuse of the linear mixed
model when evaluating risk factors of cognitive decline. Amer J Epidemiol
174(9), 1077-88

Proust-Lima and Taylor (2009). Development and validation of a dynamic
prognostic tool for prostate cancer recurrence using repeated measures of
post-treatment PSA: a joint modelling approach. Biostatistics 10, 535-49.

Proust-Lima, Sene, Taylor, Jacqmin-Gadda (2014). Joint latent class models
for longitudinal and time-to-event data: a review. Statistical Methods in
Medical Research 23, 74-90.

Proust-Lima, Amievan Jacqmin-Gadda (2013). Analysis of multivariate mixed
longitudinal data: A flexible latent process approach. Br J Math Stat
Psychol 66(3), 470-87.

Verbeke and Lesaffre (1996). A linear mixed-effects model with heterogeneity
in the random-effects population. Journal of the American Statistical
Association 91, 217-21
}
\keyword{package}
