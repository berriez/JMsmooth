myweib <-
  function(L, R, data, strata="ALL", weights=NULL, covariates=NULL, init = NULL) {
    straname <- "strata"
    ## Get strata, left time, right time
    strata <- eval(substitute(strata), data, parent.frame())
    Lt <- eval(substitute(L), data, parent.frame())
    Rt <- eval(substitute(R), data, parent.frame())
    if (length(strata) == 1) strata <- rep(strata, length(Lt))
    Lt <- as.vector(Lt)
    Rt <- as.vector(Rt)
    n <- nrow(data)
    if (!is.numeric(Lt) | !is.numeric(Rt)) stop("L and R must be numeric variables")
    strata <- as.vector(strata)

    ## Model Frame
    if (is.null(covariates)) {
      mframe <- NULL
      design <- NULL
      beta.nm <- NULL
      nbeta <- 0
    } else {
      mframe <-  model.frame(covariates, data = data, na.action = na.pass)
    }
    if (is.null(weights)) {
      w <- rep(nrow(data),1)
    } else {
      w <-  weights
    }

    ## Compare lengths
    if (length(unique(c(length(Lt), length(Rt), length(strata), dim(mframe)[1]))) != 1)
      stop("L, R, strata, and data must have same length")


    ## Design matrix
    if (!is.null(covariates)) {
      design <- model.matrix(covariates, data = mframe)[, -1, drop=F]
      beta.nm <- colnames(design)
      nbeta <- dim(design)[2]
    }

    ## Design matrix for strata
    levels <- sort(unique(strata))
    nstr <- length(levels)
    dstrata <- outer(strata, levels, "==")*1
    stra.nm <- paste(c(rep("v", nstr), rep("u", nstr)), rep(levels, 2), sep = ":")

    ## Time information
    allt <- c(Lt, Rt)
#    maxT <- max(allt[allt!=Inf])
    meanT <- mean(allt[allt!=0 & allt!=Inf])
    vini <- -log(meanT)

    ## Get event status to see if any event
    et <- Lt==Rt
    event <- sum(et) > 0
    if (event) {
      ie <- which(et)
      dstratae <- dstrata[et, ]
      designe <- design[et, ]
      Rt[et] <- Inf
    }

    ## Log-likelihood function
    LLt <- log(Lt)
    RRt <- log(Rt)
    loglik <- function(parms) {
      v <- dstrata%*%parms[1:nstr]
      ev <- exp(v)
      u <- dstrata%*%parms[(nstr+1):(2*nstr)]
      if (nbeta > 0) {prog <- design%*%parms[-(1:(2*nstr))]} else {prog <- rep(0, n)}
      lik <- sum(w*log(exp(-exp(u + prog + ev*LLt)) - exp(-exp(u + prog + ev*RRt))))
      if (event) lik <- lik + sum(w[ie]*(u[ie] + prog[ie] + v[ie] + (ev[ie]-1)*LLt[ie]))
      return(-lik)
    }

    ## Gradient function
    gradlik <-function(parms) {
      v <- dstrata%*%parms[1:nstr]
      ev <- c(exp(v))
      u <- dstrata%*%parms[(nstr+1):(2*nstr)]
      if (nbeta > 0) {prog <- design%*%parms[-(1:(2*nstr))]} else {prog <- 0}
      SL <- exp(-exp(u + prog + ev*LLt))
      SR <- exp(-exp(u + prog + ev*RRt))
      SLL <- c(SL*exp(u + prog + ev*LLt))
      SRR <- c(ifelse(Rt==Inf, 0, SR*exp(u + prog + ev*RRt)))
      Dev <- colSums(w*(cbind(ifelse(Lt==0, 0, LLt)*ev*dstrata, dstrata, design)*SLL -
                          cbind(ifelse(Rt==Inf, 0, RRt)*ev*dstrata, dstrata, design)*SRR)/c(SL-SR))
      if (event) Dev <- Dev - colSums(w[ie]*cbind((1+LLt[ie]*ev[ie])*dstratae, dstratae, designe))
      return(Dev)
    }

    ## Optimization
    if(is.null(init))
    {
      parmi <- c(rep(0, nstr), rep(vini, nstr), rep(0, nbeta))
    }
    else {parmi <- init}
    q <- optim(parmi, loglik, gradlik, method="BFGS", hessian=T, control = list(maxit=55))
    if (q$convergence!=0)	warning("Full model not converged")

    ## Covariance matrix
    covm <- solve(q$hessian)



    ## Combine test statistics

    ## report results
    rownames(covm) <- colnames(covm) <- c(stra.nm, beta.nm)
    if (nbeta > 0) {
      beta.fit <- q$par[-(1:(2*nstr))]
      beta.sd <- sqrt(diag(covm)[-(1:(2*nstr))])
      beta.z <- beta.fit/beta.sd
      p.value <- 2 * (1 - pnorm(abs(beta.z)))
      coef <- data.frame(coefficient = beta.fit, SE = beta.sd, z = beta.z, p.value = p.value)
      rownames(coef) <- beta.nm
    } else {
      coef <- NA
    }
    weib <- data.frame(straname=straname, strata=levels, gamma = exp(q$par[1:nstr]),
                       lambda = exp(q$par[(nstr+1):(2*nstr)]), stringsAsFactors=F)
    z <- list( coef=coef, weib=weib,estimate = q$par)
    class(z) <- "icweib"
    return(z)
  }
