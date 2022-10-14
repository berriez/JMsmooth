smoothJM <- function(formula, data, survival, maxiter = 200, rightskewed = c(T,T), n1 = 3, n2 = 3, eps = 0.0001, verbose = FALSE)
{
  nclas <- n1*n2
  a<-rep(1:n1,each=n2) # use for hazards
  b<-rep(1:n2,times=n1)
  survnames <- all.vars(survival)
  nsurvvar <- length(survnames)
  survtypes <- rep(FALSE,nsurvvar)
  for (i in 1:nsurvvar)
    survtypes[i] <- (is.factor(data[,survnames[i]]))
  mu1 <- matrix(rep(0:(n1-1)), n1, 1) #nakijken (of de oude nemen)
  mu2 <- matrix(rep(0:(n2-1)), n2, 1) #nakijken
  form_vars <- all.vars(formula)
  survnames <- all.vars(survival)

  allv <- c(all.vars(survival),all.vars(formula))
  data <- data[complete.cases(data[allv]),allv]

  m2 <- smoothmixed(formula, data =  data, n1 = n1, n2 = n2, maxiter = 1)
  f0 <- m2$bootstraplist$f0
  f1 <- m2$bootstraplist$f1

  nobs <- m2$nobs
  nvar <- m2$nvar
  nran <- m2$nran
  ind <- m2$ind
  z <- m2$z
  tz <- t(as.matrix(z))
  postprob <- as.matrix(m2$hlme$pprob)[,3:(2+nclas)]
  stdev <- diag(m2$l1)
  ndata   <- m2$ndata
  explain <- m2$explain
  subjmat <- m2$subjmat
  yp <- m2$yp
  dd<-data.frame(a = gl(n1,ndata*n2), b = gl(n2,ndata,n2*ndata)) #nakijken in survival
  xtil <- m2$xtil #from smoothmixed
  btil <- m2$btil
  if(rightskewed[1])
    pow <- c(1:n1)
  else pow <- c(n1:1)
  th <- .52
  p1 <- th^pow/sum(th^pow)
  if(rightskewed[2])
    pow <- c(1:n2)
  else pow <- c(n2:1)
  th <- .52
  p2 <- th^pow/sum(th^pow)
  p <- p1%x%p2


  survdata  <- unique(data[,cbind(survnames,form_vars[length(form_vars)])]) #id toevoegen same order as dataframe
  #survda = list(nsurvvar)
  #survda = lapply(1:nsurvvar,function(xt) {survda[[xt]]<- rep(survdata[,xt],drop(table(ind)))})
  #s2 <- do.call(cbind,survda) #nodig??
  nsurvdata <- nrow(survdata)
  ID <- form_vars[length(form_vars)]
  dd<-data.frame(a = gl(n1,nsurvdata*n2), b = gl(n2,nsurvdata,n2*nsurvdata))
  desmat <- cbind(model.matrix(~a -1,dd)%*%mu1, model.matrix(~b -1,dd)%*%mu2)
  alldata = list(nclas)
  alldata = lapply(1:nclas,function(xt,y) {alldata[[xt]]<- y},y=survdata)
  survdata = do.call(rbind,alldata)
  left <- survdata[,survnames[1]]
  #  time <- survdata[,survnames[1]]
  right <- ifelse(survdata[,survnames[2]]==1, survdata[,survnames[1]], Inf)
  status <- survdata[,survnames[2]]
  survdata <- cbind(survdata,left,right)
  colnames(desmat)<-c("bm_intercept","bm_slope")
  w <- as.vector(postprob)
  sd<-cbind(desmat,survdata,w,status)
  fslcmm <- as.formula(survival)
  fs <- update(fslcmm, ~ . +bm_intercept+bm_slope) #nodig in survival, maar niet in lcmm
  #
  # start of survival model
  #
  smodel <- survreg(fs, dist = "weibull", data=sd, weights = w+.001)
  # convert parameters
  gam <- 1/smodel$scale
  scoef <- -1*gam*smodel$coef #omrekenen naar ph
  params <- c(log(gam),scoef)
  nscoef <- nsurvvar-2
  vlist<-all.vars(fs)
  #
  # kan weg?
  #
  fs2 <- as.formula(paste("~",paste(vlist[-(1:2)], collapse=" + ")))
  fit1 <- myweib(L = left, R = right, init = params,  data = sd, weights=w,
                 strata = "ALL", covariates = fs2)
  gam <- log(fit1$weib[3])
  params <- fit1$estimate
  svals <- matrix(0,2*nclas,1)
  B3 <- m2$best
  # start the after party
  B2 <- rep(0.0123,5*nclas+nvar+nscoef+1)
  # a<-rep(mu1*scoef[3],each=n2) #use for hazards adjust index
  # b<-rep(mu2*scoef[4],times=n1)
  # baseline <- a+b+scoef[1]
  baseline <- rep(mu1*scoef[nscoef+2],each=n2)+rep(mu2*scoef[nscoef+3],times=n1)+scoef[1]
  for (i in 1:nclas){
    svals[1+(i-1)*2] <-  baseline[i]
    svals[2+(i-1)*2] <-  log(gam)
  }
  B2[1:(nclas-1)] <- log(p/p[nclas])[1:(nclas-1)]

  B2[nclas:(3*nclas-1)] <- unlist(svals)
  if(nscoef>0)
    B2[(3*nclas):(3*nclas+nscoef-1)]<- params[(3:(2+nscoef))]
  B2[(3*nclas+nscoef):(4*nclas+nscoef-1)] <- rep(btil[3]+btil[1]*mu1,each=n1)
  #
  #
  B2[(length(B2)-4):length(B2)] <- B3[(length(B3)-4):length(B3)] # -6 aanpassen bij random intercept only
  B2[(4*nclas+nscoef):(5*nclas+nscoef-1)] <-  rep(btil[4]+btil[2]*mu2,times=n2)
  fixpos <- c(1:((5*nclas+nscoef)-1))
  lcjmFit <- do.call("Jointlcmm",list(fixed = f0,
                                      mixture = f1, random = f1,
                                      subject = ID, ng = nclas, data = data, #was dat
                                      survival = as.formula(survival), B = B2,
                                      hazard = "Weibull", verbose = FALSE, maxiter = 120,
                                      logscale = TRUE, posfix = fixpos))
  fit1 <- myweib(L = left, R = right, init = params,  data = sd, weights=w,
                 strata = "ALL", covariates = fs2)
  #
  #
  #
  nrow<-(sqrt(8*length(m2$cholesky)+1)-1)/2
  c<-matrix(0,nrow,nrow)
  xvx <- matrix(0,(nvar+nran),(nvar+nran)) #change nclas + 1 to nclas + nvar or nvar + nran
  xvy <- matrix(0,(nvar+nran),1)
  loglik <- lcjmFit$loglik
  B2 <- lcjmFit$best  #aanpassen
  nparm <- length(B2)
  survda <- cbind(desmat,survdata,status)
  for (iter in 1:maxiter)
  {
    postprob <- as.matrix(lcjmFit$pprob)[,3:(2+nclas)]
    pold <- colMeans(postprob)
    p1<-tapply(pold,a,sum)
    p2<-tapply(pold,b,sum)
    p1 <- pestim(p1,n1,rightskewed[1])
    p2 <- pestim(p2,n2,rightskewed[2])
    p <- p1$prob%x%p2$prob
    B2<-lcjmFit$best
    B2[1:(nclas-1)] <- log(p/p[nclas])[1:(nclas-1)]
    w <- as.vector(postprob)
    for(i in 1:nrow){
      for(j in i:nrow){
        c[i,j]=lcjmFit$cholesky[(i-1)*i/2+j]
      }
    }
    se<-lcjmFit$best[nparm]                 #residual std
    VM <- solve(z%*%(Diagonal(nobs)%x%(t(c)%*%c))%*%tz+diag(ndata)*se^2)   #ZVZ'+sI
    weights <- as.vector(subjmat%*%postprob)
    xvx <- 0*xvx #matrix(0,(nvar+nran),(nvar+nran)) #change nclas + 1 to nclas + nvar or nvar + nran
    xvy <- 0*xvy #matrix(0,(nvar+nran),1)
    xv1 <- sweep(xtil, 1, weights, FUN="*")
    for(i in 1:nclas){
      idx_end <- i*ndata
      idx_start <- idx_end-ndata+1
      temp <- VM%*%xv1[idx_start:idx_end,]
      xvx <- xvx + t(xtil[idx_start:idx_end,])%*%temp
      xvy <- xvy + t(temp)%*%yp
    }
    btil<-solve(xvx,xvy)
    #
    # include myweib
    #
    #survda2 <- matrix(btil[c(1,2)],nrow(sd),2,byrow=TRUE)*survda[,c(1,2)]
    #survda2 <- cbind(survda2,survda[,-c(1,2)])
    sd<-cbind(survda,w)
    fit1 <- myweib(L = left, R = right, init = params,  data = sd, weights=w,
                   strata = "ALL", covariates = fs2)
    gam <- log(fit1$weib[3])
    params <- fit1$estimate
    baseline <-  rep(mu1*fit1$coef[nscoef+1,1], each=n2)   + rep(mu2*fit1$coef[nscoef+2,1], times=n1)
    #baseline <-  rep(mu1*fit1$coef[nscoef+1,1], each=n2) *btil[1]  + rep(mu2*fit1$coef[nscoef+2,1], times=n1)*btil[2]

    for (i in 1:nclas)
    {
      svals[1+(i-1)*2] <-  baseline[i] + log(fit1$weib[4])
      svals[2+(i-1)*2] <-  gam
    }
    B2[nclas:(3*nclas-1)] <- unlist(svals)
    if(nscoef>0)
      B2[(3*nclas):(3*nclas+nscoef-1)]<- params[(3:(2+nscoef))] #do in mixed also
    B2[(3*nclas+nscoef):(4*nclas+nscoef-1)] <-  rep(btil[3]+btil[1]*mu1,each = n2)
    B2[(4*nclas+nscoef):(5*nclas+nscoef-1)] <-  rep(btil[4]+btil[2]*mu2,times = n1)
    lcjmFit <- do.call("Jointlcmm",list(fixed = f0,
                                        mixture = f1, random = f1,
                                        subject = ID, ng = nclas, data = data,
                                        survival = as.formula(survival), B = B2, verbose = FALSE,
                                        hazard = "Weibull", maxiter = 120,
                                        logscale = TRUE, posfix = fixpos))
    loglikold <- loglik
    loglik <- lcjmFit$loglik
    if((loglik-loglikold)<eps) break
    if(verbose) print(loglik)
  }
  crosstable <- matrix(p, n1, n2, byrow = TRUE)
  res <-list(best = lcjmFit$best, n1 = n1, n2 = n2, crosstable = crosstable, data = data, f0 = f0,
             f1 = f1, pprob = lcjmFit$pprob, formula = formula, p1 = p1, p2 = p2,
             nclas = nclas,
             subjmat = subjmat, nscoef = nscoef, survival = survival, baseline = baseline,
             yp = yp, nrow = nrow, p = p, ind = ind, nobs = nobs, btil = btil, sd = sd, niter = iter,
             skew = rightskewed, jlcmm = lcjmFit, survmodel = fit1$coef, weibull = fit1$weib, loglik = loglik,
             survdata = survdata)
  class(res) <- c("smoothJM")
  return(res)
}

