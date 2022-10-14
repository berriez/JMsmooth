smoothmixedSE <- function(x, nboot = 100, maxiter = 100, eps = .001)
{
  data  <- x$data
  ind   <- x$bootstraplist$ind
  nvar  <- x$bootstraplist$nvar
  n1    <- x$n1
  n2    <- x$n2
  right <- x$bootstraplist$right
  ID <- x$bootstraplist$subj
  B2 <- x$best
  cr <- x$crosstab
  p  <- as.vector(cr)
  f0 <- x$bootstraplist$f0
  f1 <- x$bootstraplist$f1
  maxiter <- maxiter
  eps <- eps
  nparms <- length(B2)
  params <- c(x$p1$th, x$p2$th, x$btil, B2[(nparms-3) : nparms]) # check for one component
  names <- c("theta1","theta2", row.names(x$btil), names(B2[(nparms-3) : nparms]))
  a<-rep(1:n1,each=n2) #use for hazards
  b<-rep(1:n2,times=n1)
  nclas <- n1*n2
#  survival <- x$survival
  formula <- as.formula(x$bootstraplist$formula)
  nobs  <- length(table(ind))
  obslabels <- names(table(ind))
  results <- list(nboot)

  for (samp in 1:nboot){
  nums <- runif(nobs)
  alldata = list(nobs)
  alldata = lapply(1:nobs,function(xt) {alldata[[xt]] <-
    cbind(data[data[,ID]==obslabels[1+floor(nums[xt]*nobs)],],xt)})
  bootstrapdata = do.call(rbind,alldata)
  model   <- lme4::lmer(formula, data=bootstrapdata,  REML = FALSE)
  z       <- getME(model,"Z")   #random design matrices
  x       <- getME(model,"X")   #formula design matrices
  yp      <- getME(model,"y")   #response vector
  ind_b     <- as.numeric(unlist(getME(model,"flist")))
  nobs_b    <- length(table(ind_b))
  ndata   <- length(ind_b)
  nvar    <- ncol(x)
  nran    <- ncol(z)/nobs_b
  subjmat <- outer(ind_b, levels(factor(ind_b)), `==`)*1
  subjmat <- as(subjmat, 'sparseMatrix')
  # survnames <- all.vars(survival)
  # nsurvvar <- length(survnames)
  alldata = list(nclas)
  alldata = lapply(1:nclas,function(xt,y) {alldata[[xt]] <- y},y = x)
  explain = do.call(rbind,alldata)
  ID <- names(getME(model,"l_i"))
  mu1 <- matrix(rep(0:(n1-1)),n1,1)
  mu2 <- matrix(rep(0:(n2-1)),n2,1)
  dd<-data.frame(a = gl(n1,ndata*n2)) #, b = gl(n2,ndata,n2*ndata)
  dat <- data.frame(yp,x,ind_b)
  names(dat) <-c(all.vars(formula(model) )[1],"intercept",all.vars(formula(model) )[2:nvar],ID) # was ind
  ran2 <- diag(n2)%x%dat[,all.vars(lme4::findbars(formula(model))[[1]])[1]]
  alldata = list(n2)
  alldata = lapply(1:n1,function(xt,y) {alldata[[xt]]<- y},y=ran2)
  desran2 = do.call(rbind,alldata)
  desmat <- cbind(model.matrix(~a-1,dd)%*%mu1, desran2%*%mu2)
  colnames(desmat)<-c("bm_intercept","bm_slope")
  xtil <- cbind(desmat,explain)
  m2 <- do.call("hlme",list(fixed = f0,
                            mixture = f1, random = f1,
                            subject = ID, ng = nclas, data = bootstrapdata,
                            B = B2,
                            verbose = FALSE, maxiter = 120,
                            posfix = c(1:((3*nclas)-1))))
  nrow<-(sqrt(8*length(m2$cholesky)+1)-1)/2
  xvx <- matrix(0,(nvar+nran),(nvar+nran))
  xvy <- matrix(0,(nvar+nran),1)
  c<-matrix(0,nrow,nrow)
  B2<-m2$best
  loglikold <- m2$loglik
  #
  # start iterations
  #
  posfix <- c(1:((3*nclas)-1))
  posfix <- c(1:((3*nclas)+nvar-3))
  for (iter1 in 1:maxiter)
  {
  postprob <- as.matrix(m2$pprob)[,3:(2+nclas)]
  pold <- colMeans(postprob)
  p1<-tapply(pold,a,sum)
  p2<-tapply(pold,b,sum)
  p1 <- pestim(p1,n1,right[1])
  p2 <- pestim(p2,n2,right[2])
  p <- p1$prob%x%p2$prob
  lcmmp <- log(p/p[nclas])[1:(nclas-1)]
  B2[1:(nclas-1)] <- lcmmp

  #{
    for(i in 1:nrow){
      for(j in i:nrow){
        c[i,j]=m2$cholesky[(i-1)*i/2+j]
      }
    }
    se<-B2[length(B2)]
    d<-Diagonal(nobs_b)%x%(t(c)%*%c)     #kronecker product of identity times D (sparse matrix cla)
    VM <- solve(z%*%d%*%t(as.matrix(z))+Diagonal(ndata)*se^2)   #ZVZ'+sI
    weights <- as.vector(subjmat%*%postprob)
    xvx <- 0*xvx
    xvy <- 0*xvy
    xv1 <- sweep(xtil, 1, weights, FUN="*")
    for(i in 1:nclas){
      idx_end <- i*ndata
      idx_start <- idx_end-ndata+1
      temp <- VM%*%xv1[idx_start:idx_end,]
      xvx <- xvx + t(xtil[idx_start:idx_end,])%*%temp
      xvy <- xvy + t(temp)%*%yp
    }
    btil<-solve(xvx,xvy)
    B2[(nclas):(nclas+nclas-1)] <- rep(btil[3]+btil[1]*mu1,each=n2)
    B2[(nclas+nclas):(nclas+nclas+nclas-1)] <-  rep(btil[4]+btil[2]*mu2,times=n1)
 #   B2[1:(nclas-1)] <- log(p/p[nclas])[1:(nclas-1)]
    if(nvar>2)
    {
      B2[(3*nclas):(3*nclas+nvar-3)] <- btil[5:(2+nvar)]
    }
    m2 <- do.call("hlme",list(fixed = f0,
                              mixture = f1, random = f1,
                              subject = ID, ng = nclas, data = bootstrapdata,
                              B = B2,
                              verbose = FALSE, maxiter = 150, posfix = c(1:((3*nclas)-1))))
    B2<-m2$best
    print(m2$loglik)
    if((m2$loglik-loglikold)<0) print("divergence")
    if((m2$loglik-loglikold)<eps && iter1 > 3) break
    loglikold <- m2$loglik
   # }
  }
  print(samp)
  results[[samp]] <- c(p1$th,p2$th,btil,B2[(nparms-3):nparms])
} # end bootstrap iterations
    #rbind(results)
   mat <- do.call(cbind,results)
   mat <- t(mat)
 # colnames(mat) <- names
  mat <- as.data.frame(mat)
  mat <- do.call(cbind,results)
  mat <- t(mat)
  mean2 <- apply(mat,2,mean)
  var2 <- apply(mat,2,var)
  std <- sqrt(var2)
  lowv <- apply(mat,2,minval,nb=nboot)
  maxv <- apply(mat,2,maxval,nb=nboot)
  mat2 <- cbind(params,mean2,std,lowv,maxv)
  row.names(mat2) <- names
  return (mat2)
}
minval <- function(dd, nb)
{
  sort(as.numeric(dd))[max(0.025 * nb,1)]
}
maxval <- function(dd, nb)
{
  sort(as.numeric(dd))[0.975 * nb]
}
