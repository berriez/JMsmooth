smoothmixed <- function(formula, data, maxiter = 300, n1 = 3, n2 = 3, verbose = FALSE, rightskewed = c(T, T), eps = .0001)
  {
  if(class(formula)!="formula") stop("The argument formula must be a formula")
  if(missing(data)){ stop("The argument data should be specified and defined as a data.frame")}
  if(nrow(data)==0) stop("Data should not be empty")
  if(maxiter<0) stop("specify at least a few iterations")
  if(n1<2) stop("specify at least 2 more components in this random effects distribution")
  if(n2<2) stop("specify at least 2 more components in this random effects distribution")
 # check formula
  f_0 <- as.formula(formula)
  f1 <- all.vars(f_0)
  data[order(paste("data",f1[length(f1)],sep='$')),]
  idvar <- f1[length(f1)]
  f2 <- as.formula(paste("~1",all.vars(lme4::findbars(formula(f_0))[[1]])[1],sep = "+"),env = parent.frame())
  paste(f2,f1[length(f1)],sep="|")
  ik <- grep(all.vars(lme4::findbars(formula(f_0))[[1]])[1],all.vars(formula(f_0)))
  sik <- all.vars(f_0)
  tmp <- sik[ik]
  sik[ik] <- sik[2]
  sik[2] <- tmp
 # fixed <- as.formula(paste(sik[1],"~",paste(sik[2:(nvar)],collapse = "+")))
  random <- paste("(",paste(f2,f1[length(f1)],sep="|")[2],")",sep="")
  formula <- paste(paste(sik[1],"~",paste(sik[2:(length(sik)-1)],collapse = "+")), random, sep="+")
  # end check formula
  nclas <- n1*n2
  right <- rightskewed
  a<-rep(1:n1,each=n2) #design for probabilities
  b<-rep(1:n2,times=n1)
  model   <- lme4::lmer(f_0, data = data,  REML = FALSE)
  z       <- getME(model,"Z")   #random design matrices
  x       <- getME(model,"X")   #fixed design matrices
  yp      <- getME(model,"y")   #response vector
  ind     <- as.numeric(unlist(getME(model,"flist")))
  nobs    <- length(table(ind))
  ndata   <- length(ind)
  nvar    <- ncol(x)
  nran    <- ncol(z)/nobs
  subjmat <- outer(ind, levels(factor(ind)), `==`)*1
  subjmat <- as(subjmat, 'sparseMatrix')
  alldata <- list(nclas)
  alldata <- lapply(1:nclas,function(xt,y) {alldata[[xt]] <- y},y = x)
  explain <- do.call(rbind,alldata)
  #ID <- names(getME(model,"l_i"))
  #stpoints<-c(0,cumsum(rle(ind)$length)) #for later use check
  postprob <- matrix(0,nrow = nobs, ncol = nclas)
  mu1 <- matrix(rep(0:(n1-1)),n1,1)
  mu2 <- matrix(rep(0:(n2-1)),n2,1)
  dd<-data.frame(a = gl(n1,ndata*n2)) #, b = gl(n2,ndata,n2*ndata)
  # dat <- data.frame(yp,x,ind)
  # names(dat) <-c(all.vars(formula(model) )[1],"intercept",all.vars(formula(model) )[2:nvar],ID) # was ind
  # ran2 <- diag(n2)%x%dat[,all.vars(lme4::findbars(formula(model))[[1]])[1]]
  ran2 <- diag(n2)%x%((z%*%(rep(1,nobs)%x%diag(nran)))[,2]) #is dit goed? dan kan wel alle data weg, dat enz
  alldata <- list(n1)
  alldata <- lapply(1:n1, function(xt,y) {alldata[[xt]]<- y}, y=ran2)
  desran2 <- do.call(rbind,alldata)
  desmat <- cbind(model.matrix(~a-1,dd)%*%mu1, desran2%*%mu2)
  colnames(desmat)<-c("bm_intercept","bm_slope")
  xtil <- cbind(desmat,explain)
  #
  # starting values should be adjusted p uniform, sort p
  #
  if(rightskewed[1])
    pow <- c(1:n1)
  else pow <- c(n1:1)
  th <- .32
  p1 <- th^pow/sum(th^pow)
  if(rightskewed[2])
    pow <- c(1:n2)
  else pow <- c(n2:1)
  th <- .42
  p2 <- th^pow/sum(th^pow)
  p <- p1%x%p2
  lcmmp <- log(p/p[nclas])[1:(nclas-1)]
  f1 <- as.formula(paste("~1",all.vars(lme4::findbars(formula(model))[[1]])[1], sep = "+"))
                   #,env = parent.frame())
  B2 <- rep(0.000000000009123,3*nclas+nvar+3-2)#aanpassen
  #
  # pas hier de scaling van de initiele waarden aan +fixef(model)
  #
  #
  l1 <- as.matrix(getME(model,"Lambda")[1:nran,1:nran]*sigma(model))
  # l12 <- l1%*%t(l1)
  # stdev <- sqrt(diag(l12))
  stdev <- diag(l1)
  b1 <- getME(model,"beta")
  B2[(nclas):(2*nclas-1)] <- rep(b1[1]+stdev[1]*mu1*.8,each=n2)#aanpassen
  B2[(nclas+nclas):(3*nclas-1)] <-  rep(b1[2]+stdev[2]*mu2*.8,times=n1)#aanpassen
  B2[1:(nclas-1)] <- lcmmp
  if(nvar-nran>0)
  {
    B2[(3*nclas):(3*nclas+nvar-nran-1)] <- b1[(nran+1):nvar]
  }
  B2[c((3*nclas+nvar-nran),(3*nclas+nvar-nran+2))] <- diag(l1)#maybe outside this if statement
  B2[(3*nclas+nvar-nran):(3*nclas+nvar-nran+2)] <- as.vector(l1[lower.tri(l1, diag = TRUE)])
  f0 <- as.formula(lme4::nobars(formula(model)))
  B2[length(B2)] <- sigma(model)
  m2 <- do.call("hlme",list(fixed = f0,
                            mixture = f1, random = f1,
                            subject = idvar, ng = nclas, data = data,
                            B = B2,
                            verbose = FALSE, maxiter = 120,
                            posfix = c(1:((3*nclas)-1))))
  B2<-m2$best
  loglikold <- m2$loglik
  B2[1:(nclas-1)] <- lcmmp
  postprob <- as.matrix(m2$pprob)[,3:(2+nclas)]
  nrow<-(sqrt(8*length(m2$cholesky)+1)-1)/2
  xvx <- matrix(0,(nvar+nran),(nvar+nran)) #change nclas + 1 to nclas + nvar or nvar + nran
  xvy <- matrix(0,(nvar+nran),1)
  c<-matrix(0,nrow,nrow)
  posfix <- c(1:((3*nclas)+nvar-3))
  nparms <- length(B2)
  #
  # start iterations
  #
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
    for(i in 1:nrow){
      for(j in i:nrow){
        c[i,j] <- m2$cholesky[(i-1)*i/2+j]
      }
    }
    se<-B2[nparms]
    d<-Diagonal(nobs)%x%(t(c)%*%c)               #kronecker product of identity times D (sparse matrix cla)
    VM <- solve(z%*%d%*%t(as.matrix(z))+diag(ndata)*se^2)   #ZVZ'+sI
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
    B2[(2*nclas):(3*nclas-1)] <-  rep(btil[4]+btil[2]*mu2,times=n1)
    B2[1:(nclas-1)] <- log(p/p[nclas])[1:(nclas-1)]
    m2 <- do.call("hlme",list(fixed = f0,
                              mixture = f1, random = f1,
                              subject = idvar, ng = nclas, data = data,
                              B = B2,
                              verbose = FALSE, maxiter = 120,
                              posfix = c(1:((3*nclas)-1))))

    B2<-m2$best
    if(verbose)
    print(m2$loglik)
    if((m2$loglik-loglikold)<eps) break
    loglikold <- m2$loglik
  }
  crosstable <- matrix(p, n1, n2, byrow = TRUE)
  beta.fit <- btil
  beta.sd <- sqrt(diag(solve(xvx)))
  beta.z <- beta.fit/beta.sd
  p.value <- 2 * (1 - pnorm(abs(beta.z)))
  coef <- data.frame(coefficient = beta.fit, SE = beta.sd, z = beta.z, p.value = p.value)
  rownames(coef) <- rownames(btil)
  bootstraplist <- list(f0 = f0, f1 = f1, nvar = nvar, ind = ind, formula = f_0, right = right, subj = idvar)
  res <-list(best = m2$best, n1 = n1, n2 = n2, crosstable = crosstable, data = data,  bootstraplist = bootstraplist,
             nrow = nrow, p1 = p1, p2 = p2,  nobs = nobs, ndata = ndata, btil = btil,
             ll_null = logLik(model), hlme = m2,  ll_model = m2$loglik, explain = explain,
             l1=l1, cholesky = m2$cholesky, xtil = xtil, nvar = nvar,z=z, subjmat = subjmat,
             nran=nran, yp = yp,
             niter = iter1, coef = coef)
  class(res) <-c("smoothmixed")
  return(res)
}

