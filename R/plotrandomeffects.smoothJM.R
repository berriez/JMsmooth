plotrandomeffects <- function(x, which = "intercept", overlaid = TRUE, ...)
  UseMethod("plotrandomeffects")
plotrandomeffects.smoothJM <- function(x, which = "intercept", overlaid = TRUE,...){
  if(!which%in%c("intercept","slope"))
    stop ( "the argument which should be equal to intercept, slope, or joint" )
  nclas <- x$nclas
  nsvar <- x$nscoef
  n1 <- x$n1
  n2 <- x$n2
  if(which=="intercept") {
    var <- x$best[length(x$best)-3] #variance
    var <- as.numeric(var)
    s_dist <- sqrt(var)
    mus <- x$best[(3*nclas+nsvar) : (4*nclas+nsvar-1)] #n vars in survival model
    mus <- unique(mus)
    p <- rowSums(x$crosstable)
    mu_appr <- sum(p*mus)
    sigma <- 0
    for (i in 1:length(p))
    {

      sigma <- sigma + p[i]*(mus[i]-mu_appr)^2
    }
    sigma <- sigma + var

    xas <- seq(min(mus)-2.5*s_dist,max(mus)+2.5*s_dist,.001)
    plot(xas,mixnorm(p,xas,mus,s_dist),type="l", ...)
    if(overlaid)
    lines(xas,dnorm(xas,mu_appr,sqrt(sigma)),type="l",col="blue")

  } else if(which=="slope"){
    var <- x$best[length(x$best)-1] #variance
    s_dist <- sqrt(var)
    mus <- x$best[(4*nclas+nsvar):(4*nclas+n2+nsvar-1)] #n vars in survival model
    mus <- unique(mus)
    p <- colSums(x$crosstable)
    mu_appr <- sum(p*mus)
    sigma <- 0
    for (i in 1:length(p))
    {

      sigma <- sigma + p[i]*(mus[i]-mu_appr)^2
    }
    sigma <- sigma + var

    xas <- seq(min(mus)-2.5*s_dist,max(mus)+2.5*s_dist,.001)
    plot(xas,mixnorm(p,xas,mus,s_dist),type="l", ...)
    if(overlaid)
    lines(xas,dnorm(xas,mu_appr,sqrt(sigma)),type="l",col="blue")
  }
}

mixnorm <- function(p,x,mu,sd)
{
  val <- 0
  for (i in 1:length(p))
  {
    val <- val+p[i]*dnorm(x,mu[i],sd)
  }
  val
}


