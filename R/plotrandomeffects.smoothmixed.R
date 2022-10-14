plotrandomeffects.smoothmixed <- function(x, which = "intercept", overlaid = TRUE, ...){
  if (!which %in% c("intercept", "slope"))
    stop ("the argument which should be equal to intercept or slope")
  n1 <- x$n1
  n2 <- x$n2
  nclas <- n1 * n2
  if (which == "intercept") {
     var <- x$best[length(x$best) - 3] #variance
     var <- as.numeric(var)
     s_dist <- sqrt(var)
     mus <- x$best[(nclas):(2 * nclas - 1)]
     mus <- unique(mus)
     p <- rowSums(x$crosstable)
     mu_appr <- sum(p  * mus)
     sigma <- 0
    for (i in 1:length(p)){
      sigma <- sigma + p[i] * (mus[i] - mu_appr) ^ 2
    }
    sigma <- sigma + var
    xas <- seq(min(mus) - 2.5 * s_dist, max(mus) + 2.5 * s_dist, .001)
    plot(xas, mixnorm(p, xas, mus, s_dist), type = "l", ...)
    if (overlaid)
    lines(xas, dnorm(xas, mu_appr, sqrt(sigma)), type = "l", col = "blue")
    } else if (which == "slope"){
    var <- x$best[length(x$best) - 1] #variance
    s_dist <- sqrt(var)
    mus <- x$best[(2 * nclas):(3 * nclas - 1)]
    mus <- unique(mus)
    p <- colSums(x$crosstable)
    mu_appr <- sum(p * mus)
    sigma <- 0
    for (i in 1:length(p)) {
      sigma <- sigma + p[i] * (mus[i] - mu_appr) ^ 2
    }
    sigma <- sigma + var
    xas <- seq(min(mus) - 2.5 * s_dist, max(mus) + 2.5 * s_dist, .001)
    plot(xas, mixnorm(p, xas, mus, s_dist), type = "l", ...)
    if (overlaid)
    lines(xas, dnorm(xas, mu_appr, sqrt(sigma)), type = "l", col = "blue")
    }
}

mixnorm <- function(p, x, mu, sd){
  val <- 0
  for (i in 1:length(p)){
    val <- val + p[i] * dnorm(x, mu[i], sd)
  }
  val
}
plotrandomeffects <- function(x, which = "intercept", overlaid = TRUE, ...)
  UseMethod("plotrandomeffects")
