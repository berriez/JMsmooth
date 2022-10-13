pestim <- function(dat, nclas, right = TRUE)
{
  if(right)
     pow <- c(1 : nclas)
  else pow <- c(nclas : 1)
     xmin <- optimize(function(x,a) sum(a*log(x^pow/sum(x^pow))), c(0, 1), maximum = TRUE,
                      tol = 0.0001, a = dat)
     prob <- xmin$maximum^pow/sum(xmin$maximum^pow)
     theta <- xmin$maximum
     result <- list(prob = prob, theta = theta)
  return(result)
}
