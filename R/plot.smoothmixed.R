plot.smoothmixed <- function(x, ...){
  hlme <- x$hlme
  hlme$call$subject <- x$subj
  hlme$data <- x$dd2
  plot(hlme, ...)
}
