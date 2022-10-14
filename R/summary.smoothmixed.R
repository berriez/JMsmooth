summary.smoothmixed <- function(object,...){
    x <- object
    if (!inherits(x, "smoothmixed")) stop("use only with \"smoothmixed\" objects")
    cat("Linear mixed model with a smooth random effects distribution", "\n")
    cat("     fitted by maximum likelihood method", "\n")
    cat("Summary of iterations:", "\n")
    cat(paste("     Number of subjects:", x$nobs),"\n")
    cat(paste("     Number of observations:", x$ndata),"\n")
    cat(paste("     Number of iterations:", x$niter,"\n"))
    print(x$coef)
    return(invisible(x$coef))
}


