################################## 
#   glasso related functions
##################################
library(glassoFast)
rglasso.path <- function(S, lambda, weight = TRUE, nlambda, ratio){
  
  if(weight){
    W <- sqrt(diag(S) %o% diag(S))
  }else{
    W <- matrix(1, nrow(S), ncol(S))
  }
  
  diag(W) <- 0
  
  if(missing(lambda)){
    lambda.max <- max(abs((S/W)[upper.tri(S)]))
    lambda.min <- ratio * lambda.max
    lambda <- exp(seq(log(lambda.min), log(lambda.max), length = nlambda))
  }
  
  f <- function(rho){
    res <- glassoFast::glassoFast(S, rho * W)
    return(list(w = res$w, wi = res$wi))
  }
  
  out <- lapply(lambda, FUN = f)
  w <- lapply(out, function(x) x$w)
  wi <- lapply(out, function(x) x$wi)
  
  return(list(w = w, wi = wi, lambda = lambda))
}

IC.rglasso <- function(out, s, n){
  
  f <- function(wi){
    nll <- -log(det(wi)) + sum(wi * s)
    J <- sum(wi[upper.tri(wi, diag = TRUE)] != 0)
    return(c(aic = nll + (2/n)*J, bic = nll + (log(n)/n)*J))
  }
  
  res <- sapply(out$wi, FUN = f)
  
  return(res)
}


# rglasso(covU,weight = TRUE, N=nrow(set1)-1)
# S <- covU
# weight = TRUE
# N=nrow(set1)-1
# nlambda = 10
# ratio <- sqrt(log(nrow(S))/N)

rglasso <- function(S, lambda = "aic", weight = TRUE, nlambda = 10, ratio, N){
  
  if(lambda == "aic" | lambda == "bic"){
    
    if(missing(ratio)) ratio <- sqrt(log(nrow(S))/N)
    
    out <- rglasso.path(S, weight = weight, nlambda = nlambda, ratio = ratio)
    ic <- IC.rglasso(out, S, N)
    
    ##For referees##
    bics <- ic[2, ]
    lambdas <- out[["lambda"]]
    precisions <- out$wi
    precisions2 <- out$wi
    for (a in 1: length(precisions)) {
      diag(precisions[[a]])<- 1
    }
    zeroes <- NULL
    for (i in 1:(length(precisions))) {
      zeroes = append(zeroes, as.numeric(length(precisions[[i]])-nnzero(precisions[[i]])) / as.numeric (length(precisions[[i]]) ),after = length(zeroes))
    }
    ################
    
    idx <- which.min(ic[lambda, ])
    #w <- out$w[[idx]]
    wi <- out$wi[[idx]]
    w <- solve(wi)
    lambda <- out$lambda[idx]
    
  }else{
    out <- rglasso.path(S, lambda = lambda, weight = weight)
    #w <- out$w[[1]]
    wi <- out$wi[[1]]
    w <- solve(wi)
  }
  
  attr(w, "lambda") <- lambda
  attr(wi, "lambda") <- lambda
  
  dimnames(w) <- dimnames(S)
  dimnames(wi) <- dimnames(S)
  return(list(cmat = w, pmat = wi, bics = bics, lambdas = lambdas, precisions = precisions, precisions_diag = precisions2, idx = idx,lambda = lambda, zeroes = zeroes))
  # return(list(cmat = w, pmat = wi, bics = bics, lambdas = lambdas, precisions = precisions, zeroes = zeroes))
}