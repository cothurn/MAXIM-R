#' Binomial Cumulative Distribution Function
#'
#' @param n a non-negative integer.
#' @param p A value between 0 and 1, inclusive.
#' @return y(i) = P(X <= i-1), i=1,2,...,n+1, where X is a  Binomial(n,p) rv.
#' @examples
#' bincdf(n = 10, p = 0.5)
bincdf <- function(n=1,p=0.5){
  if(p < 0 || p >1){
    print("error with value p")
  }
  else if(n < 0){
    print("error with value n")
  }
  else{
    y <- c()
    for(k in 0:n)
      y = c(y,pbinom(k, size=n, prob=p))
  }
  return(y)
}

#' Binomial Probability Mass Function
#'
#' @param n a non-negative integer.
#' @param p A value between 0 and 1, inclusive.
#' @return y(i) = P(X = i-1),  i=1,...n+1, where X  is a  Binomial(n,p) rv.
#' @examples
#' bincdf(n = 10, p = 0.5)
binpmf <- function(n=1,p=0.5){
  if(p < 0 || p >1){
    print("error with value p")
  }
  else if(n < 0){
    print("error with value n")
  }
  else{
    y <- c()
    for(k in 0:n)
      y = c(y,dbinom(k, size=n, prob=p))
  }
  return(y)
}

#' Erlang Cumulative Distribution Function
#' @param k ingeter >= 1
#' @param l integer >= 0
#' @param x a row vector
#' @return y(i) =  F(x(i)), where F is an Erlang(k,l) cdf.
#' @examples
#' erlangcdf(1, 1, c(0.4,0.3,0.2,0.1))
erlangcdf <- function(k=1,l=1,x=c(0.4,0.3,0.2,0.1)){
  if(k < 1){
    print("error with value k")
  }
  else if(l < 0){
    print("error with value l")
  }
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,pgamma(x[j],k,l))
    }
    return(y)
  }
}

#' Plot of Erlang Cumulative Distribution Function
#' @param k ingeter >= 1
#' @param l integer >= 0
#' @param x a row vector
#' @return Plot of the cdf of Eralng(k,l) (y vs. x).
#' @examples
#' erlangcdfplot(1, 1, c(0.4,0.3,0.2,0.1))
erlangcdfplot <- function(k=1,l=1,x=c(0.4,0.3,0.2,0.1)){
  y = erlangcdf(k,l,x)
  plot(x, y)
  lines(x, y)
}

#' Eralng Probability Density Function
#' @param k ingeter >= 1
#' @param l >= 0
#' @param x a non negative row vector
#' @return y(i) =f(x(i)), where f is an Erlang(k,l) pdf.
#' @examples
#' erlangpdf(1, 1, c(0.4,0.3,0.2,0.1))
erlangpdf <- function(k=1,l=1,x=c(0.4,0.3,0.2,0.1)){
  if(k < 1){
    print("error with value k")
  }
  else if(l < 0){
    print("error with value l")
  }
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,dgamma(x[j],k,l))
    }
    return(y)
  }
}

#' Plot of Erlang Probability Density Function
#' @param k ingeter >= 1
#' @param l >= 0
#' @param x a non negative row vector
#' @return Plot of the pdf of Eralng(k,l) (y vs. x).
#' @examples
#' erlangpdfplot(1, 1, c(0.4,0.3,0.2,0.1))
erlangpdfplot <- function(k=1,l=1,x=c(0.4,0.3,0.2,0.1)){
  y = erlangpdf(k,l,x)
  plot(x, y)
  lines(x, y)
}

#' Exponential Cumulative Distribution Function
#' @param l >= 0
#' @param x a non negative row vector
#' @return y(i) =  F(x(i)), where F is an Exp(l) cdf
#' @examples
#' expcdf(1,c(0.4,0.3,0.2,0.1))
expcdf <- function(l=1,x=c(0.4,0.3,0.2,0.1)){
  if(l < 0){
    print("error with value l")
  }
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,pexp(x[j],l))
    }
    return(y)
  }
}

#' Plot of Exponential Cumulative Distribution Function
#' @param l >= 0
#' @param x a non negative row vector
#' @return Plot of the cdf of exp(l) (y vs. x).
#' @examples
#' expcdfplot(1,c(0.4,0.3,0.2,0.1))
expcdfplot <- function(l=1,x=c(0.4,0.3,0.2,0.1)){
  y = expcdf(l,x)
  plot(x,y)
  lines(x,y)
}

#' Title Exponential Probability Density Function
#' @param l >= 0
#' @param x a non negative row vector
#' @return y(i) =f(x(i)), where f is an Exp(l) pdf.
#' @examples
#' exppdf(1,c(0.4,0.3,0.2,0.1))
exppdf <- function(l=1,x=c(0.4,0.3,0.2,0.1)){
  if(l < 0){
    print("error with value l")
  }
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,dexp(x[j],l))
    }
    return(y)
  }
}


exppdfplot <- function(l=1,x=c(0.4,0.3,0.2,0.1)){
  y = exppdf(l,x)
  plot(x,y)
  lines(x,y)
}
geometriccdf <- function(p,k){
  if(!(k>=1))
    print("error with k")
  else if(p > 1 || p < 0)
    print("erro with p")
  else{
    y <- c()
    for(j in 1:k)
      y = c(y,pgeom(j,p))
    return(y)
  }
}
geometricpmf <- function(p,k){
  if(!(k>=1))
    print("error with k")
  else if(p > 1 || p < 0)
    print("error with p")
  else{
    y <- c()
    for(j in 0:k+1)
      #print(list(j-1,dgeom(j-1,p)))
      y = c(y,dgeom(j-1,p))
    return(y)
  }
}
negbincdf <- function(r,p,k){
  if(r < 1)
    print("Issue with r")
  else if(k < 1)
    print("issue with k")
  else{
    y <- c()
    for(j in r:k)
      y = c(y,pnbinom(j,r,p))
    return(y)
  }
}
negbinpmf <- function(r,p,k){
  if(r < 1)
    print("Issue with r")
  else if(k < 1)
    print("issue with k")
  else{
    y <- c()
    for(j in r:k)
      y = c(y,dnbinom(j,r,p))
    return(y)
  }
}
normalcdf <- function(m,s,x){
  if(s < 0)
    print("Issue with s")
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,pnorm(x[j],m,s))
    }
    return(y)
  }
}
normalcdfplot <- function(m,s,x){
  y = normalcdf(m,s,x)
  plot(x,y)
  lines(x,y)
}
normalpmf <- function(m,s,x){
  if(s < 0)
    print("Issue with s")
  else{
    y <- c()
    for(j in 1:length(x)){
      y = c(y,dnorm(x[j],m,s))
    }
    return(y)
  }
}
normalpmfplot <- function(m,s,x){
  y = normalpmf(m,s,x)
  plot(x,y)
  lines(x,y)
}
poissoncdf <- function(l,k){
  if(l < 0 || l > 700)
    print("Issue with l")
  else{
    y <- c()
    for(j in 0:k){
      y = c(y,ppois(j,l))
    }
    return(y)
  }
}
poissonpmf <- function(l,k){
  if(l < 0 || l > 700)
    print("Issue with l")
  else{
    y <- c()
    for(j in 0:k){
      y = c(y,dpois(j,l))
    }
    return(y)
  }
}
