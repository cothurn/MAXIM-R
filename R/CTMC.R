subdiag <- function(vec, size, offset=0){ 
  M <- matrix(0, size, size)
  M[row(M)-offset == col(M)] <- vec
  return(M)}
ctmctm <- function(R,t){
  ep = 0.00001
  q = max(colSums(t(R)))
  s = dim(R)
  con = exp(-q*t)
  c = con
  B = c*subdiag(1,s[1],0)
  if(q == 0)
    return()
  P = (R + q * subdiag(1,s[1],0) - subdiag(colSums(t(R)),length(colSums(t(R))),0))/q
  A = P
  sumc = c
  k = 1
  while(sumc < 1-ep){
    c = c * q * t/k
    B = B + c * A
    A = A %*% P
    sumc = sumc + c
    k = k+1
  }
  return(B)
}
ctmctd <- function(a,R,t){
  #TODO: check input validity
  y = ctmctm(R,t)
  y = a * y
  return(y)
}
ctmctdplot <- function(a,R,t,j){
  #TODO: check input validity
  mR = dim(R)[1]
  nR = dim(R)[2]
  y = c(a[j])
  NN = 20
  B = ctmctm(R,t/NN)
  for(n in 1:NN){
    a = a %*% B
    y = c(y,a[j])
  }
  plot(seq(0,t,by = t/NN),y, pch = ".")
  lines(seq(0,t,by = t/NN),y)
}
ctmcom <- function(R,t){
  #TODO: check input validity
  ep = 0.00001
  q = max(colSums(t(R)))
  s = dim(R)
  if(q == 0){
    return(t * subdiag(1,s[1],0))
  }
  P = (R + q * subdiag(1,s[1],0) - subdiag(colSums(t(R)),length(colSums(t(R))),0))/q
  A = P
  k = 0
  yek = exp(-q*t)
  ygk = 1 - yek
  sumy = ygk
  B = ygk * subdiag(1,s[1],0)
  while((sumy/q) < t-ep){
    k = k+1
    yek = yek*q*t/k
    ygk = ygk - yek
    B = B + ygk * A
    A = A %*% P
    sumy = sumy + ygk
  }
  B = B / q
  return(B)
}
ctmcod <- function(R){
  d = colSums(t(R))
  Q = R
  m = dim(Q)
  for(i in 1:m[1]){
    Q[i,i] = -d[i]
  }
  Q[,1] = t(t(rep(1,m[1])))
  y = solve(t(Q),c(1,rep(0,m[1]-1)))
  return(y)
}
ctmctc <- function(R,c,t){
  y = ctmcom(R,t)
  y = y * c
  return(y)
}
ctmclrc <- function(R,c){
  y = ctmcod(R) * c
}
ctmcfpt <- function(t,R){
  d = colSums(t(R))
  Q = R
  m = dim(Q)
  for(i in 1:m[1]){
    Q[i,i] = -d[i]
  }
  s = length(t)
  Q = Q[-t,] 
  Q = Q[,-t]
  NT = c(1:m[1])
  NT = NT[-t] 
  y1 = solve(-Q,t(t(rep(1,m[2]-s))))
  y2 = solve(-Q,2*y1)
  returnVal = data.frame(NT = NT,y1 = y1,y2 = y2)
  return(returnVal)
}
#EX6
ex6ar <- function(l){
  #?
  seq1 = c(l,l,l,2*l,2*l,2*l)
  seq2 = c(l,2*l,0,l,2*l,0,l,2*l)
  return(subdiag(seq1,9,3) + subdiag(seq2,9,1))
}
ex6cc <- function(l,m,M,H){
  R = matrix(rep(0,(M+H+1)^2),nrow = M+H+1)
  for (i in 0:(M+H-1)){
    R[i+1,i+2] = l
  }
  for(i in 1:M){
    R[i+1,i] = i * m
  }
  for(i in (M+1):(M+H)){
    R[i+1,i] = M*m
  }
  return(R)
}
ex6cccost <- function(l,m,M,H,r,h){
  c = rep(0,(M+H+1))
  for(i in 0:M){
    c[i+1] = -r * l
  }
  for(i in (M+1):(M+H-1)){
    c[i+1] = h * (i-M) - r * l
  }
  c[M+H+1] = h * H
}
ex6fbd <- function(l,m){
  return(subdiag(l,length(l)+1,-1) + subdiag(m,length(m)+1,1))
}
ex6gms <- function(l,m,N,M){
  R = matrix(0,nrow = N+1,ncol = N+1)
  for(i in 0:(N-1)){
    r[i+1,i+2] = l * min((N-i),M)
  }
  for(i in 1:N){
    R[i+1,i] = i * m
  }
  return(R)
}
ex6gmscost <- function(l,m,N,M,r,dc,rc){
  c = rep(0,N+1)
  for(i in 0:N){
    c[i+1] = rc * min(N-i,M) - i * r + (N-i) * dc
  }
}
ex6inv <- function(l,m,k,r){
  R = subdiag(rep(m,k+r),k+r+1,1) + subdiag(rep(l,k+1),k+r+1,-r)
  return(R)
}
ex6invcost <- function(l,m,k,r,hc,oc,rev){
  c = rep(0,k+r+1)
  c[1] = l * oc
  for(i in 1:(k-r-1)){
    c[i+1] = l * oc + hc * i - rev * m
  }
  for(i in max(k-r,1):(k+r)){
    c[i+1] = hc * i - rem * m
  }
}
ex6lb <- function(l,m,D,M){
  K = M + D
  return(ex66fbd(rep(l,K),rep(m,K)))
}
ex6mfg <- function(l,m,k,K){
  R = subdiag(c(rep(l,k),rep(m,K-k-1)),(2 * K - k),-1) + subdiag(c(rep(m,K-1),rep(0,K-k)),(2 * K - k),1)
  R[2*K-k,k+1] = m
  return(R)
}
ex6mfgcost <- function(l,m,k,K,rev,hc,du){
  c = rep((-m*rev),(2*K-k))
  c[1] = 0
  for(i in 0:K){
    c[i+1] = c[i+1] + i * hc
  }
  for(i in (k+1):(K-1)){
    c[2*K-i+1] = c[2*K-i+1] + i * hc
  }
  c[2*K-k] = c[2*K-k] + m * du
  return(c)
}
ex6ssq <- function(l,m,K){
  return(ex6fbd(rep(l,K),rep(m,K)))
}
ex6ssqcost <- function(l,m,K,hc,lc,bc){
  c = rep(bc,K+1)
  c[1] = 0
  for(i in 0:K){
    c[i+1] = c[i+1] + i * hc
  }
  c[K+1] = c[K+1] + l * lc
  return(c)
}
ex6tel <- function(l,m,K){
  R = matrix(0,nrow = K+1,ncol = K+1)
  R[1,2] = l
  R[K+1,K] = K * m
  for(i in 2:K){
    R[i,i-1] = (i-1) * m
    R[i,i+1] = l
  }
  return(R)
}
ex6telcost <- function(l,m,K,hc,lc){
  c = rep(0,K+1)
  for(i in 0:K){
    c[i+1] = i * hc
  }
  c[K+1] = c[K+1] + l * lc
  return(c)
}