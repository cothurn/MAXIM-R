library(expm)
library(matrixcalc)
dtmctd <- function(a,P,n){
  #TODO: check validity
  #P must be  a square stochastic matrix of size N by N.
  #a must be  a row vector of length N, representing the distribution of the initial state.
  #n is a non-negative integer. 
  return(a%*%(P%^%n))
}
dtmctdplot <- function(a, j, NN, p){
  #TODO: check validity
  #Input: P must be a square stochastic matrix of size N by N.
  #a must be  a row vector of length N, representing the distribution of the initial state. 
  #NN is a non-negative integer. 
  #1 <= j <= N, an integer. 
  y = c()
  for(k in 1:NN)
  {
    y = c(y, dtmctd(a, p, k)[j])
  }
  plot(x = c(1:NN),y = y)
  lines(x = c(1:NN),y = y)
}
dtmcot <- function(P, n){
  #TODO: check validity
  returnVal = diag(dim(P)[1])
  for(k in 1:n){
    returnVal = returnVal + P%^%k
  }
  return(returnVal)
}
dtmcod <- function(P){
  #TODO: check validity
  temp = diag(dim(P)[1]) - P
  temp[,1] = c(rep(1,dim(P)[1]))
  A = c(1,rep(0,dim(P)[1]-1))
  B = temp
  return(solve(t(B),A))
}
dtmctc <- function(P,c,n){
  #TODO: check validity
  tempVal = dtmcot(P, n)
  return(tempVal %*% c)
}
dtmclrc <- function(P,c){
  #TODO: check validity
  return(dtmcod(P) %*% c)
}
dtmcfpt <- function(T,P){
  #TODO: check validity
  s = dim(P)
  sx = s[1]
  x = c(1:sx)
  s = dim(T)
  for( k in 1:s[2]){
    x[T[1,k]] <-  0
  }
  NT = c()
  for(k in 1:sx){
    if(x[k] != 0)
      NT = c(NT, k)
  }
  s = length(NT)
  top = diag(s) - P[NT,NT]
  bot1 = matrix(c(rep(1,s-1),ncol = 1))
  y1 = solve(top,bot1)
  bot2 = 2 * (y1 - bot1)
  y2 = solve(top,bot2)
  return(matrix(c(NT,y1,y1+y2),nrow = s))
}
ex5mr <- function(uu,dd,k,r){
  #TODO: check validity
  r = min(c(k,r))
  U = diag(k+1) * (1-uu) + rbind(cbind(matrix(rep(0,k)),diag(k) * uu),matrix(rep(0,k+1),nrow = 1))
  D = diag(k+1) * (1-dd) + rbind(cbind(matrix(rep(0,k)),diag(k) * dd),matrix(rep(0,k+1),nrow = 1))
  PU = matrix(rep(0,(k+1)^2),nrow = k+1)
  PD = matrix(rep(0,(k+1)^2),nrow = k+1)
  PU[1,1] = 1
  PD[1,1] = 1
  for(i in 1:k+1){
    PU[i,] = PU[i-1,] %*% U
    PD[i,] = PD[i-1,] %*% D
  }
  P = matrix(rep(0,(k+1)^2),nrow = k+1)
  for(i in 0:k){
    for(j in 0:k){
      P[i+1,j+1] = 0
      if(!max(0,j-r,j+i-k)>min(i,j))
      {
        for(l in max(0,j-r,j+i-k):min(i,j)){
          print(c(i+1,j+1,PU[i+1,l+1],PD[min(k-i,r)+1,min(k-i,r)-j+l+1],PU[i+1,m+1] * PD[min(k-i,r)+1,min(k-i,r)-j+l+1]))
          P[i+1,j+1] = P[i+1,j+1] + (PU[i+1,l+1] * PD[min(k-i,r)+1,min(k-i,r)-j+l+1])
        }
      }
    }
  }
  return(P)
}
ex5mrcost <- function(uu,dd,k,r,ru,cd,cbr){
  #TODO: check validity
  #TODO: validate
  c = rep(0,k+1)
  for(i in 0:k){
    c[i+1] = -ru * i + cd * (k-i) + cbr * min(k-i,r)
  }
  return(c)
}
#required for the next function
subdiag <- function(vec, size, offset=0){ 
  M <- matrix(0, size, size)
  M[row(M)-offset == col(M)] <- vec
  return(M)}
ex5inv <- function(s,S,y){
  #TODO: check validity
  #Y is vector
  P = matrix(rep(0,(S-s+1)^2),nrow = S-s+1)
  si = length(y)
  for(i in 0:min(si-1,S-s)){
    P = P + y[i+1] * subdiag(1,S-s+1,i)
  }
  P[S-s+1,S-s+1] = 0
  x = colSums(t(P))
  P[,S-s+1] = t(t(rep(1,S-s+1))) - t(t(x))
  return(P)
}
ex5invcost <- function(s,S,y,hc,ps,oc){
  #TODO: check validity
  #TODO: test this
  c = c(rep(0,S-s))
  for(i in s:S){
    c[i-s+1] = hc*i
    sales = 0
    for(j in 0:(ny-1)){
      sale = sale + min(i,j) * y[j+1]
    }
    c[i-s+1] = c[i-s+1] + ps * sale + oc * sum(y[i-s+1,ny])
  }
  return(c)
}
ex5manp <- function(p,l,a){
  #TODO: check validity
  #TODO: test(should be fine)
  #a is a vector
  si = length(a)
  P = matrix(rep(0,si^2),nrow = si)
  P = subdiag(p[1:si-1],si,-1) + subdiag(1-l-p,si,0) + t(t(l)) %*% t(a)
  return(P)
}
ex5manpcost <- function(p,l,a,s,b,d,t){
  #TODO: check validity
  #TODO: test
  na = length(a)
  c = c(rep(0,si))
  for(i in i:na){
    c[i] = s[i] + b[i]*p[i] + (d[i] + t(t)*t(t(a))) * l[i]
  }
  return(c)
}
ex5mfg <- function(A,B,a1,a2){
  #TODO: check validity
  P = matrix(c(rep(0,(A+B+1)^2)),nrow = A+B+1)
  P = (1-a1) * a2 * subdiag(rep(1,A+B),A+B+1,1) + ((a1*a2) + (1-a1)*(1-a2)) * subdiag(rep(1,A+B+1),A+B+1,0) + (1-a2) * a1 * subdiag(rep(1,A+B),A+B+1,-1)
  P[1,1] = 1-a1
  P[1,2] = a1
  P[A+B+1,A+B+1]=1-a2
  P[A+B+1,A+B] = a2
  return(P)
}
ex5mfgcost <- function(A,B,a1,a2,r,du,hA,hB){
  #TODO: check validity
  c = rep(0,A+B+1)
  c[1] = hB*B-r*a1+du*a1
  for(i in 2:B){
    c[i] = hB*(B+1-i)-r*a1
  }
  c[B+1] = -r*a1*a2
  for(j in (B+2):(A+B)){
    c[j] = hA*(j-B-1)-r*a2
  }
  c[A+B+1] = hA*A - r*a2 + du*a2
  return(c)
}
ex5stock <- function(L,U){
  #TODO: check validity
  d = U - L + 1
  P = matrix(0,nrow = d, ncol = d)
  P = 0.2 * (subdiag(1,d,0) + subdiag(1,d,1) + subdiag(1,d,-1) + subdiag(1,d,2) + subdiag(1,d,-2))
  P[1:2,1] = c(0.6,0.4)
  P[((d-1):d),d] = c(0.4,0.6)
  return(P)
}
ex5tel <- function(K,a){
  #TODO:check validity
  P = subdiag(a[1],K+1,1)
  si = length(a)
  for(i in 0:(min(K,si-2))){
    P = P + a[i+2] * subdiag(1,k+1,-i)
  }
  P[1,1:(min(K+1,si))] = a[1:(min(K+1,si))]
  P[,K+1] = 0
  x = colSums(t(P))
  P[,K+1] = 1 - t(x)
  return(P)
}
ex5telcost <- function(K,a,rt,cl){
  #TODO: check validity
  c = rep(0,K+1)
  na = length(a)
  for(r in K:(na-1)){
    c[1] = c[1] + (r-K)*a[r+1]
  }
  c[1] = cl * c[1]
  for(i in 1:K){
    for(r in (K+1-i):(na-1)){
      c[i+1] = c[i+1] + (r-K-1+i) * a[r+1]
    }
    c[i+1] = cl * c[i+1] - rt
  }
  return(c)
}
ex5wea <- function(){
  return(matrix(c(.5,.3,.2,.5,.2,.3,.4,.5,.1),nrow = 3,byrow = TRUE))
}