smpod <- function(P,w){
  #TODO: not tested
  pi = dtmcod(P)
  y = pi * w
  y = y %*% solve(colSums(y))
  return(y)
}
smplrc <- function(P,w,c,d){
  pi = smpod(P,w)
  y = pi %*% (c + d/w)
}
ex7ser <- function(N,v,tau){
  P = matrix(rep(0,(N+1)^2),nrow = N+1)
  P[1,2:(N+1)] = t(v)/colSums(v)
  P[2:(N+1),1]=1
  w = matrix(c(1/(colSums(v)),tau),ncol = 1)
  return(list(P,w))
}

mm1k <- function(l,m,K,i){
  R = ex6ssq(l,m,K)
  yy = ctmcod(R)
  y =switch(i,yy%*%t(t(c(0:K))),
            yy%*%t(t(c(0,0:(K-1)))),
            yy%*%t(t(c(0:K)))/(l*(1-yy[K+1])),
            yy%*%t(t(c(0,0:(K-1))))/(l*(1-yy[K+1])),
            yy,
            yy[1:K]/(1-yy[K+1])
  )
  return(y)
}
mmsk <- function(l,m,s,K,i){
  R = ex6cc(l,m,s,K-s)
  yy = ctmcod(R)
  y = switch(i,
             yy%*%t(t(c(0:K))),
             yy[(s+1):(K+1)]%*%t(t(c(0:(K-s)))),
             yy%*%t(t(c(0:K)))/(l*(1-yy[K+1])),
             yy[(s+1):(K+1)]%*%t(t(c(0:(K-s))))/(l*(1-yy[K+1])),
             yy,
             yy[1:K]/(1-yy[K+1]),
             l * (1-yy[K+1])/m,
             yy%*%t(t(c(rep(0,s),rep(1,K-s+1)))))
  return(y)
}
mm1 <- function(l,m,i){
  rho = l/m
  if(rho >= 1)
  {
    print("Queue unstable")
    return()
  }
  y = switch(i,
         rho/(1-rho),
         (rho^2)/(1-rho),
         1/(m-l),
         rho/(m-l),
         geometricpmf((1-rho), max(floor((-10-log(1-rho))/log(rho)),5))
         )
  return(y)
}
mms <- function(l,m,s,i){
  rho = l/(s*m)
  yy = exp(l/m) * poissonpmf(l/m,s)
  normf = sum(yy[1:s]) + yy[s+1]/(1-rho)
  yy = yy/normf
  y = switch(i,
             l/m + yy[s+1]*rho/(1-rho)^2,
             yy[s+1]*rho/(1-rho)^2,
             1/m + yy[s+1]/(s*m*(1-rho)^2),
             yy[s+1]/(s*m*(1-rho)^2),
             c(yy[1:s], yy[s+1]*rho^c(0:max(0,(max(s,floor((-10-log(yy[s+1]))/log(rho))))-s))),
             l/m
             )
  return(y)
}
mminf <- function(l,m,i){
  y = switch(i,
             l/m,
             0,
             1/m,
             0,
             poissonpmf(l/m,floor(-15/log(l/m))),
             l/m
             )
  return(y)
}
mg1 <- function(l,m,s2,i){
  rho = l * m
  y = switch(i,
             rho + .5*l^2*(s2+m^2)/(1-rho),
             .5*l^2*(s2+m^2)/(1-rho),
             m+ .5*l*(s2+m^2)/(1-rho),
             .5*l*(s2+m^2)/(1-rho),
             rho
             )
  return(y)
}
funeq <- function(i,x,m){
  y = "error"
  if(i == 1){
    if(x < 0)
      print('invalid value of x')
    else if(m <= x)
      print('Queue is unstable')
    else
      y = x/m
  }
  if(i == 2){
    mx = dim(x)[0]
    nx = dim(x)[1]
    if (mx !=1||nx != 2)
      print('invalid entry for x or i')
    else if(x[1] < 1 || x[1] - floor(x[1]) != 0)
      print('x(1) must be a positive integer')
    else if(x[2] < 0)
      print('x(2) must be positive')
    else if(m <= x[2]/x[1])
      print('Queue is unstable')
    else
      y0 = 0
      y1 = (x[2]/(m+x[2]))^x[1]
      while(abs(y1-y0) > .00001){
        y0 = y1;
        y1 = (x[2]/(m*(1-y0)+x[2]))^x[1]
      }
      y=y1
  }
  if(i == 3){
    if( x[1] < 1 || x[1] - floor(x[1]) != 0 )
      print('x[1] must be a positive integer')
    mx = dim(x)[0]
    nx = dim(x)[1]
    if (mx != 1 || nx != 2*x[1]+1)
      print('invalid entry for x')
    l=x[2:(x[1]+1)]
    p=x[(x[1]+2):(2*x[1]+1)];
    if (any(l < 0)||any(p < 0)||abs(sum(p) - 1.0) > 10^(-8))
      msgbox('invalid entry for x')
    else if(m <= 1/sum(p/l))
      print('Queue is unstable')
    else
      y0 = 0
      y1 = sum(p*(l/(m+l)));
      while(abs(y1-y0) > .00001){
        y0 = y1;
        y1 = sum(p*(l/(m*(1-y0)+l)))
      }
      y=y1
  }
  if(i == 4){
    mx = dim(x)[0]
    nx = dim(x)[1]
    if(mx !=1 || nx != 1)
      print('invalid entry for x or i')
    if(x < 0)
      print('x must be  a positive scalar')
    else if(m <= 1/x)
      print('Queue is unstable')
    else
      y0 = 0
      y1 = exp(-m*x)
      while(abs(y1-y0) > .00001){
        y0 = y1;
        y1 = exp(-m*(1-y0)*x)
      }
      y=y1
  }
  if(i == 5){
    mx = dim(x)[0]
    nx = dim(x)[1]
    if(mx != 1)
      print('x must be row vector')
    if(any(x < 0) || abs(sum(x) - 1.0) > 10^(-10))
      print('x must be a valid pmf')
    else if (m <= 1/sum(x*c(0:(nx-1))))
      print('Queue is unstable')
    else
      y0 = 0
      y1 = sum(x*exp(-m*c(0:(nx-1))));
      while(abs(y1-y0) > .00001){
        y0 = y1
        y1 = sum(x*exp(-m*(1-y0)*c(0:(nx-1))));}
      y=y1
  }
  if(i == 6){
    if(x[1] < 0||x[1] > x[2])
      print('invalid entry for x')
    else if(m <= 2/(x[1]+x[2]))
      print('Queue is unstable')
    else
      y0 = 0
      y1 = (exp(-x[1]*m) - exp(-x[2]*m))/(m*x[2]-x[1])
      while(abs(y1-y0) > .00001){
        y0 = y1
        y1 = (exp(-x[1]*m*(1-y0)) - exp(-x[2]*m*(1-y0)))/(m*(1-y0)*(x[2]-x[1]))
      }
      y=y1
  }
  return(y)
}
gm1 <- function(l,m,a,i){
  rho = l/m
  y = switch(i,
             rho/(1-a),
             rho*a/(1-a),
             1/(m*(1-a)),
             a/(m*(1-a)),
             c(1-rho,rho*geometricpmf(1-a,floor((-10-log(1-a))/log(a)))),
             geometricpmf(1-a,floor((-10-log(1-a))/log(a)))
  )
  return(y)
}
qn <- function(l,m,s,RM,i){
  si = dim(RM)
  a = l * solve(subdiag(1,si[2],0)-RM)
  if(i == 0)
    return(a)
  if(any(a > s%*%m))
    print("network is unstable")
    return()
  y = c()
  for(j in 1:si[2]){
    y = c(y,mms(a[j],m[j],s[]))
  }
  return(y)
}
