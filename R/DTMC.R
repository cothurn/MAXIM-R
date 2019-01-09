checkaP <- function(a,P){
  sA = length(a)
  sP = dim(P)
  if(sum(a>0)!=sA)
    print("initial probabilities must be non-negative")
  if(abs(sum(a))-1 > 10^-12)
    print("initial probabiliies must sum to 1")
  if(sum(rowSums(P)!=1)!=0)
    print("row sums of the transition matrix must be 1")
  if(sP[1]!=sP[2])
    print("transition matrix must be square")
  if(sum(P>0)!=sP[1]*sP[2])
    print("transition probabilities must be non-negative")
  if(sA != sP[1])
    print("the initial distribution is not compatible with the transition matrix")
}
checkP <- function(P){
  sP = dim(P)
  if(sp[1]!=sp[2])
    print("transition matrix must be square")
  if(sum(P>0)!=sP[1]*sP[2])
    print("transition probabilities must be non-negative")
  if(sum(rowSums(P)!=1)!=0)
    print("row sums of the transition matrix must be 1")
  
}
checkcP <- function(c,P){
  sC = length(c)
  sP = dim(P)
  if(sC!=sP[2])
    print("cost vector not compatible witht he transition matrix")
}

#' Transient Distribution
#' @param a a row vector of length N
#' @param P a square stochastic matrix of size N by N.
#' @param n a non-negative integer
#' @return y(i) = P(X_n = i),  1 <= i <= N, where X_n, n >= 0 is a DTMC with transition probability matrix P and initial distribution a.
#' @examples
#'dtmctd(c(1,1),matrix(c(0.5,0.5,0.5,0.5),nrow = 2),5)
dtmctd <- function(a,P,n){
  if(n<0 || n%%1!=0)
    print("Invalid entry of n")
  checkaP(a,P)
  return(a%*%(P%^%n))
}

#' Time Plot of Transient Distribution
#' @param a a row vector of length N
#' @param j a integer 1 <= j <= N
#' @param NN a non-negative integer
#' @param P a square stochastic matrix of size N by N.
#' @return y(i) = P(X_n = i),  1 <= i <= N, where X_n, n >= 0 is a DTMC with transition probability matrix P and initial distribution a.
#' @examples
#'dtmctdplot(c(1,1),2,5,matrix(c(0.5,0.5,0.5,0.5),nrow = 2))
dtmctdplot <- function(a, j, NN, p){
  if(j < 1 || j > length(a))
    print('invalid entry for j')
  if(NN < 0 || NN%%1 != 0)
    print("invalid entry for NN")
  checkaP(a,P)
  y = c()
  for(k in 1:NN)
  {
    y = c(y, dtmctd(a, p, k)[j])
  }
  plot(x = c(1:NN),y = y)
  lines(x = c(1:NN),y = y)
}

#' Occupancy Times
#' @param P a square stochastic matrix of size N by N.
#' @param n a nonnegative integer
#' @return M(i,j) = the expected number of visits to state j starting from state i by a DTMC X_n, n >= 0 with transition probability matrix P.
#' @examples
#' dtmcot(matrix(c(0.5,0.5,0.5,0.5),nrow = 2), 5)
dtmcot <- function(P, n){
  if(n<0 || n %% 1 != 0)
    print("invalid entry for n")
  checkP(P)
  returnVal = diag(dim(P)[1])
  for(k in 1:n){
    returnVal = returnVal + P%^%k
  }
  return(returnVal)
}

#' Limiting Distribution and Occupancy Distribution.
#' @param P an aperiodic irreducible square stochastic matrix of size N by N. 
#' @return y(i) = limit as n goes to infinity of P(X_n = i),  1 <= i <= N, where X_n, n >= 0 is an irreducible aperiodic DTMC with transition probability matrix P. If the DTMC is periodic y is the occupancy distribution of X
#' @examples
#' dtmcod(matrix(c(0.5,0.5,0.5,0.5),nrow = 2))
dtmcod <- function(P){
  checkP(P)
  temp = diag(dim(P)[1]) - P
  temp[,1] = c(rep(1,dim(P)[1]))
  A = c(1,rep(0,dim(P)[1]-1))
  B = temp
  return(solve(t(B),A))
}

#' Total Cost
#' @param P a square stochastic matrix of size N by N
#' @param c vector of length N
#' @param n a nonnegative integer
#' @return y(i) is the total expected cost incurred over time 0 through n staring in state i, for a DTMC X_n, n >= 0 with a transition probability matrix P, and that incurs an expected cost c(i) every time it visits state i.
#' @examples
#' dtmctc(matrix(c(0.5,0.5,0.5,0.5),nrow = 2),c(1,2),5)
dtmctc <- function(P,c,n){
  if(n<0 || n %% 1 != 0)
    print("invalid entry for n")
  checkcP(c,P)
  tempVal = dtmcot(P, n)
  return(tempVal %*% c)
}

#' Long-run Cost Rate
#' @param P an irreducible square stochastic matrix of size N by N
#' @param c vector of length N
#' @return y is the long-run expected cost per unit time for a DTMC X_n, n >= 0 with a transition probability matrix P, and that incurs an expected cost c(i) every time it visits state i.
#' @examples
#' dtmclrc(matrix(c(0.5,0.5,0.5,0.5),nrow = 2),c(1,2))
dtmclrc <- function(P,c){
  checkcP(c,P)
  return(dtmcod(P) %*% c)
}

#' First Passage Times
#' @param P a stochastic matrix of size N by N 
#' @param T row vector representing the set of target states
#' @return y=[y0 y1 y2], where y0, y1, and y2 are column vectors.  [y1(i) y2(i)] is the [mean,  second moment] of the first passage time to visit any of the states in the target set of states T, starting in a non-target state y0(i) for a DTMC X_n, n >= 0 with a transition probability matrix P.
#' @examples
#' dtmcfpt(c(1,2),matrix(c(rep(1/3,9)),nrow = 3))
dtmcfpt <- function(Tr,P){
  checkP(P)
  s = dim(P)
  sx = s[1]
  x = c(1:sx)
  s = length(Tr)
  for(k in 1:s){
    x[Tr[k]] <-  0
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

#' Machine Reliability
#' @param uu P(up|up)
#' @param dd P(down|down)
#' @param k number of machines
#' @param r number of repair persons
#' @return P = transition  probability matrix for the machine reliability problem. (See Example 5.4)
#' @examples
#' ex5mr(0.5,0.5,3,1)
ex5mr <- function(uu,dd,k,r){
  if(uu < 0 || uu > 1)
    print("invalid entry for uu")
  if(dd < 0 || dd > 1)
    print("invalid entry for dd")
  if(k<0 || k%%1 != 0)
    print("invalid entry for k")
  if(r < 0 || r%%1 != 0)
    print("invalid entry for r")
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
          #print(c(i+1,j+1,PU[i+1,l+1],PD[min(k-i,r)+1,min(k-i,r)-j+l+1],PU[i+1,j+1] * PD[min(k-i,r)+1,min(k-i,r)-j+l+1]))
          P[i+1,j+1] = P[i+1,j+1] + (PU[i+1,l+1] * PD[min(k-i,r)+1,min(k-i,r)-j+l+1])
        }
      }
    }
  }
  return(P)
}

#' Machine Reliability Cost Model
#' @param uu P(machine is up tomorrow|it is up today)
#' @param dd P(machine is down tomorrow|it is down today)
#' @param k number of machines
#' @param r number of repair persons
#' @param ru per day revenue of an up machine
#' @param cd per day cost of a down machine
#' @param cbr per day cost of a busy repair person. 
#' @return c(i) = one day cost if i machines are working at the beginning of the day.
#' @examples
#' ex5mrcost(0.5,0.5,3,1,10,5,3)
ex5mrcost <- function(uu,dd,k,r,ru,cd,cbr){
  if(uu < 0 || uu > 1)
    print("invalid entry for uu")
  if(dd < 0 || dd > 1)
    print("invalid entry for dd")
  if(k<0 || k%%1 != 0)
    print("invalid entry for k")
  if(r < 0 || r%%1 != 0)
    print("invalid entry for r")
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

#' Inventory Systems
#' @param s base stock level (stock not allowed to go below this level)
#' @param S Restocking level
#' @param y row vector of the pmf of weekly demand
#' @return P = transition probability matrix for the inventory system problem. (See Example 5.6)
#' @examples ex5inv(5,10,c(0.4,0.3,0.2,0.1))
ex5inv <- function(s,S,y){
  if(s<0 || s%%1 != 0)
    print("invalid entry for s")
  if(S<0 || S%%1 != 0)
    print("invalid entry for S")
  if(sum(y<0) != 0 || sum(y) > 1.000000000000001)
    print("invalid entry for y")
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

#' Inventory Systems Cost Model
#' @param s base stock level (stock not allowed to go below this level); 
#' @param S Restocking level;
#' @param y row vector of the pmf of weekly demand;
#' @param hc cost of holding one item for one unit of time;
#' @param ps profit from selling one item;
#' @param oc cost of placing an order;
#' @return c(i) = expected cost in the current period if the inventory at the beginning is i.
#' @examples
#' ex5invcost(5,10,c(0.4,0.3,0.2,0.1),1,5,2)
ex5invcost <- function(s,S,y,hc,ps,oc){
  if(s<0 || s%%1 != 0)
    print("invalid entry for s")
  if(S<0 || S%%1 != 0)
    print("invalid entry for S")
  if(sum(y<0) != 0 || sum(y) > 1.000000000000001)
    print("invalid entry for y")
  
  c = c(rep(0,S-s))
  for(i in s:S){
    c[i-s+1] = hc*i
    sales = 0
    for(j in 0:(length(y)-1)){
      sales = sales + min(i,j) * y[j+1]
    }
    c[i-s+1] = c[i-s+1] - ps * sales
    if(!is.na(sum(y[(i-s+1):length(y)])))
      c[i-s+1] = c[i-s+1] + oc * sum(y[(i-s+1):length(y)])
  }
  return(c)
}

#' Manpower Systems
#' @param p a row vector, p(i) =  probability of promotion from grade i to i+1
#' @param l a row vector of the same length as p,  l(i) =  probability of leaving from grade i
#' @param a row vector, a(i) = probability that a new employee joins grade i.
#' @return P = transition  probability matrix for the manpower system problem. (See Example 5.8)
#' @examples
#' ex5manp(c(0.4,0.4,0.2,0),c(0.4,0.3,0.2,0.1),c(0.4,0.3,0.2,0.1))
ex5manp <- function(p,l,a){
  if(length(p)!=length(l) || length(l)!=length(a))
    print("a, l, and p must have the same length")
  if(sum(p<0)!=0 || sum(p>1) != 0)
    print("invalid value for p")
  if(p[length(p)]!=0)
    print("last entry of p must be zero")
  if(sum(l<0)!=0 || sum(l>1) != 0)
    print("invalid value for l")
  if(sum((l+p)>1)!=0){
    print("invalid value for l or p")
  if(sum(a<0)!=0 || abs(sum(a)-1) > 10^(-12))
    print("invalid value for a")
  }
  si = length(a)
  P = matrix(rep(0,si^2),nrow = si)
  P = subdiag(p[1:si-1],si,-1) + subdiag(1-l-p,si,0) + t(t(l)) %*% t(a)
  return(P)
}

#' Manpower Systems Cost Model
#' @param p a row vector, p(i) =  probability of promotion from grade i to i+1. The last element of p must be zero. 
#' @param l a row vector of the same length as p, l(i) =  probability of leaving from grade i.
#' @param a row vector, a(i) = probability that a new employee joins grade i. Must be a valid pmf.
#' @param s a row vector, s(i) = salary of an employee in grade i.
#' @param b a row vector, b(i) = bonus for promotion from grade i to i+1. The last element of b must be zero.
#' @param d a row vector, d(i) = cost of an employee departing from grade i.
#' @param t a row vector, t(i) = cost of training an employee staring in grade i.
#' @return c = a row vector, c(i) = expected one-period cost in state i.
#' @examples
#' ex5manpcost(c(0.4,0.4,0.2,0),c(0.4,0.3,0.2,0.1),c(0.4,0.3,0.2,0.1),c(1,2,3,4),c(1,2,3,0),c(1,2,3,4),c(1,2,3,4))
ex5manpcost <- function(p,l,a,s,b,d,tr){
  if(max(length(p),length(l),length(a),length(s),length(b),length(d),length(t)) != min(length(p),length(l),length(a),length(s),length(b),length(d),length(t)))
    print("a, l, p, s, b, d, and t must have the same length")
  if(sum(p<0)!=0 || sum(p>1) != 0)
    print("invalid value for p")
  if(sum(l<0)!=0 || sum(l>1) != 0)
    print("invalid value for l")
  if(sum(a<0)!=0 || abs(sum(a)-1) > 10^(-12))
    print("invalid value for a")
  na = length(a)
  c = c(rep(0,na))
  for(i in 1:na){
    c[i] = s[i] + b[i]*p[i] + (d[i] + t(tr)%*%t(t(a))) * l[i]
  }
  return(c)
}

#' Manufacturing Systems
#' @param A size of the bin for machine 1;
#' @param B size of the bin for machine 2;
#' @param a1 prob(non-defective) for machine 1;
#' @param a2 prob(non-defective) for machine 2.
#' @return P = transitrion probability matrix for the manufacturing system. (See Example 5.7.)
#' @examples
#' ex5mfg(5,3,0.5,0.3)
ex5mfg <- function(A,B,a1,a2){
  if(A<0 || A %% 1 != 0)
    print("invalid value for A")
  if(B<0 || B %% 1 != 0)
    print("invalid value for B")
  if(a1 >1 || a1 < 0)
    print("invalid value for a1")
  if(a2 >1 || a2 < 0)
    print("invalid value for a2")
  P = matrix(c(rep(0,(A+B+1)^2)),nrow = A+B+1)
  P = (1-a1) * a2 * subdiag(rep(1,A+B),A+B+1,1) + ((a1*a2) + (1-a1)*(1-a2)) * subdiag(rep(1,A+B+1),A+B+1,0) + (1-a2) * a1 * subdiag(rep(1,A+B),A+B+1,-1)
  P[1,1] = 1-a1
  P[1,2] = a1
  P[A+B+1,A+B+1]=1-a2
  P[A+B+1,A+B] = a2
  return(P)
}

#' Manufacturing Systems Cost Model
#' @param A size of bin for machine 1;
#' @param B size of bin for machine 2;
#' @param a1 prob(non-defective) for machine 1;
#' @param a2 prob(non-defective) for machine 2.
#' @param r revenue from a complete assembly;
#' @param du cost of turning a machine on;
#' @param hA cost of holding an item in bin A for one period;
#' @param hB cost of holding an item in bin B for one period;
#' @return c = row vector, c(i) = expected one-period cost in state i.
#' @examples
#' ex5mfgcost(5,3,0.5,0.3,10,1,2,1)
ex5mfgcost <- function(A,B,a1,a2,r,du,hA,hB){
  if(A<0 || A %% 1 != 0)
    print("invalid value for A")
  if(B<0 || B %% 1 != 0)
    print("invalid value for B")
  if(a1 >1 || a1 < 0)
    print("invalid value for a1")
  if(a2 >1 || a2 < 0)
    print("invalid value for a2")
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

#' Stock Market
#' @param L Lower bound for the stock price,
#' @param U Upper bound for the stock price.
#' @return P = transition  probability matrix for the stock market.
#' @examples
#' ex5stock(1,10)
ex5stock <- function(L,U){
  if(L < 0 || L %% 1 > 0)
    print("invalid value for L")
  if(U < 0 || U %% 1 > 0)
    print("invalid value for U")
  #TODO: check validity
  d = U - L + 1
  P = matrix(0,nrow = d, ncol = d)
  P = 0.2 * (subdiag(1,d,0) + subdiag(1,d,1) + subdiag(1,d,-1) + subdiag(1,d,2) + subdiag(1,d,-2))
  P[1:2,1] = c(0.6,0.4)
  P[((d-1):d),d] = c(0.4,0.6)
  return(P)
}

#' Telecommunications
#' @param K buffer capacity
#' @param a = row vector, a(i) =  p(i-1 packets arrive during one time slot)
#' @return P = transition  probability matrix for the Telecommunications system. (See  Example 5.10.)
#' @examples
#' ex5tel(5,c(0.1,0.2,0.3,0.3,0.1))
ex5tel <- function(K,a){
  if(K < 0 || K %% 1 != 0)
    print("invalid entry for K")
  if(any(a<0) || sum(a) > 1)
    print("invalid entry for a")
  P = subdiag(a[1],K+1,1)
  si = length(a)
  for(i in 0:(min(K,si-2))){
    P = P + a[i+2] * subdiag(1,K+1,-i)
  }
  P[1,1:(min(K+1,si))] = a[1:(min(K+1,si))]
  P[,K+1] = 0
  x = colSums(t(P))
  P[,K+1] = 1 - t(x)
  return(P)
}

#' Telecommunications Cost Model
#' @param K = buffer capacity,
#' @param a = row vector, a(i) =  p(i-1 packets arrive during one time slot). a must be a valid pmf.
#' @param rt =  revenue from transmitting a single packet,
#' @param cl =  cost of losing a single packet.
#' @return c = column vector, c(i) =  expected cost in one slot if there are i-1 packets in the buffer at the end of the previous slot, 1 <= i <= K+1.
#' @examples
#' ex5telcost(5,c(0.1,0.2,0.3,0.3,0.1),10,1)
ex5telcost <- function(K,a,rt,cl){
  if(K < 0 || K %% 1 != 0)
    print("invalid entry for K")
  if(any(a<0) || sum(a) > 1.0001)
    print("invalid entry for a")
  c = rep(0,K+1)
  na = length(a)
  if(K <= (na-1)){
    for(r in K:(na-1)){
      c[1] = c[1] + (r-K)*a[r+1]
    }
  }
  c[1] = cl * c[1]
  for(i in 1:K){
    for(r in (K+1-i):(na-1)){
      c[i+1] = c[i+1] + (r-K-1+i) * a[r+1]
    }
    if(is.na(c[i+1])){
      c[i+1] = 0
    }
    c[i+1] = cl * c[i+1] - rt
  }
  return(c)
}

#' Weather Model
#' @return The 3X3 matrix of the weather model of Example 5.5
#' @examples
#' ex5wea()
ex5wea <- function(){
  return(matrix(c(.5,.3,.2,.5,.2,.3,.4,.5,.1),nrow = 3,byrow = TRUE))
}