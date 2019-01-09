
## Probability
MAXIM can compute the distributions of discrete and continuous random variables as described below.
### Binomial Cumulative Distribution Function
Usage: y=bincdf(n,p)
Input: n a non-negative integer;  0 <= p <= 1.
Output: y(i) = P(X <= i-1), i=1,2,...,n+1, where X is a  Binomial(n,p) rv.

### Binomial Probability Mass Function
Usage:  y=binpmf(n,p)
Input: n: a non-negative integer; 0 <= p <= 1.
Output: y(i) = P(X = i-1),  i=1,...n+1, where X  is a  Binomial(n,p) rv.

### Erlang Cumulative Distribution Function
Usage:  y=erlangcdf(k,l,x)
Input: k >= 1, integer; l >= 0, x is a row vector.
Output: y(i) =  F(x(i)), where F is an Erlang(k,l) cdf. 

### Plot of Erlang Cumulative Distribution Function
Usage:  erlangcdfplot(k,l,x),
Input: k >= 1, integer; l >= 0; x is a row vector.
Output: Plot of the cdf of Eralng(k,l) (y vs. x).
 
### Eralng Probability Density Function
Usage:  y=erlangpdf(k,l,x)
Input: k >= 1, integer; l >= 0, x is a row vector.
Output: y(i) =f(x(i)), where f is an Erlang(k,l) pdf.

### Plot of Erlang Probability Density Function
Usage:  erlangpdfplot(k,l,x),
Input: k >= 1, integer; l >= 0, x is a row vector.
Output: Plot of the pdf of Eralng(k,l) (y vs. x).

### Exponential Cumulative Distribution Function
Usage:  y=expcdf(l,x)
Input: l >= 0, x is a non-negative row vector.
Output: y(i) =  F(x(i)), where F is an Exp(l) cdf. 

### Plot of Exponential Cumulative Distribution Function
Usage:  expcdfplot(l,x),
Input: l >= 0, x is a row vector.
Output: Plot of the cdf of exp(l) (y vs. x). 

### Exponential Probability Density Function
Usage:  y=exppdf(l,x)
Input:  l >= 0 , x is a row vector.
Output: y(i) =f(x(i)), where f is an Exp(l) pdf.

### Plot of Exponential Probability Density Function
Usage:  exppdfplot(l,x),
Input: l >= 0, x is a row vector.
Output: Plot of the pdf of exp(l) (y vs. x). 

### Geometric Cumulative Distribution Function
Usage:  y=geometriccdf(p,k)
Input: k >= 1, integer; 0 <= p <= 1.
Output: y(i) = P(X <= i),  i=1,2,...,k, where X  is a Geometric(p) rv.

### Geometric Probability Mass Function
Usage:  y=geometricpmf(p,k)
Input: k >= 1, integer; 0 <= p <= 1.
Output: y(i) = P(X = i),  i=1,...k, where X is a  Geometric(p) rv.

### Negative Binomial Cumulative Distribution Function
Usage:  y=negbincdf(r,p,k)
Input: r,k >= 1, integer; 0 <= p <= 1.
Output: y(i) = P(X <= i),  i=r,r+1,...,r+k, where X  is a Negative Binomial(r,p) rv.

### Negative Binomial Probability Mass Function
Usage:  y=negbinpmf(r,p,k)
Input: r,k >= 1, integer; 0 <= p <= 1.
Output: y(i) = P(X = i),  i=r,r+1,...,r+k, where X  is a  Negative Binomial(r,p) rv.

###  Normal Cumulative Distribution Function
Usage:  y=normalcdf(m,s,x)
Input: s >= 0, x is a row vector.
Output: y(i) =  F(x(i)), where F is an Normal(m,s) cdf, with mean  m and variance s. 

### Plot of Normal Cumulative Distribution Function
Usage:  normalcdfplot(m,s,x),
Input: s >= 0, x is a row vector.
Output: Plot of the cdf of Normal(m,s) (y vs. x). 

### Normal Probability Density Function
Usage:  y=normalpdf(m,s,x)
Input:  s >= 0 , x is a row vector.
Output: y(i) =f(x(i)), where f is an Normal(m,s) pdf, with mean m and variance s. 

### Plot of Normal Probability Density Function
Usage:  normalpdfplot(m,s,x)
Input: s >= 0, x is a row vector
Output: Plot of the pdf of Normal(m,s) (y vs. x). 

### Poisson Cumulative Distribution Function
Usage:  y=poissoncdf(l,k)
Input: 0 <= l <= 700.
Output: y(i) = P(X <= i-1),  i=1,2,...,k+1, where X is a  Poisson(l) rv.
 
### Poisson Probability Mass Function
Usage:  y=poissonpmf(l,k)
Input:  0 <= l <= 700.
Output: y(i) = P(X = i-1),  i=1,...k+1, where X  is a Poisson(l) rv.
