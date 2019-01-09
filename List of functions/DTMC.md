## Discrete Time Markov Models

MAXIM includes 7 functions for the computational analysis of Discrete Time Markov Chains. It also includes 12 functions that generate special Discrete Time Markov Models.

### Transient Distribution.
Usage: y = dtmctd(a,n,P).

Input: 

P must be a square stochastic matrix of size N by N.

a must be a row vector of length N, representing the distribution of the initial state.

n is a non-negative integer. 

Output: y(i) = P(X_n = i), 1 <= i <= N, where X_n, n >= 0 is a DTMC with transition probability matrix P and initial distribution a.

## Time Plot of Transient Distribution.
Usage: y = dtmctdplot(a,j,NN,P).

Input: 

P must be a square stochastic matrix of size N by N.

a must be a row vector of length N, representing the distribution of the initial state.

NN is a non-negative integer. 1 <= j <= N, an integer. 

Output: y(n) = P(X_n-1 = j), 1 <= n <= NN+1, where X_n, n >= 0 is a DTMC with transition probability matrix P and initial distribution a. Also produces a plot y against the time axis.

## Occupancy Times.
Usage: M = dtmcot(P,n).

Input: 

P must a square stochastic matrix of size N by N. 

n >= 0 an integer. 

Output: M(i,j) = the expected number of visits to state j starting from state i by a DTMC X_n, n >= 0 with transition probability matrix P.

## Limiting Distribution and Occupancy Distribution.
Usage: y = dtmcod(P).


Input:

P must an aperiodic irreducible square stochastic matrix of size N by N. 

Output: y(i) = limit as n goes to infinity of P(X_n = i), 1 <= i <= N, where X_n, n >= 0 is an irreducible aperiodic DTMC with transition probability matrix P. If the DTMC is periodic y is the occupancy distribution of \X.

## Total Cost.
Usage: y = dtmctc(P,c,n).

Input: 

P must a square stochastic matrix of size N by N.

c is row vector of length N.

n >= 0 is an integer.  

Output: y(i) is the total expected cost incurred over time 0 through n staring in state i, for a DTMC X_n, n >= 0 with a transition probability matrix P, and that incurs an expected cost c(i) every time it visits state i.
 
## Long-run Cost Rate.
Usage: y = dtmclrc(P,c).

Input: 
P must an irreducible square stochastic matrix of size N by N

c is row vector of length N. 

Output: y is the long-run expected cost per unit time for a DTMC X_n, n >= 0 with a transition probability matrix P, and that incurs an expected cost c(i) every time it visits state i.
 
## First Passage Times.
Usage: y = dtmcfpt(T,P).

Input: 

P must a stochastic matrix of size N by N

T is row vector representing the set of target states. 

Output: y=[y0 y1 y2], where y0, y1, and y2 are column vectors. [y1(i) y2(i)] is the [mean, second moment] of the first passage time to visit any of the states in the target set of states T, starting in a non-target state y0(i) for a DTMC X_n, n >= 0 with a transition probability matrix P.

## Machine Reliability
Usage: P=ex5mr(uu,dd,k,r).
Input: 
uu = P(up|up); 

dd = P(down|down);

k = number of machines;

r = number of repair persons.

Output: P = transition probability matrix for the machine reliability problem. (See Example 5.4)
 
## Machine Reliability Cost Model
Usage: c=ex5mrcost(uu,dd,k,r,ru,cd,cbr)

Input: 

uu = P(machine is up tomorrow|it is up today); 

dd = P(machine is down tomorrow|it is down today);

k = number of machines;

r = number of repair persons;

ru = per day revenue of an up machine;

cd = per day cost of a down machine;

cbr = per day cost of a busy repair person. 

Output: c(i) = one day cost if i machines are working at the beginning of the day.

## Inventory Systems
Usage: P=ex5inv(s,S,y).

Input: 
s = base stock level (stock not allowed to go below this level); 

S = Restocking level;

y = row vector of the pmf of weekly demand.

Output: P = transition probability matrix for the inventory system problem. (See Example 5.6)
 
## Inventory Systems Cost Model
Usage: c=ex5invcost(s,S,y,hc,ps,oc)
Input: 

s = base stock level (stock not allowed to go below this level); 

S = Restocking level;

y = row vector of the pmf of weekly demand;

hc = cost of holding one item for one unit of time;

ps = profit from selling one item;

oc = cost of placing an order;

Output: c(i) = expected cost in the current period if the inventory at the beginning is i.

## Manpower Systems
Usage: P=ex5manp(p,l,a).

Input: 

p a row vector, p(i) = probability of promotion from grade i to i+1. 

l = a row vector of the same length as p, l(i) = probability of leaving from grade i.

a = row vector, a(i) = probability that a new employee joins grade i.

Output: P = transition probability matrix for the manpower system problem. (See Example 5.8)

## Manpower Systems Cost Model
Usage: c=ex5manpcost(p,l,a,s,b,d,t)

Input: 
p a row vector, p(i) = probability of promotion from grade i to 

i+1. The last element of p must be zero. 

l = a row vector of the same length as p, l(i) = probability of leaving from grade i.

a = row vector, a(i) = probability that a new employee joins grade i. Must be a valid pmf.

s = a row vector, s(i) = salary of an employee in grade i.

b = a row vector, b(i) = bonus for promotion from grade i to i+1. The last element of b must be zero.

d = a row vector, d(i) = cost of an employee departing from grade i.

t = a row vector, t(i) = cost of training an employee staring in grade i.

Output: c = a row vector, c(i) = expected one-period cost in state i. 

## Manufacturing Systems
Usage: P=ex5mfg(A,B,a1,a2).

Input: 

A = size of the bin for machine 1;

B = size of the bin for machine 2;

a1 = prob(non-defective) for machine 1;

a2 = prob(non-defective) for machine 2.

Output: P = transitrion probability matrix for the manufacturing system. (See Example 5.7.)

## Manufacturing Systems Cost Model
Usage: c=ex5mfgcost(A,B,a1,a2,r,du,hA,hB)

Input:

A = size of bin for machine 1;

B = size of bin for machine 2;

a1 = prob(non-defective) for machine 1;

a2 = prob(non-defective) for machine 2.

r = revenue from a complete assembly;

du = cost of turning a machine on;

hA = cost of holding an item in bin A for one period;

hB = cost of holding an item in bin B for one period;

Output: c = row vector, c(i) = expected one-period cost in state i.
 
## Stock Market
Usage: P=ex5stock(L,U).

Input: 

L = Lower bound for the stock price,

U = Upper bound for the stock price.

Output: P = transition probability matrix for the stock market. (See Example 5.9.)

## Telecommunications.
Usage: P=ex5tel(K,a).

Input: 
K = buffer capacity,

a = row vector, a(i) = p(i-1 packets arrive during one time slot).

Output: P = transition probability matrix for the Telecommunications system. (See Example 5.10.)

## Telecommunications Cost Model
Usage: c=ex5telcost(K,a,rt,cl)

Input: 
K = buffer capacity,

a = row vector, a(i) = p(i-1 packets arrive during one time slot). a must be a valid pmf.

rt = revenue from transmitting a single packet,

cl = cost of losing a single packet.

Output: c = column vector, c(i) = expected cost in one slot if there are i-1 packets in the buffer at the end of the previous slot, 1 <= i = K+1.

## Weather Model.
Usage: y = ex5wea.

Output: The 3X3 matrix of the weather model of Example 5.5. 
