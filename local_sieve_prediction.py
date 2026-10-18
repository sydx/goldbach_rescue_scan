from math import log,exp
from sympy import primerange
J=120
L=26
primes=list(primerange(7,10001))
log_weights={j:0 for j in range(1,J+1) if j%5}
for j in log_weights:
 for ell in primes:
  nu=1 if j%ell==0 else 2
  log_weights[j]+=log(1-nu/ell)-2*log(1-1/ell)
w={j:exp(log_weights[j]) for j in log_weights}
print('j marginal weights',[(j,round(w[j],4)) for j in (1,2,3,4,6,7,8,11,14,21)])
for lag in range(1,L+1):
 pairs=[(j,j+lag) for j in log_weights if j+lag in log_weights]
 numer=0; denom=0
 for j,h in pairs:
  logratio=0
  for ell in primes:
   nu_j=1 if j%ell==0 else 2
   nu_h=1 if h%ell==0 else 2
   nu_jh=len({j%ell,(-j)%ell,h%ell,(-h)%ell})
   logratio+=log(1-nu_jh/ell)-log(1-nu_j/ell)-log(1-nu_h/ell)
  numer+=w[j]*w[h]*exp(logratio)
  denom+=w[j]*w[h]
 print(lag, round(numer/denom,6))
