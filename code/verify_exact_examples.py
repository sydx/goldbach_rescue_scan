#!/usr/bin/env python3
"""Independent integer arithmetic: exact r=92 configurations, CRT claims,
finite corrections to the general four-linear-forms singular series."""
from fractions import Fraction
from math import isqrt
from sympy import factorint


def prime_flags(n):
    f = bytearray(b'\x01')*(n+1)
    f[0:2] = b'\x00\x00'
    for p in range(2,isqrt(n)+1):
        if f[p]: f[p*p::p] = b'\x00'*(((n-p*p)//p)+1)
    return f

P=prime_flags(2000)

def pairs(N):
    return [(p,N-p) for p in range(3,N//2+1,2) if P[p] and P[N-p]]

r=92
v=[]
for k in range(r-1):
    count=0
    for i in range(1,r-k):
        j=r-k-i
        if j>=1 and all(P[z] for z in (6*i-1,6*(i+k)+1,6*j+1,6*(j+k)+5)):
            count+=1
    v.append(count)
assert sum(v)==231
assert sum(c>0 for c in v)==72
assert max(k for k,c in enumerate(v) if c)>0
assert [v[k] for k in range(8)]==[1,5,8,6,1,2,7,7]
assert len(pairs(554))==11 and len(pairs(556))==11
assert set(sum(([a,b] for a,b in pairs(554)),[])).isdisjoint(set(sum(([a,b] for a,b in pairs(556)),[])))
print('r=92: total types',sum(c>0 for c in v),'total oriented rescues',sum(v),'last nonempty',max(k for k,c in enumerate(v) if c))
print('Pairs for 554:',pairs(554))
print('Pairs for 556:',pairs(556))
print('First 13 type counts:',v[:13])


def coefficient_C0(a,b,N):
    assert N>a+b
    coeff=Fraction(0,1)
    for eps in (-1,1):
        for delta in (-1,1):
            residues={0,-eps*a,N,N+delta*b}
            f=Fraction(3-len({x%3 for x in residues}),1)
            D=abs(a*b*N*(N+eps*a)*(N+delta*b)*(N+eps*a+delta*b))
            for ell in factorint(D):
                if ell>=5:
                    f *= Fraction(ell-len({x%ell for x in residues}), ell-4)
            coeff+=f
    return coeff*(Fraction(27,4) if a==b else Fraction(27,2))

for a,b,N in ((2,4,104),(6,6,122),(10,10,104),(2,30,104)):
    print((a,b,N),coefficient_C0(a,b,N))
assert coefficient_C0(2,4,104)==Fraction(1500,91)
assert coefficient_C0(6,6,122)==Fraction(18271,475)
for r in range(5,500,5):
    F1={(x%5) for x in (1,-7,6*(r-1)+1,6*r+5)}
    F2={(x%5) for x in (1,-13,6*(r-2)+1,6*r+5)}
    assert F1|F2 == set(range(5))
print('CRT type-1/type-2 mod-5 union verified')
