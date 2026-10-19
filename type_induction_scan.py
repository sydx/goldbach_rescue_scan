#!/usr/bin/env python3
"""Exact type-count scan (FFT rounded to integers and spot-checked by direct sums).
R_k(r) = sum_{i+j=r-k} 1_P(6i-1) 1_P(6(i+k)+1)
                        1_P(6j+1) 1_P(6(j+k)+5), i,j>=1.
"""
import numpy as np
from scipy.signal import fftconvolve

def sieve(mx):
    a=np.ones(mx+1, dtype=np.bool_); a[:2]=False
    for p in range(2,int(mx**.5)+1):
        if a[p]:a[p*p::p]=False
    return a

R=300000; K=12
prime=sieve(6*(R+K+2)+5)
ix=np.arange(1,R+1)
counts=[]
for k in range(K+1):
    U=np.zeros(R+1,dtype=np.bool_);V=U.copy()
    U[1:]=prime[6*ix-1] & prime[6*(ix+k)+1]
    V[1:]=prime[6*ix+1] & prime[6*(ix+k)+5]
    f=fftconvolve(U.astype(np.float64),V.astype(np.float64))
    x=np.zeros(R+1,dtype=np.int64)
    x[k:]=np.rint(f[:R+1-k]).astype(np.int64)
    residual=np.max(np.abs(f-np.rint(f)))
    assert residual<1e-4, (k, residual)
    counts.append(x)
    # Direct arithmetic confirmation at a few representative r
    for r in [92,93,100,277,1000,9999,R]:
        if 2+k<=r<=R:
            direct=sum(int(U[i] and V[r-k-i]) for i in range(1,r-k))
            assert direct==x[r],(k,r,direct,x[r])
    z=np.flatnonzero(x[2+k:]==0)+2+k
    print(f'k={k:2d} gaps=({6*k+2:2d},{6*k+4:2d}) R_k(92)={x[92]:2d} total zeros r>={k+2}: {z.size:5d} last_zero={z[-1] if z.size else "none":>6}  Rk(100)={x[100]:2d} Rk(1000)={x[1000]:3d}')

c=np.array(counts)
for a,b in [(0,0),(1,1),(1,2),(1,3),(2,4),(3,6),(1,6),(1,12),(0,12)]:
    v=np.sum(c[a:b+1],axis=0)
    begin=max(2,b+2)
    z=np.flatnonzero(v[begin:]==0)+begin
    print(f'group k={a}..{b}: joint zero count={z.size}, last={z[-1] if z.size else "none"}, first={z[:15].tolist()}')

for lo,hi in [(100,1000),(1000,10000),(10000,100000),(100000,R+1)]:
    X=c[:,lo:hi]
    rates=np.mean(X>0,axis=1)
    print('range',lo,hi,'positivity',(rates[:8]).round(5).tolist(),'mean r_0',(X[0].mean()),'r_1',(X[1].mean()))

np.savez_compressed('/mnt/data/type_induction_scan_300k.npz', counts=c)
