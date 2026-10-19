#!/usr/bin/env python3
"""Verify twin–cousin four-sum tiles up to 6*R (default R=1666666).

A[a] = 1 iff (6a-1,6a+1) are both prime.
B[b] = 1 iff (6b+1,6b+5) are both prime.
H[r] = sum_{a=1}^{r-1} A[a] B[r-a].
The FFT yields counts; an independent arbitrary-precision bitset OR verifies
positive support exactly, without floating-point arithmetic.
"""
import argparse
import numpy as np
from scipy.signal import fftconvolve


def sieve(limit):
    isprime = np.ones(limit + 1, dtype=np.bool_)
    isprime[:2] = False
    for p in range(2, int(limit**0.5)+1):
        if isprime[p]:
            isprime[p*p::p] = False
    return isprime


def main(R):
    p = sieve(6*(R+1)+5)
    ix = np.arange(1,R+1)
    A = np.r_[False, p[6*ix-1] & p[6*ix+1]]
    B = np.r_[False, p[6*ix+1] & p[6*ix+5]]
    floatcounts=fftconvolve(A.astype(np.float64),B.astype(np.float64))[:R+1]
    H=np.rint(floatcounts).astype(np.int64)
    err=float(np.max(abs(floatcounts-np.rint(floatcounts))))
    assert err<1e-3, f'FFT rounding error too large: {err}'
    # Exact coverage certificate using arbitrary precision integers:
    bits=int.from_bytes(np.packbits(B,bitorder='little').tobytes(),'little')
    covered=0
    for a in np.flatnonzero(A):
        covered |= bits<<int(a)
    mask=(1<<(R-1))-1
    assert ((covered>>2)&mask).bit_count()==int(np.count_nonzero(H[2:])), 'Count/support disagreement'
    missing=np.flatnonzero(H[2:]==0)+2
    print('Max r:',R,'max even 6r:',6*R)
    print('Twin-index count:',int(A.sum()),'cousin-index count:',int(B.sum()))
    print('Max FFT residual from integer:',err)
    print('Exactly certified missing r:',len(missing),'first:',missing[:20].tolist())
    print('Exactly certified covered r:',((covered>>2)&mask).bit_count())
    print('Indices r with H(r)=1:',np.flatnonzero(H==1).tolist())
    for r in [2,3,4,5,92,93,1000,10000,100000,R]:
        if r<=R:
            print('H(%d) = %d' %(r,H[r]))
            if r==92:
                pairs=[(6*a-1,6*(r-a)+1) for a in range(1,r) if A[a] and B[r-a]]
                print('r=92 lower twin / lower cousin witnesses:',pairs)
    if R>=1000000:
        L=500000;U=min(R+1,1600000)
        for mod in range(5):
            v=np.arange(L+(mod-L)%5,U,5)
            score=float(np.mean(H[v]/(v/np.log(6*v)**4)))
            print('r mod5=%d: mean H/(r/log^4(6r))=%.3f'%(mod,score))

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--R',type=int,default=1666666)
    args=ap.parse_args()
    main(args.R)
