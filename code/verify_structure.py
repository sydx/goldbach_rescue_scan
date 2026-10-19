"""Check the mod-3 decomposition (Theorem 4.2) and Corollary 4.3 numerically.
For even a, b with 3 not dividing ab:
  L_a = {p > 3 : p, p+a prime},  U_a = L_a + a,  E_a = P_a minus (L_a u U_a)
  r*(n) = #{(p,q) in (L_a u U_a) x (L_b u U_b) : p+q = n}  (ordered)
  N(m)  = #{(x,y) in L_a x L_b : x+y = m}
Claims: E_a subset {3, a+3};  r*(n) = N(n)+N(n-a)+N(n-b)+N(n-a-b);  0 <= r(n)-r*(n) <= |E_a|+|E_b|.
"""
import numpy as np, sys
from numpy.fft import rfft, irfft
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**6
s = np.ones(N + 100, bool); s[:2] = False
for i in range(2, int((N + 100)**.5) + 1):
    if s[i]: s[i*i::i] = False
def Pd(d):
    up = s[d:d+N+1]; dn = np.zeros(N+1, bool); dn[d:] = s[:N+1-d]; return s[:N+1] & (up | dn)
def Ld(d):
    x = s[:N+1] & s[d:d+N+1]; x[:4] = False; return x
def shift(x, d):
    y = np.zeros_like(x); y[d:] = x[:len(x)-d]; return y
L_ = 1 << 22
conv = lambda x, y: np.rint(irfft(rfft(x.astype(float), L_) * rfft(y.astype(float), L_), L_)[:N+1]).astype(np.int64)
gaps = [2, 4, 8, 10, 14, 16, 20, 22, 26, 28]
allok = True
for i, a in enumerate(gaps):
    for b in gaps[i:]:
        Pa, Pb, La, Lb = Pd(a), Pd(b), Ld(a), Ld(b)
        Ua, Ub = shift(La, a), shift(Lb, b)
        Ea, Eb = np.nonzero(Pa & ~(La | Ua))[0], np.nonzero(Pb & ~(Lb | Ub))[0]
        ok1 = set(Ea) <= {3, a+3} and set(Eb) <= {3, b+3} and not (La & Ua).any() and not (Lb & Ub).any()
        rstar = conv(La | Ua, Lb | Ub); Nm = conv(La, Lb); r = conv(Pa, Pb)
        n = np.arange(a + b, N + 1)
        ok2 = np.array_equal(rstar[n], Nm[n] + Nm[n-a] + Nm[n-b] + Nm[n-a-b])
        d = r - rstar; ok3 = d.min() >= 0 and d.max() <= len(Ea) + len(Eb)
        allok &= ok1 and ok2 and ok3
        if not (ok1 and ok2 and ok3): print("FAIL", a, b, ok1, ok2, ok3)
print("all checks passed" if allok else "some checks failed", "for gaps", gaps, "up to", N)
