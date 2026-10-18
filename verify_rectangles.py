#!/usr/bin/env python3
"""Reproduce extinction-index crossed-bridge counts and correlations.

Requires numpy and sympy (sympy only to generate small primes for correction).
All prime values are certified by a basic sieve. The FFT coefficients are
rounded, and 8 named cases are cross-checked by an independent direct loop.
"""
import math
import numpy as np
from sympy import primerange

MAX_N = 1_000_000
TMAX = (MAX_N - 1) // 3
limit = 2 * MAX_N + 100
prime = np.ones(limit + 1, dtype=bool)
prime[:2] = False
for d in range(2, math.isqrt(limit) + 1):
    if prime[d]:
        prime[d * d::d] = False

f = prime[6 * np.arange(TMAX + 5, dtype=np.int64) + 1]
f[0] = False
fft_length = 1 << (2 * (len(f) + 4) - 1).bit_length()
Cs = []
for j in (1, 2, 3):
    g = (f[:-j] & f[j:]).astype(np.float64)
    g = np.pad(g, (0, len(f) - len(g)))
    spectrum = np.fft.rfft(g, n=fft_length)
    convolution = np.rint(np.fft.irfft(spectrum ** 2, n=fft_length)).astype(np.int64)
    count = np.zeros(TMAX + 1, dtype=np.int64)
    count[j:] = convolution[: TMAX + 1 - j]
    Cs.append(count)

def direct_count(n, j):
    t = (n - 1) // 3
    total = 0
    for a in range(1, t - j):
        if all(prime[v] for v in (6 * a + 1,
                                  6 * (a + j) + 1,
                                  6 * (t - a) + 1,
                                  6 * (t - a - j) + 1)):
            total += 1
    return total

for nn in (19, 34, 82, 163, 229, 277, 619, 100_000):
    t = (nn - 1) // 3
    for j in (1, 2, 3):
        assert direct_count(nn, j) == Cs[j - 1][t], (nn, j)
print("Selected direct-enumeration checks: PASS")

t = np.arange(1, TMAX + 1)
n = 3 * t + 1
ext = (n >= 7) & (~prime[6 * t - 1])
print("Extinction indices:", int(ext.sum()))
for j, c in enumerate(Cs, 1):
    zz = n[ext & (c[1:] == 0)]
    print(f"C{j} zero ({len(zz)}):", zz.tolist())
for i, j in ((0, 1), (0, 2), (1, 2)):
    zz = n[ext & (Cs[i][1:] == 0) & (Cs[j][1:] == 0)]
    print(f"C{i+1},C{j+1} joint zero ({len(zz)}):", zz.tolist())
allzero = ext & np.logical_and.reduce([c[1:] == 0 for c in Cs])
print("TRIPLE ZEROS:", n[allzero].tolist())
assert int(allzero.sum()) == 0

sel = ext & (n >= 100_000) & (n <= 1_000_000)
ns = n[sel]
ts = t[sel]
y = np.vstack([c[1:][sel] for c in Cs]).astype(float)
print("Correlation sample size:", len(ns))
print("Raw correlation matrix:", np.round(np.corrcoef(y), 6))

local_factors = np.ones_like(y)
for p in primerange(5, 998):
    rem = ns % p
    local_factors[:, rem == 0] *= (p - 2) / (p - 4)
    for j in (1, 2, 3):
        hits = (rem == (3 * j) % p) | (rem == (-3 * j) % p)
        local_factors[j - 1, hits] *= (p - 3) / (p - 4)
normalized = y / (local_factors * (ts / np.log(2 * ns)**4))
print("Normalized correlation matrix:", np.round(np.corrcoef(normalized), 6))
