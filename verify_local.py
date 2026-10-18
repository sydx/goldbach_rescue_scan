"""Exact CRT verification of the first 16 sexy-rescue candidates.
Requires numpy. All operations use integer modular arithmetic.
"""
import numpy as np
from math import prod

Q = [5, 11, 17, 23, 41, 47, 53, 83,
     101, 107, 131, 167, 173, 191, 227, 233]
primes = [5, 7, 11, 13, 17, 19]
P19 = prod(primes)
t = np.arange(P19, dtype=np.int32)
mask = np.zeros(P19, dtype=np.uint16)
for ell in primes:
    lut = np.zeros(ell, dtype=np.uint16)
    for x in range(ell):
        root = (6*x - 2) % ell
        lut[x] = sum(1 << i for i, q in enumerate(Q) if q % ell == root)
    mask |= lut[t % ell]

popcount = np.fromiter((i.bit_count() for i in range(65536)),
                       dtype=np.uint8, count=65536)
survivors = 16 - popcount[mask.astype(np.int32)]
assert survivors.min() == 1
assert np.count_nonzero(survivors == 1) == 2
print('P19', P19, 'min survivors', survivors.min())
print('singleton residue classes:', np.where(survivors == 1)[0].tolist())

empty_at_23 = []
for t0 in np.where(survivors == 1)[0]:
    survivor = [q for i, q in enumerate(Q) if not (int(mask[t0]) & (1 << i))]
    assert len(survivor) == 1
    for k in range(23):
        lifted = int(t0 + k*P19)
        if (6*lifted - 2 - survivor[0]) % 23 == 0:
            empty_at_23.append(lifted)
assert sorted(empty_at_23) == [24004568, 34450388]
print('P23:', 23*P19, 'empty residue classes:', sorted(empty_at_23))
