"""Exact (integer-only) determination of the exceptional set
E_{a,b}(N) = {even n in [4,N] : n is not p+q with p in P_a, q in P_b}
and of the largest minimal witness max_n min{p in P_a : n-p in P_b}.
No floating point is used: every non-exception gets an explicit witness.
"""
import numpy as np, sys, json, time

N = int(float(sys.argv[1]))
pairs = [tuple(int(x) for x in p.split("-")) for p in sys.argv[2].split(",")]
out_path = sys.argv[3]

maxgap = max(max(p) for p in pairs)
M = N + maxgap + 2
s = np.ones(M, bool); s[:2] = False
for i in range(2, int(M**.5) + 1):
    if s[i]: s[i*i::i] = False

def Pd(d):
    up = s[d:d + N + 1]
    down = np.zeros(N + 1, bool); down[d:] = s[:N + 1 - d]
    return s[:N + 1] & (up | down)

res = {}
for a, b in pairs:
    t0 = time.time()
    Pa, Pb = Pd(a), Pd(b)
    A = np.nonzero(Pa)[0]
    U = np.arange(4, N + 1, 2, dtype=np.int64)   # still-uncovered even n
    exceptions = []
    maxwit = (0, 0)      # (witness p, n)
    for p in A:
        # n < p + min(P_b) can never be covered by p or any later element
        qmin = 3 if Pb[3] else int(np.argmax(Pb))
        done = U < p + qmin
        if done.any():
            exceptions.extend(int(x) for x in U[done]); U = U[~done]
        if U.size == 0: break
        hit = Pb[U - p]
        if hit.any():
            maxwit = (int(p), int(U[hit][-1]))
            U = U[~hit]
    exceptions.extend(int(x) for x in U)   # anything left (only if P_a exhausted)
    exceptions.sort()
    res[f"{a}-{b}"] = {"N": N, "exceptions": exceptions,
                       "threshold": (max(exceptions) + 2) if exceptions else 4,
                       "max_min_witness": maxwit[0], "at_n": maxwit[1]}
    print(f"{a}-{b}: {len(exceptions)} exc, threshold {res[f'{a}-{b}']['threshold']}, "
          f"max min-witness p={maxwit[0]} (n={maxwit[1]})  [{time.time()-t0:.0f}s]", flush=True)
json.dump(res, open(out_path, "w"))
