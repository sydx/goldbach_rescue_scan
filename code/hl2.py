"""Hardy-Littlewood / Bateman-Horn prediction for r_{a,b}(n) with full inclusion-exclusion.
1_{P_a}(p) = sum_{S nonempty subset of {+a,-a}} (-1)^{|S|+1} 1[p+h prime for h in S]  (times 1[p prime])
so r_{a,b}(n) = sum_{S,T} (-1)^{|S|+|T|} #{x : x+h prime (h in {0} u S), n-x+k prime (k in {0} u T)}
and each term ~ C(H1,H2;n) * sum_m (log m)^{-|H1|} (log(n-m))^{-|H2|}.
Evaluated on a set of sample n only (S is per-n; integrals by FFT for all n).
"""
import numpy as np, sys, json, itertools
from scipy.fft import rfft, irfft

N = int(float(sys.argv[1]))
pairs = [tuple(int(x) for x in p.split("-")) for p in sys.argv[2].split(",")]
B = 2000
M = N + 64
s = np.ones(M, bool); s[:2] = False
for i in range(2, int(M**.5) + 1):
    if s[i]: s[i*i::i] = False
primes = np.nonzero(s)[0]
small = primes[primes <= B]

def Pd(d):
    up = s[d:d + N + 1]
    down = np.zeros(N + 1, bool); down[d:] = s[:N + 1 - d]
    return s[:N + 1] & (up | down)

L = 1 << int(np.ceil(np.log2(2 * N + 2)))
m = np.arange(N + 1, dtype=float); lg = np.zeros(N + 1); lg[2:] = 1 / np.log(m[2:])
Fk = {k: rfft(lg**k, L, workers=-1) for k in (1, 2, 3)}
Ik = {}
def I(j, k):
    if (j, k) not in Ik: Ik[(j, k)] = irfft(Fk[j] * Fk[k], L, workers=-1)[:N + 1]
    return Ik[(j, k)]

ev = np.arange(4, N + 1, 2)
sample = np.unique(np.concatenate([ev[ev <= 200000]] +
         [ev[(ev >= lo) & (ev < lo + 20000)] for lo in (10**5, 10**6, N - 20000)]))

def const(H1, H2, n):
    tail = np.prod([(1 - (len(H1) + len(H2)) / l) / (1 - 1/l)**(len(H1) + len(H2)) for l in primes if B < l < 10**6])
    C = np.full(len(n), tail)
    for l in small:
        tab = np.empty(l)
        for r in range(l):
            bad = {(-h) % l for h in H1} | {(r + k) % l for k in H2}
            tab[r] = ((l - len(bad)) / l) / (1 - 1/l)**(len(H1) + len(H2))
        C *= tab[n % l]
    return C

subsets = lambda d: [S for k in (1, 2) for S in itertools.combinations((d, -d), k)]
res = {}
for a, b in pairs:
    Pa, Pb = Pd(a), Pd(b)
    r = np.rint(irfft(rfft(Pa.astype(float), L, workers=-1) * rfft(Pb.astype(float), L, workers=-1), L, workers=-1)[:N + 1])
    act = r[sample].astype(float)
    pred = np.zeros(len(sample)); lead = np.zeros(len(sample))
    for S in subsets(a):
        for T in subsets(b):
            H1, H2 = (0,) + S, (0,) + T
            term = const(H1, H2, sample) * I(len(H1), len(H2))[sample]
            sign = (-1) ** (len(S) + len(T))
            pred += sign * term
            if len(S) == 1 and len(T) == 1: lead += term
    if a == b:
        half = Pa[sample // 2].astype(float)
        act = (act + half) / 2; pred /= 2; lead /= 2
    rec = {}
    for lo in (10**3, 10**4, 10**5, 10**6, N - 20000):
        w = (sample >= lo) & (sample < lo + 20000) if lo >= 10**5 else (sample >= lo) & (sample < 10 * lo)
        q = act[w] / pred[w]; ql = act[w] / lead[w]
        rec[str(lo)] = dict(ratio=float(q.mean()), sd=float(q.std()), ratio_lead=float(ql.mean()),
                            corr=float(np.corrcoef(act[w], pred[w])[0, 1]))
    big = sample >= N - 20000
    Snorm = pred[big] / I(2, 2)[sample[big]]
    rec["S_min_mean_max"] = [float(Snorm.min()), float(Snorm.mean()), float(Snorm.max())]
    p0 = np.exp(-np.clip(pred, 0, None))
    sm = sample <= 200000
    for t in (1000, 10000):
        rec[f"E_exc_gt_{t}"] = float(p0[sm & (sample > t)].sum())
        rec[f"obs_exc_gt_{t}"] = int(((act == 0) & sm & (sample > t)).sum())
    res[f"{a}-{b}"] = rec
    print(a, b, json.dumps(rec), flush=True)
json.dump(res, open(sys.argv[3], "w"))
