# Goldbach midpoint transitions and twin-cousin rescue

Article draft dated 15 October 2026, prepared with *Experimental Mathematics*
as the target journal.

## Files

- goldbach_transitions.tex and references.bib: LaTeX manuscript and bibliography.
- goldbach_transitions.pdf: compiled article.
- verify.mjs and verification.json: independent Node.js checks through one
  million pivots, plus original transition and six-step checks.
- goldbach_rescue_scan.cpp: complete C++17 billion-pivot scanner.
- verification_billion.txt: full output of the four-thread billion-pivot run.

## Build the article

Keep the source and bibliography together and run:

    latexmk -pdf goldbach_transitions.tex

## Reproduce the computations

    node verify.mjs > verification.json
    g++ -O3 -march=native -std=c++17 -pthread goldbach_rescue_scan.cpp -o goldbach_rescue_scan
    ./goldbach_rescue_scan 1000000000 --threads 4 --progress 250000000 --samples 8

The C++ program defaults to an inclusive upper pivot of 1,000,000,000. It
checks both the general conditional transition conjecture and the stronger
residue-specific formulation. A complete fallback search rules out a false
counterexample caused by the initial small-pattern cutoff. Use --help for
other limits and diagnostic options. The bit sieve at one billion occupies
about 119.2 MiB.

The Node.js implementation uses a separate dense sieve and a cousin-first
search; it reproduces the 258,814 successful residue-specific targets
through one million pivots. The C++ run found no failures among 284,220,751
such targets through one billion, and no direct-extinction pivots in the
other two classes.

## Mathematical status

The article proves exact transition and sieve statements and conditional
implications of the rescue conjecture. The billion-range computation is
finite evidence. Neither Goldbach nor the twin-prime conjecture is proved.
The strictly right-hand exchange is disproved already at pivot 277.

# Sexy-prime rescue working note — source package

## Contents

- `sexy_prime_rescue.tex` — LaTeX research note, with precise status of each claim.
- `references.bib` — verified BibTeX citations.
- `verify_convolution.cpp` — exact integer-NTT computation of `R(t)` for `2 <= t <= 1000000`.
- `verify_mid.cpp` — exact direct restricted-convolution test for `t <= 100000` and selected larger inputs.
- `verify_quadruples.cpp` — exact integer-NTT check of the four-prime strengthening for `12 <= n <= 1000000`.
- `verify_local.py` — exact CRT bit-mask certificate for 16 lower sexy primes through sieving prime 23.

## Compiling the paper

```sh
pdflatex sexy_prime_rescue.tex
bibtex8 sexy_prime_rescue  # or bibtex, if available
pdflatex sexy_prime_rescue.tex
pdflatex sexy_prime_rescue.tex
```

## Running computational checks

```sh
g++ -O3 -std=c++17 verify_convolution.cpp -o verify_convolution
./verify_convolution

g++ -O3 -std=c++17 verify_mid.cpp -o verify_mid
./verify_mid

g++ -O3 -std=c++17 verify_quadruples.cpp -o verify_quadruples
./verify_quadruples

python verify_local.py  # requires numpy
```

All finite results in this note are separate from the unproved universal conjectures. The integer-NTT modulus 998244353 exceeds every coefficient relevant to the specified ranges, making the convolutions exact.

Working draft dated 17 October 2026.

# Goldbach Intersection Geometry and Prime-Pair Correlations

Working research manuscript, 18 October 2026.

## Files

- `goldbach_intersections.tex`: main self-contained LaTeX manuscript.
- `references.bib`: BibTeX entries; compiled with biblatex+biber.
- `goldbach_intersections.pdf`: compiled 12-page reading copy.
- `verify_rectangles.py`: independent prime-sieve and FFT computation for three crossed-bridge counts through n = 1,000,000, direct spot-checks, and local-adjusted correlations.
- `computational_results.txt`: reproducible output from the supplied verification script.

## Build

From this directory:

    latexmk -pdf -interaction=nonstopmode goldbach_intersections.tex

Uses LaTeX packages including `biblatex` and the `biber` backend.

To reproduce the calculations:

    python verify_rectangles.py

Python dependencies: numpy, sympy.

## Caution

The paper distinguishes mathematical proofs, finite computations, and Hardy--Littlewood-type heuristics. It does not claim to prove Goldbach or the twin-, cousin-, or sexy-prime conjectures.

# Nearest central Goldbach pairs after complete resets

This directory supports an exact numerical investigation of the first-arrival half-gap
`d(N) = min |p-N/2|` over Goldbach prime pairs `p+q=N`, for reset indices defined by
`N=6t+4`, `t>=2`, and `N-5` composite. In the special guaranteed-reset progression
`N=60s+10`, `s>=1`, the condition is automatic because `N-5=5(12s+1)`.

## Files
- `first_arrival.cpp`: C++17 prime sieve, complete enumeration of `N=60s+10` through 30,000,000, recording nearest pairs and singular-series enhancement factor `H(N)`.
- `reset_first_arrivals_30m.csv.gz`: All 499,999 guaranteed-reset targets between 70 and 29,999,950, columns `N,u,half_gap,j,p,q,H,div7,div11,div13,div17,div19`, with `d(N)=6j`.
- `first_arrival_all_sample.cpp`: Prime sieve across **all** reset indices through 30 million, deterministic sampling by `t % 101 < 20` (~19.8%), recording 813,709 rows. The executable counts but does not check Goldbach pairs on unsampled reset indices.
- `reset_first_arrivals_all_sample_30m.csv.gz`: Broad sample including variable factors of 5, as well as rows with `d=0` (diagonal prime pairs).
- `analyze.py`, `analyze_all.py`: Calculations of aggregate statistics and divisor effects.
- `plot_results.py`: Rebuilds three plots.
- `survival_vs_exponential.png`, `scaling_with_N.png`, `successive_reset_dynamics.png`: Plots.

## Build / run

```
g++ -O3 -std=c++17 first_arrival.cpp -o first_arrival
./first_arrival 30000000 > reset_first_arrivals_30m.csv
python analyze.py
python plot_results.py

g++ -O3 -std=c++17 first_arrival_all_sample.cpp -o first_arrival_all_sample
./first_arrival_all_sample 30000000 | gzip -1 > reset_first_arrivals_all_sample_30m.csv.gz
python analyze_all.py
```

`analyze.py` currently reads the uncompressed guaranteed-reset CSV by default; if the archive has only the `.csv.gz`, unpack it first, or edit the file path to `.csv.gz` (pandas supports both).

## Analytic model

The classical Hardy–Littlewood enhancement factor for even N is

`H(N) = product_{ell|N, ell>2} (ell-1)/(ell-2)`.

With twin-prime constant `C2 = 0.6601618158468696`, the **heuristic** local arrival
intensity per unit positive half-gap `d` is

`lambda_N ≈ 2*C2*H(N) / log(N/2)^2`.

The model predicts, **heuristically**,

`Z(N) = 2*C2*H(N)*d(N) / log(N/2)^2`

has a roughly exponential survival curve. This is not an established theorem,
and local congruence conditions and correlations prevent exact independence.

For `N=60s+10`, the midpoint is `m=30s+5`, so the nonexceptional candidate pair is `m-6j, m+6j`, with `j>=1`. The sieve excludes every `j ≡ 0 (mod 5)` exactly. For any odd prime `ell | N`, `j ≡ 0 (mod ell)` is likewise excluded.

All compiled figures and numerical conclusions are empirical and confined to the stated finite ranges.

# Joint local sieve and Goldbach first-arrival dynamics

This analysis studies the guaranteed-reset progression N=60s+10, N between 1,000,000 and 30,000,000 inclusive. The midpoint is m=N/2, and the candidate event A_j(N) is that both m-6j and m+6j are prime. Candidate indices j divisible by 5 are identically impossible in this progression.

## Reproduction

```
g++ -O3 -std=c++17 joint_sieve_test.cpp -o joint_sieve_test
./joint_sieve_test > joint_sieve_results.txt
python local_sieve_prediction.py > local_sieve_predictions.txt

g++ -O3 -std=c++17 factorial_cumulants.cpp -o factorial_cumulants
./factorial_cumulants > factorial_cumulant_results.txt

g++ -O3 -std=c++17 window_variance.cpp -o window_variance
./window_variance
```

Python singular-series prediction script uses `sympy` to obtain primes through 10000. The Euler products are truncated: they are heuristic predictions for correlations of actual primes, not theorems about prime tuple frequencies. The local residue-counting identities are exact. The sample tests all 483,333 targets N in the stated progression and range. The C++ computation uses an Eratosthenes sieve through 30,000,000.

`joint_sieve_results.txt`: `lag, eligible_trials, both_hits, rate_joint, p_independent, ratio` and first-arrival-window moments for the first 120 eligible/noneligible positions.

`local_sieve_predictions.txt`: singular-series marginal weights and predicted correlations for lags 1 through 26.

`factorial_cumulant_results.txt`: rows are width, mean count, second factorial cumulant, third factorial cumulant, observed void frequency, naive Poisson void frequency, second-cumulant approximation, third-cumulant approximation.

The finite-sieve product for a set J of candidate indices and prime ell>=7 is `1 - |{+j,-j (mod ell): j in J}|/ell`, when N ranges over the specified progression. Do not interpret finite-sieve survival as certification that candidates are prime.

# Prime-gap types, crossed rescues, and multitype induction

Working mathematical manuscript, 19 October 2026.

## Files
- `multitype_goldbach.tex`: complete LaTeX source
- `references.bib`: BibTeX bibliography
- `multitype_goldbach.pdf`: compiled 10-page PDF
- `multitype_goldbach.bbl`: resolved bibliography for compiling without BibTeX
- `figures/type_zero_heatmap.pdf` and `.png`: figure for the manuscript
- `code/type_induction_scan.py`: FFT-based 13-type scan through r=300,000
- `code/type_induction_scan_300k.npz`: archived numerical counts, shape (13,300001)
- `code/verify_exact_examples.py`: direct integer primality and CRT verification for r=92, and exact rational correction factors.

## Compile
From the package root:

    pdflatex -interaction=nonstopmode multitype_goldbach.tex
    bibtex8 multitype_goldbach
    pdflatex -interaction=nonstopmode multitype_goldbach.tex
    pdflatex -interaction=nonstopmode multitype_goldbach.tex

Standard `bibtex` works instead of `bibtex8` where available. If no BibTeX implementation is installed, the provided `.bbl` file lets pdflatex use the already resolved bibliography.

## Reproduce numerical checks
Requires Python 3, numpy, scipy, sympy and matplotlib (matplotlib only for figure regeneration).

    python code/verify_exact_examples.py
    python code/type_induction_scan.py

The FFT-based count scan rounds near-integer convolution results and spot-checks selected coefficients by integer summation. It is finite computational evidence, not a formal proof of all scanned coefficients. The direct verifier uses exact integer arithmetic.

## Mathematical status
The transition identities, local singular-series factors, cross-type congruence exclusions, and finite-sieve theorem are presented with proofs. The generalized fixed-gap Goldbach conjecture, four-prime asymptotic, and bounded-type renewal conjecture remain unproved. In particular, the manuscript does **not** prove Goldbach's conjecture.

# Twin--Cousin Goldbach Representations

Paper author: Paul A. Bilokon
Date: 21 October 2026

Files:
- `twin_cousin_goldbach.tex`: complete, editable LaTeX source.
- `twin_cousin_goldbach.bib`: BibTeX bibliography.
- `twin_cousin_goldbach.pdf`: compiled eight-page paper.
- `verify_twin_cousin.py`: standalone Python 3 script (standard library only) that reproduces the finite calculations.
- `verification_output.txt`: output from the default verification run.

To compile with a standard TeX installation from this directory:

```sh
pdflatex twin_cousin_goldbach.tex
bibtex twin_cousin_goldbach
pdflatex twin_cousin_goldbach.tex
pdflatex twin_cousin_goldbach.tex
```

To recompute the finite evidence:

```sh
python3 verify_twin_cousin.py
```

The manuscript sharply distinguishes exact theorems and finite numerical verification from the still-open Twin--Cousin Goldbach (TCG) and Twin--Cousin Block (TCB) conjectures. In particular, local admissibility alone is not a prime-existence proof.
