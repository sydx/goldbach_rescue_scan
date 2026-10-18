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
