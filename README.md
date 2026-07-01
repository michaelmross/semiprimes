# Exact semiprimes in short intervals — computation

Code and data for the computational section of

> M. M. Ross, *Exact Semiprimes in Short Intervals: Prime Supply, Factor-Pair
> Multiplicity, and the Parity Barrier* (2026).
> https://doi.org/110.5281/zenodo.21089081

An exact enumeration of Ω(n) for all n < 10¹¹, producing the global-density,
worst-case-gap, square-interval, dispersion, and gap-self-similarity statistics
reported in the paper. The computation is fully deterministic and the figures and
tables are reproducible end to end.

## Contents

```
code/
  semiprime_production.py   main segmented Ω-sieve + all six experiments
  export_npz.py             write results.npz from a checkpoint (no re-run)
  finalize.py               re-emit the full report from a checkpoint
  plot_pub.py               produce the three paper figures + summary.csv
data/
  semi_results.npz          archived results to 10¹¹ (arrays for every table/figure)
  summary.csv               per-decade table (count, Landau term, ratio, gap, …)
figures/
  fig1_gap_vs_sqrt.png      G_E2(x) vs √x and (log x)²
  fig2_square_detrend.png   detrended E2(n)/Ê2(n), odd vs even
  fig3_gap_selfsim.png      normalized-gap distribution and CV per decade
```

## Method

A segmented sieve computes Ω(n) on each block [A, B) by dividing out every prime
p ≤ √B and recording the single residual large prime, so primes (Ω=1), exact
semiprimes (Ω=2), and 3-almost primes (Ω=3) come from one pass. Memory is O(segment),
so the height is bounded only by patience. The run checkpoints to a pickle every few
segments and resumes automatically if interrupted.

## Reproduce

```bash
pip install -r requirements.txt

# full run to 10^11 (overnight, single core; see runtime/RAM below)
python code/semiprime_production.py 1e11

# if the run is interrupted, or to read results without re-running:
python code/export_npz.py            # writes semi_results.npz from semi_ckpt.pkl

# regenerate the three figures + summary.csv from the archived results:
python code/plot_pub.py
```

To skip the long run entirely, point `plot_pub.py` at the archived
`data/semi_results.npz` (copy it next to the script) and every figure and table
regenerates in seconds.

A quick sanity run finishes in under a minute and prints the validation block:

```bash
python code/semiprime_production.py 1e8
```

## Runtime and memory

Single-threaded; ~3 million n/s on one modern core.

| height | wall time (1 core) | peak RAM (default segment) |
|--------|--------------------|----------------------------|
| 10¹⁰   | ~1 hour            | ~0.5 GB                    |
| 10¹¹   | ~8–16 hours        | ~0.5 GB                    |

Lower `SEG` (or pass it as the second argument) to reduce RAM. The work is
per-segment independent, so 10¹² is best parallelized rather than run serially.

## Determinism and validation

The sieve is exact and deterministic. The only randomness is the random-window
experiment (`[F]`), which uses a fixed seed (`np.random.default_rng(0)`), so the
full output is bit-reproducible. Each run self-validates:

- counts checked against the known values π_{E₂}(10ᵏ) for k ≤ 8;
- a second, independent least-prime-factor lane sum
  Σ_{p≤√X}(π(X/p) − π(p) + 1) agrees to the digit at 10⁶, 10⁷, 10⁸;
- a direct `sympy` factorization spot check confirms Ω near the top of the range.

Reference: π_{E₂}(x) (exact semiprimes, p² included) below 10¹¹ is
**13,959,990,342**.

## Using the archived results

```python
import numpy as np
d = np.load("data/semi_results.npz")
d["E2arr"][n]        # semiprime count in (n², (n+1)²], valid for n ≤ d["nmax_valid"]
d["Gcum"][k]         # worst-case gap G_E2(10^k)
d["gap_hist"][k]     # raw gap histogram for decade [10^k, 10^{k+1})
d["dec_semi"], d["dec_tri"], d["dec_prime"]   # per-decade Ω=2, Ω=3, Ω=1 counts
```

## Notes

- The progress line prints an ETA; after a resume it shows a transient hump
  (the new process's elapsed clock restarts while the percent-complete does not)
  before settling into a countdown. Cosmetic only — the results are unaffected.

## License

Code is released under the MIT License (see `LICENSE`). The data files and figures
are released under CC BY 4.0.
