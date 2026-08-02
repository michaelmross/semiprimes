# One-step-Fermat semiprimes in J_n

An exhaustive census, for 2 ≤ n ≤ 10⁸, of the one-step-Fermat semiprimes in the
square-centred intervals

```
J_n = [4n² − n, 4n² + n]
```

Call a semiprime N = pq **one-step** if Fermat's method splits it on the first
trial: with m = ⌈√N⌉, the quantity m² − N is a perfect square. One-step
semiprimes are simultaneously the semiprimes most resistant to trial division
and the ones Fermat's method breaks instantly.

## The reduction

One-step semiprimes in J_n are in bijection with Goldbach representations
having **almost equal summands**:

| half of J_n | representation | window on the summand |
|---|---|---|
| lower, [4n² − n, 4n²] | 4n = p + q | \|p − 2n\| ≤ √n |
| upper, (4n², 4n² + n] | 4n + 2 = p + q | √(3n) < \|p − (2n+1)\| ≤ √(4n) |

Writing the prime gap as q − p = 2x, the two halves are counted by

```
R(n) = #{x odd,  1 ≤ x ≤ √n        : 2n−x,   2n+x   both prime}
U(n) = #{x even, 3n+1 ≤ x² ≤ 4n    : 2n+1−x, 2n+1+x both prime}
```

and R(n) + U(n) is the number of one-step semiprimes in J_n. The gap lies in
the class 2 mod 4 in the lower half and 0 mod 4 in the upper half — the entire
congruence structure lives on the gap. The window exponent is θ = 1/2 with no
logarithmic slack, which places the pointwise question below every result in
the almost-equal-summands literature.

The parity condition on x is not an extra hypothesis: if x has the wrong
parity, both factors are even and N is divisible by 4, so N is not a semiprime
above 4. Nor does x = 0 contribute, since J_n contains no prime square for
n ≥ 2.

## Results

| list | meaning | zeros | first | last | beyond 10⁶ | beyond 10⁷ |
|---|---|---|---|---|---|---|
| `zeros_R.csv` | lower half of J_n empty | 8,660 | 4 | 4,409,947 | 154 | 0 |
| `zeros_U.csv` | upper half of J_n empty | 444,749 | 2 | 99,919,133 | 260,111 | 49,942 |
| `zeros_J.csv` | all of J_n empty (R = U = 0) | 3,226 | 5 | 1,884,296 | 9 | 0 |

The headline statistic: **n = 1,884,296 is the last n ≤ 10⁸ whose interval J_n
contains no one-step-Fermat semiprime.** The lower half alone is last empty at
n = 4,409,947. The upper half, whose admissible window is shorter by a factor
of 2 − √3 ≈ 0.268, is still empty for 49,942 values of n above 10⁷ and gives no
sign of a last zero within the census range.

Late zeros are arithmetically structured exactly as the singular series
predicts. Of the 154 R-zeros beyond 10⁶, 43% have ω(n) ≤ 1 counting distinct
odd prime divisors, and their mean 𝔖(4n) is 1.358 against a global mean of
2.001: the hard indices are precisely those with small singular series.

## Contents

```
census_final.py     the census: counts, Hardy–Littlewood band ratios, Poisson
                    zero predictions, singular-series stratification, CSV dumps
zeros_R.csv         n with R(n) = 0,        one column, header "n"
zeros_U.csv         n with U(n) = 0,        one column, header "n"
zeros_J.csv         n with R(n) + U(n) = 0, one column, header "n"
verify_zeros.py     independent verifier for the three CSVs (see below)
verify_reduction.py exact identity check: the census count equals the
                    prime-pair (shift-form) count, integer for integer
CHECKSUMS.sha256    SHA-256 of every file in this directory
```

All three CSVs are strictly increasing, LF-terminated, and ASCII.
`zeros_J.csv` is exactly the intersection of the other two.

## Reproduce

```bash
pip install numpy

python census_final.py 1e6     # minutes
python census_final.py 1e8     # the published census; memory-bound, see below
```

The script self-checks on every run: it compares its counts against a direct
interval scan (factor every N in J_n, test Fermat's first trial) for n ≤ 300,
and against the known value Σ_{n≤10⁵} R(n) = 319,013.

At X = 10⁸ the run needs roughly 400 MB for the prime sieve, 800 MB for the
smallest-prime-factor sieve and 400 MB for the count arrays — about 1.6 GB, so
a 4 GB machine is enough. It is single-threaded and memory-bound.

## Verify

`verify_zeros.py` checks the artifacts without trusting `census_final.py`:

```bash
python verify_zeros.py               # checks 1-4, completeness to N = 10⁶
python verify_zeros.py 1e7           # completeness to N = 10⁷
python verify_zeros.py 1e7 --bpsw 200  # adds check 5
```

It runs four independent checks, plus an optional fifth:

1. **checksums** — every file present matches `CHECKSUMS.sha256`.
2. **structure** — headers, strict monotonicity, and `zeros_J = zeros_R ∩ zeros_U`.
3. **soundness** — every one of the 456,635 listed zeros is confirmed to be a
   zero, by enumerating the admissible x for that n directly rather than by the
   sliding-slice method the census uses, so a slicing error cannot hide in both.
4. **completeness** — the zero sets are rebuilt for n ≤ N by an independently
   written recount and must match the CSVs on that range, i.e. no zero is
   missing.
5. **bpsw** (`--bpsw K`) — each record plus K random indices per list are
   re-tested with sympy's Baillie–PSW primality test, a primality path sharing
   nothing with the sieve.

Checks 1–3 always cover the full lists. Check 4 covers n ≤ N and costs
O(N^1.5); at N = 10⁷ it certifies `zeros_R.csv` and `zeros_J.csv` end to end,
since both lists terminate below 10⁷. Completeness of the U tail on (10⁷, 10⁸]
is certified only by the census run itself. The sieve for checks 3–4 is sized
by the largest listed zero (≈ 2 × 10⁸) regardless of N, so the run needs about
300 MB; expect a few minutes at N = 10⁷, machine-dependent.

`verify_reduction.py` checks the mathematics rather than the artifacts: at
X = 10⁵ it computes Σ R(n) both as the census does and as a count of prime
pairs indexed by the gap (the exact reduction of the paper), and requires
integer-for-integer agreement (both sides equal 319,013); it also confirms the
two singular-series means (→ 2) and the Hardy–Littlewood band ratios. No
arguments; runs in seconds.

Cross-platform note: `zeros_J.csv` has been regenerated from scratch on an
independent Linux machine at X = 2 × 10⁶ (the list is complete below that
bound) and is byte-identical, SHA-256
`5638e60447998006…`, to the file produced by the X = 10⁸ census on Windows.

## Citation

If you use this census, please cite the repository archive on Zenodo — see the
badge in the top-level `README.md` — and the accompanying paper.

## License

Code MIT, data CC BY 4.0, as for the rest of the repository.
