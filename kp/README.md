# kp_semiprimes.py

Code for [*Semiprimes 2p and 3p in Quadratic Intervals: A Legendre-Class Persistence Problem*](https://doi.org/10.5281/zenodo.19682906) (Ross, M. M., 2026).

Compute and verify the existence of semiprimes of the form *kp* (with *p* prime) in quadratic intervals *(n², (n+1)²)*.

## Background

For a prime *k*, define *S_k(n)* as the count of semiprimes *kp* in the interval *(n², (n+1)²)*, where *p* is prime. Equivalently, *S_k(n)* counts primes in *(n²/k, (n+1)²/k)*.

This script verifies that *S_2(n) ≥ 1* and *S_3(n) ≥ 1* for all *n* up to a given bound (existence mode), or computes exact counts *S_k(n)* for each interval (count mode).

## Requirements

- Python 3.8+
- NumPy (optional, improves sieve performance in count mode)

## Usage

### Existence mode (default)

Checks whether *S_k(n) ≥ 1* for each *n*, using deterministic Miller-Rabin. Fast enough for *n* up to 10⁷ or beyond.

```bash
# Verify k=2,3 to 10 million (the paper's main result)
python kp_semiprimes.py --nmax 10000000

# Verify k=2,3,5,7,11 to 100,000
python kp_semiprimes.py --nmax 100000 --kvals 2,3,5,7,11

# Adjust progress reporting interval
python kp_semiprimes.py --nmax 10000000 --progress 500000
```

### Count mode

Computes exact *S_k(n)* for every interval via segmented sieve. Slower; practical to *n* ≈ 100,000–500,000 depending on hardware.

```bash
# Full counts for k=2,3,5,7,11 to 100,000 with CSV output
python kp_semiprimes.py --nmax 100000 --kvals 2,3,5,7,11 --output counts_100k.csv

# Counts for k=2,3 only
python kp_semiprimes.py --nmax 200000 --output counts_200k.csv
```

### Resuming or splitting runs

Use `--nmin` to start from a specific *n*. Useful for splitting large runs across multiple processes.

```bash
python kp_semiprimes.py --nmin 5000001 --nmax 10000000
```

## Options

| Flag | Description | Default |
|------|-------------|---------|
| `--nmax N` | Upper bound for *n* | 1,000,000 |
| `--nmin N` | Lower bound for *n* | 2 |
| `--kvals K1,K2,...` | Comma-separated list of *k* values | `2,3` |
| `--count` | Use count mode (segmented sieve) instead of existence mode | off (existence) |
| `--progress N` | Report progress every *N* intervals | 100,000 |
| `--output FILE` | Write per-interval CSV (count mode only) | none |

## Output

### Terminal

Both modes print periodic progress reports and a final summary showing failures (if any), mean and minimum counts (count mode), and block minimums.

### CSV (count mode with `--output`)

One row per interval:

```
n,S_2,S_3,S_5,S_7,S_11
2,1,1,0,0,0
3,2,1,2,1,0
4,1,1,0,1,1
...
```

## Primality testing

Existence mode uses a deterministic Miller-Rabin test with 12 witnesses (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37), which is provably correct for all *n* < 3.317 × 10²⁴. Results are not probabilistic.

## Performance

Benchmarks on a single core (results will vary):

| Mode | Range | k values | Time |
|------|-------|----------|------|
| Existence | *n* ≤ 10⁷ | 2, 3 | ~2 hours |
| Count | *n* ≤ 10⁵ | 2, 3, 5, 7, 11 | ~24 minutes |

Existence mode runs at roughly 1,000–5,000 *n*/s depending on *n* (larger *n* is slightly slower per step but finds the first prime quickly due to high density). Count mode is slower because it sieves the entire sub-interval.

## Results summary

As reported in the paper:

- **k = 2:** *S_2(n) ≥ 1* for all *n* ≤ 10⁷ (zero failures). Intervals up to 10¹⁴.
- **k = 3:** *S_3(n) ≥ 1* for all *n* ≤ 10⁷ (zero failures).
- **k = 5:** 5 failures, last at *n* = 24.
- **k = 7:** 6 failures, last at *n* = 25.
- **k = 11:** 15 failures, last at *n* = 121.

## License

[Repository license file](https://github.com/michaelmross/semiprimes/blob/main/LICENSE.md).
