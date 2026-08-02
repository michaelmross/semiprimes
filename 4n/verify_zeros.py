"""
verify_zeros.py -- independent re-verification of the deposited zero lists.

Usage: python verify_zeros.py [zeros_J.csv] [zeros_R.csv]

Independence: uses sympy.isprime (deterministic BPSW/MR path) directly on the
candidates 2n +/- x and 2n+1 +/- x; shares NO code with the sieve pipeline
(census_final.py) that produced the lists.

Checks:
  1. Every n in zeros_J.csv is genuinely one-step-empty: no admissible x in
     either half of J_n yields a prime pair.
  2. Every n in zeros_R.csv has empty lower half.
  3. Spot check: the 50 indices immediately above the last entry of each list
     are NOT in the empty condition (i.e., the record is a boundary, not an
     artifact of truncation), except where they legitimately appear in the
     other list.
  4. Reports SHA256 of each input file for cross-checking against CHECKSUMS.
"""
import sys, csv, hashlib
from math import isqrt
from sympy import isprime

def lower_empty(n):
    m = 2 * n
    for x in range(1, isqrt(n) + 1, 2):          # only odd x can give prime pairs
        if isprime(m - x) and isprime(m + x):
            return False
    return True

def upper_empty(n):
    m = 2 * n + 1
    lo = isqrt(3 * n)                             # x even, 3n+1 <= x^2 <= 4n
    if lo * lo < 3 * n + 1:
        lo += 1
    hi = isqrt(4 * n)
    x = lo + (lo % 2)                             # first even x >= lo
    while x <= hi:
        if isprime(m - x) and isprime(m + x):
            return False
        x += 2
    return True

def sha(fn):
    return hashlib.sha256(open(fn, "rb").read()).hexdigest()

def load(fn):
    with open(fn) as f:
        r = csv.reader(f); next(r)                # header 'n'
        return [int(row[0]) for row in r]

def main():
    fj = sys.argv[1] if len(sys.argv) > 1 else "zeros_J.csv"
    fr = sys.argv[2] if len(sys.argv) > 2 else "zeros_R.csv"
    for name, fn, empt in [("J", fj, lambda n: lower_empty(n) and upper_empty(n)),
                           ("R", fr, lower_empty)]:
        ns = load(fn)
        print(f"[{name}] {fn}: {len(ns)} rows, sha256 {sha(fn)[:16]}..., last {ns[-1]}")
        bad = [n for n in ns if not empt(n)]
        print(f"[{name}] all listed indices verified empty: {'OK' if not bad else f'FAIL {bad[:5]}'}")
        other = set(load(fj)) | set(load(fr))
        wrongly_missing = [m for m in range(ns[-1] + 1, ns[-1] + 51)
                           if empt(m) and m not in other]
        print(f"[{name}] spot check above record ({ns[-1]}+1..+50): "
              f"{'OK (no unlisted empties)' if not wrongly_missing else f'FAIL {wrongly_missing}'}")

if __name__ == "__main__":
    main()
