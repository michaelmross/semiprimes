#!/usr/bin/env python3
"""
kp-semiprimes in quadratic intervals (n^2, (n+1)^2).

For each prime k, counts or verifies existence of semiprimes k*p (p prime).
Two modes:
  --mode existence  Fast check S_k(n) >= 1 via Miller-Rabin. Handles n to 10^7+.
  --mode count      Full count via segmented sieve. Handles n to ~500K.

Usage:
    python kp_semiprimes.py --nmax 10000000 --mode existence
    python kp_semiprimes.py --nmax 100000 --mode count --output results.csv
    python kp_semiprimes.py --nmax 500000 --klist 2,3,5,7,11 --mode count
"""

import argparse
import math
import time
import csv
import sys

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


# ── Miller-Rabin (deterministic for n < 3.317e24) ─────────────────────────

_SMALL_PRIMES = frozenset({2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37})
_WITNESSES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)

def is_prime(n):
    """Deterministic primality test for n < 3.317 * 10^24."""
    if n < 2:
        return False
    if n in _SMALL_PRIMES:
        return True
    if n % 2 == 0 or n % 3 == 0 or n % 5 == 0:
        return False
    if n < 49:
        return True

    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    for a in _WITNESSES:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# ── Sieve utilities ────────────────────────────────────────────────────────

def sieve_primes(limit):
    """Return list of primes up to limit."""
    if limit < 2:
        return []
    if HAS_NUMPY:
        s = np.ones(limit + 1, dtype=np.bool_)
        s[:2] = False
        for i in range(2, int(math.isqrt(limit)) + 1):
            if s[i]:
                s[i*i::i] = False
        return np.nonzero(s)[0].tolist()
    s = bytearray(b'\x01') * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if s[i]:
            s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(2, limit + 1) if s[i]]


def sieve_count(lo, hi, base_primes):
    """Count primes in [lo, hi) using segmented sieve."""
    if hi <= lo or hi <= 2:
        return 0
    lo = max(lo, 2)
    size = hi - lo

    if HAS_NUMPY:
        s = np.ones(size, dtype=np.bool_)
        for p in base_primes:
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p) - lo
            s[start::p] = False
        return int(np.count_nonzero(s))
    else:
        s = bytearray(b'\x01') * size
        for p in base_primes:
            if p * p >= hi:
                break
            start = max(p * p, ((lo + p - 1) // p) * p) - lo
            s[start::p] = bytearray(len(s[start::p]))
        return sum(s)


# ── Core computation ───────────────────────────────────────────────────────

def kp_bounds(n, k):
    """Return (p_lo, p_hi) such that primes p in [p_lo, p_hi)
    correspond to kp-semiprimes k*p in (n^2, (n+1)^2)."""
    nn1 = n * n
    nn2 = (n + 1) * (n + 1)
    p_lo = nn1 // k + 1            # smallest p with k*p > n^2
    p_hi = (nn2 - 1) // k + 1      # largest valid p + 1
    return p_lo, p_hi


def existence_check(n, k):
    """Check if at least one kp-semiprime exists in (n^2, (n+1)^2)."""
    p_lo, p_hi = kp_bounds(n, k)
    for p in range(p_lo, p_hi):
        if is_prime(p):
            return True
    return False


def count_kp(n, k, base_primes):
    """Count kp-semiprimes in (n^2, (n+1)^2) via segmented sieve."""
    p_lo, p_hi = kp_bounds(n, k)
    return sieve_count(p_lo, p_hi, base_primes)


# ── Run modes ──────────────────────────────────────────────────────────────

def run_existence(nmin, nmax, ks, report_every):
    stats = {k: {"failures": []} for k in ks}
    t_start = time.time()

    for n in range(nmin, nmax + 1):
        for k in ks:
            if not existence_check(n, k):
                stats[k]["failures"].append(n)

        if (n - nmin + 1) % report_every == 0 or n == nmax:
            elapsed = time.time() - t_start
            rate = (n - nmin + 1) / elapsed if elapsed > 0 else 0
            eta = (nmax - n) / rate if rate > 0 else 0
            pct = 100 * (n - nmin + 1) / (nmax - nmin + 1)
            print(f"n = {n:>12,d}  ({pct:5.1f}%)  "
                  f"rate: {rate:,.0f}/s  "
                  f"elapsed: {elapsed:,.0f}s  "
                  f"ETA: {eta:,.0f}s")
            for k in ks:
                nf = len(stats[k]["failures"])
                last = stats[k]["failures"][-1] if nf else None
                s = f"    k={k:>2}: "
                s += "no failures" if nf == 0 else f"{nf} failures (last n={last})"
                print(s)
            print(flush=True)

    return stats


def run_count(nmin, nmax, ks, report_every, output_file):
    sieve_limit = nmax + 2
    print(f"Building base sieve to {sieve_limit}...", end=" ", flush=True)
    t0 = time.time()
    base_primes = sieve_primes(sieve_limit)
    print(f"done ({len(base_primes)} primes, {time.time()-t0:.1f}s)\n")

    stats = {}
    for k in ks:
        stats[k] = {
            "total": 0, "min_count": float("inf"), "min_at": 0,
            "failures": [], "block_mins": [],
            "_bmin": float("inf"), "_bstart": nmin,
        }

    csv_out, writer = None, None
    if output_file:
        csv_out = open(output_file, "w", newline="")
        writer = csv.writer(csv_out)
        writer.writerow(["n"] + [f"S_{k}" for k in ks])

    t_start = time.time()

    for n in range(nmin, nmax + 1):
        counts = []
        for k in ks:
            c = count_kp(n, k, base_primes)
            st = stats[k]
            st["total"] += c
            if c < st["min_count"]:
                st["min_count"] = c
                st["min_at"] = n
            if c == 0:
                st["failures"].append(n)
            if c < st["_bmin"]:
                st["_bmin"] = c
            counts.append(c)

        if writer:
            writer.writerow([n] + counts)

        if (n - nmin + 1) % report_every == 0 or n == nmax:
            for k in ks:
                st = stats[k]
                st["block_mins"].append((st["_bstart"], n, st["_bmin"]))
                st["_bmin"] = float("inf")
                st["_bstart"] = n + 1

            elapsed = time.time() - t_start
            rate = (n - nmin + 1) / elapsed if elapsed > 0 else 0
            eta = (nmax - n) / rate if rate > 0 else 0
            pct = 100 * (n - nmin + 1) / (nmax - nmin + 1)
            print(f"n = {n:>12,d}  ({pct:5.1f}%)  "
                  f"rate: {rate:,.0f}/s  "
                  f"elapsed: {elapsed:,.0f}s  "
                  f"ETA: {eta:,.0f}s")
            for k in ks:
                st = stats[k]
                nf = len(st["failures"])
                last = st["failures"][-1] if nf else None
                avg = st["total"] / (n - nmin + 1)
                s = f"    k={k:>2}: avg={avg:>8.1f}, min={st['min_count']}, "
                s += "no failures" if nf == 0 else f"{nf} failures (last n={last})"
                print(s)
            print(flush=True)

    if csv_out:
        csv_out.close()
        print(f"Results written to {output_file}")

    return stats


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="kp-semiprimes in quadratic intervals (n², (n+1)²)")
    parser.add_argument("--nmin", type=int, default=2)
    parser.add_argument("--nmax", type=int, default=1_000_000)
    parser.add_argument("--klist", type=str, default="2,3",
                        help="Comma-separated k values (default: 2,3)")
    parser.add_argument("--mode", choices=["existence", "count"],
                        default="existence")
    parser.add_argument("--report-every", type=int, default=100_000)
    parser.add_argument("--output", type=str, default=None,
                        help="CSV output file (count mode only)")
    args = parser.parse_args()

    ks = [int(x) for x in args.klist.split(",")]
    t_global = time.time()

    print(f"kp-semiprimes in (n², (n+1)²)")
    print(f"  range:  n ∈ [{args.nmin}, {args.nmax}]")
    print(f"  k:      {ks}")
    print(f"  mode:   {args.mode}")
    print(f"  numpy:  {'yes' if HAS_NUMPY else 'no'}")
    print()

    if args.mode == "existence":
        stats = run_existence(args.nmin, args.nmax, ks, args.report_every)
    else:
        stats = run_count(args.nmin, args.nmax, ks,
                          args.report_every, args.output)

    # Final summary
    total_time = time.time() - t_global
    print("=" * 72)
    print("FINAL SUMMARY")
    print("=" * 72)
    print(f"Range: n = {args.nmin:,} to {args.nmax:,} "
          f"({args.nmax - args.nmin + 1:,} intervals)")
    print(f"Mode:  {args.mode}")
    print()

    for k in ks:
        st = stats[k]
        nf = len(st["failures"])
        print(f"k = {k}")

        if nf == 0:
            print(f"  S_{k}(n) >= 1 for ALL n in [{args.nmin}, {args.nmax}]")
        else:
            print(f"  {nf} failures (S_{k}(n) = 0)")
            print(f"  Last failure: n = {st['failures'][-1]}")
            if nf <= 30:
                print(f"  All: {st['failures']}")
            else:
                print(f"  First 15: {st['failures'][:15]}")
                print(f"  Last  5:  {st['failures'][-5:]}")

        if args.mode == "count":
            avg = st["total"] / (args.nmax - args.nmin + 1)
            print(f"  Mean S_{k}(n): {avg:.2f}")
            print(f"  Min S_{k}(n):  {st['min_count']} (at n={st['min_at']})")
            bm = st.get("block_mins", [])
            if bm:
                print(f"  Block minimums:")
                show = bm[:3] + [None] + bm[-3:] if len(bm) > 7 else bm
                for entry in show:
                    if entry is None:
                        print(f"    ...")
                    else:
                        bs, be, bv = entry
                        print(f"    n in [{bs:>10,}, {be:>10,}]: min = {bv}")
        print()

    print(f"Total runtime: {total_time:,.1f}s")


if __name__ == "__main__":
    main()
