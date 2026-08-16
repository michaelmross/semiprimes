#!/usr/bin/env python3
"""
property_scan.py -- property-only scanner for the A226525 / d2 conjecture.

Definitions
-----------
Let s_1 < s_2 < ... be the semiprimes (Omega = 2), and

    Q_n = product_{k=1}^n s_k.

Define

    d2(n) = min { m > 1 : Q_n + m is semiprime }.

The purpose of this program is NOT to determine d2(n) exactly.  It asks only
whether d2(n) must be prime.

Key stopping rule
-----------------
Scan offsets m = 2,3,... upward and classify Q_n+m whenever that can be done
cheaply.  Some rough shifts cannot be distinguished (without factoring) between
semiprime and Omega >= 3; these are marked unresolved.

As soon as a *resolved* semiprime occurs at a PRIME offset p, the property is
certified (at PRP confidence) provided every unresolved earlier offset is also
prime.  The true d2(n) is then either p or one of those earlier unresolved prime
offsets; either way it is prime.

If the first resolved semiprime occurs at a COMPOSITE offset c, then only the
unresolved offsets below c matter.  If none exist, c is a PRP-level
counterexample.  If unresolved prime offsets exist, the row is a CHALLENGE:
any one of them being semiprime saves the conjecture at that n.

Structural acceleration
-----------------------
Put y = floor(s_n/2).  Every prime <= y divides Q_n.  Hence:

  * if m contains two distinct prime divisors <= y, Q_n+m already has at least
    those two prime divisors and a nontrivial cofactor, so Omega >= 3;
  * if m has exactly one such prime divisor p, divide Q_n+m by p exactly.
    If p occurs once, semiprimality is equivalent to primality of the cofactor;
    if it occurs at least twice, the shift is normally dead immediately;
  * only y-rough offsets require the more expensive small-factor hunt.  In the
    usual race range (< q0^2, q0 the least prime > y), every y-rough m > 1 is
    prime, which is exactly why unresolved candidates normally cannot hide a
    composite counterexample.

For y-rough shifts, the scanner looks for a factor <= --trial by taking one gcd
against primorial(trial).  If it finds a factor, it strips all factors from that
gcd and PRP-tests the remaining cofactor.  If no factor <= trial is found and
the whole shifted number is composite, the shift remains unresolved.

Confidence
----------
Large cofactors are tested with gmpy2.is_prime(..., MR_EXTRA).  Therefore
PASS-PRP and COUNTEREXAMPLE-PRP are computational/PRP-level statements, not
formal primality certificates.  The offset m itself is a small integer and its
prime/composite status is exact.

Typical use
-----------
    python3 property_scan.py 1 999

or, to continue the existing repository's range:

    python3 property_scan.py 251 999 --trial 1000000 --cap-mult 8

Outputs
-------
  d2_property_<lo>_<hi>.csv       one compact row per n
  d2_property_<lo>_<hi>.jsonl     full unresolved-offset lists
  d2_challenges_<lo>_<hi>.txt     n,m pairs needing deeper work

Dependencies: gmpy2.  (The existing Fortunate scripts already use it.)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import gmpy2
from gmpy2 import gcd, is_prime, mpz

MR_EXTRA = 24


# ---------------------------------------------------------------------------
# Small integer sieves

def primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    sieve = bytearray(b"\x01") * (B + 1)
    sieve[0:2] = b"\x00\x00"
    r = int(B ** 0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            sieve[p * p : B + 1 : p] = b"\x00" * (((B - p * p) // p) + 1)
    return [i for i in range(2, B + 1) if sieve[i]]


def spf_sieve(B: int) -> List[int]:
    """Smallest-prime-factor table for 0..B."""
    spf = list(range(B + 1))
    if B >= 0:
        spf[0] = 0
    if B >= 1:
        spf[1] = 1
    r = int(B ** 0.5)
    for p in range(2, r + 1):
        if spf[p] == p:
            start = p * p
            for x in range(start, B + 1, p):
                if spf[x] == x:
                    spf[x] = p
    return spf


def distinct_prime_factors_small(x: int, spf: Sequence[int]) -> List[int]:
    out: List[int] = []
    while x > 1:
        p = spf[x]
        out.append(p)
        while x % p == 0:
            x //= p
    return out


def first_semiprimes(nmax: int) -> List[int]:
    """Return the first nmax semiprimes, by an Omega sieve."""
    if nmax <= 0:
        return []
    B = max(64, 4 * nmax)
    while True:
        # Omega(x), counted with multiplicity.
        omega = bytearray(B + 1)
        for p in primes_upto(B):
            pk = p
            while pk <= B:
                for x in range(pk, B + 1, pk):
                    omega[x] += 1
                if pk > B // p:
                    break
                pk *= p
        sps = [x for x in range(4, B + 1) if omega[x] == 2]
        if len(sps) >= nmax:
            return sps[:nmax]
        B *= 2


# ---------------------------------------------------------------------------
# PRP / smooth-factor utilities

def is_prp(x: mpz) -> bool:
    return bool(is_prime(x, MR_EXTRA))


def primorial_upto(B: int, primes: Sequence[int]) -> mpz:
    """Use GMP's primorial if available; otherwise multiply once."""
    try:
        return mpz(gmpy2.primorial(B))
    except AttributeError:
        P = mpz(1)
        for p in primes:
            P *= p
        return P


def strip_gcd_factors(
    h: mpz,
    d: mpz,
    trial_primes: Sequence[int],
    stop_after: int = 3,
) -> Tuple[mpz, int]:
    """
    d = gcd(h, primorial(B)), so d is squarefree and B-smooth.
    Identify its prime divisors and strip their FULL multiplicity from h.
    Return (remaining h, number of prime factors stripped with multiplicity).

    We stop once 'stop_after' factors have been stripped, because for the
    semiprime question Omega >= 3 is already dead.
    """
    dd = int(d)
    stripped = 0

    for p in trial_primes:
        if dd == 1:
            break
        if p * p > dd:
            p = dd
        if dd % p == 0:
            while h % p == 0:
                h //= p
                stripped += 1
                if stripped >= stop_after:
                    return h, stripped
            # d is squarefree, so remove p once from dd.
            dd //= p
        if p * p > dd and dd > 1:
            # dd is now prime.
            p2 = dd
            while h % p2 == 0:
                h //= p2
                stripped += 1
                if stripped >= stop_after:
                    return h, stripped
            dd = 1
            break

    if dd > 1:
        # Defensive fallback; should only happen if trial_primes did not reach B.
        p = dd
        while h % p == 0:
            h //= p
            stripped += 1
            if stripped >= stop_after:
                break

    return h, stripped


# ---------------------------------------------------------------------------
# Candidate classification

# Return values for classify_shift:
#   "SEMIPRIME"    resolved Omega == 2
#   "DEAD"         resolved Omega != 2
#   "UNRESOLVED"   composite rough number, no factor <= trial bound found


def classify_shift(
    Q: mpz,
    m: int,
    y: int,
    spf: Sequence[int],
    trial_primorial: mpz,
    trial_primes: Sequence[int],
) -> str:
    """Classify Q+m only as far as needed for the d2 prime-property test."""
    factors_m = distinct_prime_factors_small(m, spf)
    small = [p for p in factors_m if p <= y]

    # Two distinct inherited prime divisors already force Omega(Q+m) >= 3.
    if len(small) >= 2:
        return "DEAD"

    g = Q + m

    if len(small) == 1:
        # p divides both Q and m, hence p divides g.  Strip it exactly.
        p = small[0]
        h = g
        vp = 0
        while h % p == 0:
            h //= p
            vp += 1

        if vp == 0:
            raise AssertionError(f"inherited divisor p={p} failed for m={m}")

        if h == 1:
            return "SEMIPRIME" if vp == 2 else "DEAD"

        if vp >= 2:
            # At least p^2 times a nontrivial cofactor -> Omega >= 3.
            return "DEAD"

        # Exactly one inherited p.  Then Q+m is semiprime iff h is prime.
        return "SEMIPRIME" if is_prp(h) else "DEAD"

    # No prime divisor of m is <= y.  This is the only genuinely rough case.
    # First look for any factor <= trial bound in one GMP gcd.
    d = gcd(g, trial_primorial)
    if d > 1:
        h, stripped = strip_gcd_factors(g, d, trial_primes, stop_after=3)
        if stripped >= 3:
            return "DEAD"
        if h == 1:
            return "SEMIPRIME" if stripped == 2 else "DEAD"
        if is_prp(h):
            return "SEMIPRIME" if stripped + 1 == 2 else "DEAD"
        # A composite cofactor contributes at least two more prime factors.
        return "DEAD"

    # No factor <= trial bound.  If g itself is prime, it is not semiprime.
    # Otherwise it may be a product of two large primes, or have Omega >= 3.
    return "DEAD" if is_prp(g) else "UNRESOLVED"


# ---------------------------------------------------------------------------
# Per-n property scan

def scan_property_n(
    n: int,
    sp_n: int,
    Q: mpz,
    lnQ: float,
    spf: Sequence[int],
    small_primes: Sequence[int],
    trial_primorial: mpz,
    trial_primes: Sequence[int],
    cap: int,
) -> Dict:
    y = sp_n // 2
    q0 = next(p for p in small_primes if p > y)

    unresolved_prime: List[int] = []
    unresolved_composite: List[int] = []
    witness: Optional[int] = None
    witness_kind: Optional[str] = None

    for m in range(2, cap + 1):
        verdict = classify_shift(
            Q, m, y, spf, trial_primorial, trial_primes
        )

        m_prime = (spf[m] == m)

        if verdict == "UNRESOLVED":
            if m_prime:
                unresolved_prime.append(m)
            else:
                unresolved_composite.append(m)
            continue

        if verdict == "SEMIPRIME":
            witness = m
            witness_kind = "prime-m" if m_prime else "composite-m"
            break

    if witness is None:
        status = "NO-WITNESS"
        reason = f"no resolved semiprime through m={cap}"
    elif witness_kind == "prime-m" and not unresolved_composite:
        status = "PASS-PRP"
        reason = (
            "resolved semiprime at prime offset; every unresolved earlier "
            "offset is also prime"
        )
    elif witness_kind == "composite-m" and not unresolved_prime and not unresolved_composite:
        status = "COUNTEREXAMPLE-PRP"
        reason = "resolved composite semiprime offset with no unresolved earlier shifts"
    else:
        status = "CHALLENGE"
        if witness_kind == "composite-m":
            reason = (
                "composite resolved witness; earlier unresolved shifts must be resolved"
            )
        else:
            reason = (
                "prime resolved witness, but an earlier unresolved composite offset exists"
            )

    return {
        "n": n,
        "sp_n": sp_n,
        "digits_Q": len(str(Q)),
        "lnQ": lnQ,
        "y": y,
        "q0": q0,
        "two_q0": 2 * q0,
        "q0_squared": q0 * q0,
        "cap": cap,
        "status": status,
        "reason": reason,
        "witness_m": witness,
        "witness_kind": witness_kind,
        "unresolved_prime_count": len(unresolved_prime),
        "unresolved_composite_count": len(unresolved_composite),
        "unresolved_prime": unresolved_prime,
        "unresolved_composite": unresolved_composite,
    }


# ---------------------------------------------------------------------------
# Output / driver

def csv_row(r: Dict) -> Dict:
    return {
        "n": r["n"],
        "sp_n": r["sp_n"],
        "digits_Q": r["digits_Q"],
        "y": r["y"],
        "q0": r["q0"],
        "two_q0": r["two_q0"],
        "cap": r["cap"],
        "status": r["status"],
        "witness_m": "" if r["witness_m"] is None else r["witness_m"],
        "witness_kind": "" if r["witness_kind"] is None else r["witness_kind"],
        "open_prime": r["unresolved_prime_count"],
        "open_composite": r["unresolved_composite_count"],
    }


def write_challenge_targets(path: str, records: Iterable[Dict]) -> None:
    with open(path, "w") as fh:
        fh.write("# Property-only d2 challenge targets.\n")
        fh.write("# Resolve Q_n+m only for these offsets; exact d2 values are unnecessary.\n")
        fh.write("# n,m,offset_kind,witness_m,witness_kind\n")
        for r in records:
            if r["status"] != "CHALLENGE":
                continue
            for m in r["unresolved_prime"]:
                fh.write(
                    f"{r['n']},{m},prime,{r['witness_m']},{r['witness_kind']}\n"
                )
            for m in r["unresolved_composite"]:
                fh.write(
                    f"{r['n']},{m},composite,{r['witness_m']},{r['witness_kind']}\n"
                )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Property-only scanner: test whether d2(n) is prime without computing d2 exactly."
    )
    ap.add_argument("n_lo", nargs="?", type=int, default=1)
    ap.add_argument("n_hi", nargs="?", type=int, default=999)
    ap.add_argument(
        "--trial",
        type=int,
        default=10**6,
        help="small-factor bound for rough shifts (default: 1,000,000)",
    )
    ap.add_argument(
        "--cap-mult",
        type=float,
        default=8.0,
        help="scan cap = floor(cap_mult*log(Q_n))+50 unless --cap is given",
    )
    ap.add_argument(
        "--cap",
        type=int,
        default=None,
        help="fixed maximum offset for every n (overrides --cap-mult)",
    )
    ap.add_argument(
        "--out-prefix",
        default=None,
        help="output prefix (default: d2_property_<lo>_<hi>)",
    )
    ap.add_argument(
        "--resume",
        action="store_true",
        help="skip n already present in the CSV output and append new rows",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.n_lo < 1 or args.n_hi < args.n_lo:
        raise SystemExit("Require 1 <= n_lo <= n_hi")
    if args.trial < 2:
        raise SystemExit("--trial must be >= 2")

    prefix = args.out_prefix or f"d2_property_{args.n_lo}_{args.n_hi}"
    csv_path = prefix + ".csv"
    jsonl_path = prefix + ".jsonl"
    challenge_path = prefix.replace("d2_property", "d2_challenges") + ".txt"

    print(f"Generating first {args.n_hi} semiprimes...", file=sys.stderr)
    sps = first_semiprimes(args.n_hi)

    # Compute log(Q_nmax) first so we know the largest small-offset SPF table needed.
    lnQmax = sum(math.log(s) for s in sps)
    global_cap = args.cap if args.cap is not None else int(args.cap_mult * lnQmax) + 50
    global_cap = max(global_cap, 100)

    # q0 is only a little above sp(n)/2, but include enough headroom in this prime list.
    small_prime_bound = max(global_cap, sps[-1] + 1000)
    print(f"Sieving small integers through {small_prime_bound:,}...", file=sys.stderr)
    spf = spf_sieve(small_prime_bound)
    small_primes = [i for i in range(2, small_prime_bound + 1) if spf[i] == i]

    print(f"Building primorial({args.trial:,}) for rough-shift gcds...", file=sys.stderr)
    trial_primes = primes_upto(args.trial)
    trial_primorial = primorial_upto(args.trial, trial_primes)
    print(
        f"trial primes: {len(trial_primes):,}; max scan cap: {global_cap:,}",
        file=sys.stderr,
    )

    completed = set()
    old_records: List[Dict] = []
    if args.resume and os.path.exists(csv_path):
        with open(csv_path, newline="") as fh:
            for row in csv.DictReader(fh):
                completed.add(int(row["n"]))
        # Preserve full old JSONL records if available, so challenge file can be rebuilt.
        if os.path.exists(jsonl_path):
            with open(jsonl_path) as fh:
                for line in fh:
                    if line.strip():
                        old_records.append(json.loads(line))
        print(f"Resume: {len(completed)} indices already present.", file=sys.stderr)

    fieldnames = [
        "n", "sp_n", "digits_Q", "y", "q0", "two_q0", "cap", "status",
        "witness_m", "witness_kind", "open_prime", "open_composite",
    ]

    csv_mode = "a" if (args.resume and os.path.exists(csv_path)) else "w"
    json_mode = "a" if (args.resume and os.path.exists(jsonl_path)) else "w"

    new_records: List[Dict] = []
    t0 = time.time()
    Q = mpz(1)
    lnQ = 0.0

    with open(csv_path, csv_mode, newline="") as cf, open(jsonl_path, json_mode) as jf:
        cw = csv.DictWriter(cf, fieldnames=fieldnames)
        if csv_mode == "w":
            cw.writeheader()
            cf.flush()

        for n, sp_n in enumerate(sps, 1):
            Q *= sp_n
            lnQ += math.log(sp_n)
            if n < args.n_lo:
                continue
            if n > args.n_hi:
                break
            if n in completed:
                continue

            cap = args.cap if args.cap is not None else int(args.cap_mult * lnQ) + 50
            cap = max(cap, 100)
            if cap > global_cap:
                raise AssertionError("internal cap exceeded SPF table")

            t1 = time.time()
            r = scan_property_n(
                n=n,
                sp_n=sp_n,
                Q=Q,
                lnQ=lnQ,
                spf=spf,
                small_primes=small_primes,
                trial_primorial=trial_primorial,
                trial_primes=trial_primes,
                cap=cap,
            )
            r["seconds"] = time.time() - t1
            new_records.append(r)

            cw.writerow(csv_row(r))
            cf.flush()
            jf.write(json.dumps(r, separators=(",", ":")) + "\n")
            jf.flush()

            w = "-" if r["witness_m"] is None else str(r["witness_m"])
            print(
                f"n={n:4d} sp={sp_n:5d} dig={r['digits_Q']:4d} "
                f"2q0={r['two_q0']:5d} witness={w:>6} "
                f"{(r['witness_kind'] or '-'):>11} "
                f"openP={r['unresolved_prime_count']:4d} "
                f"openC={r['unresolved_composite_count']:2d} "
                f"{r['status']:<18} {r['seconds']:.2f}s",
                flush=True,
            )

            # A PRP-level counterexample deserves an immediate, conspicuous stop.
            if r["status"] == "COUNTEREXAMPLE-PRP":
                print("*** PRP-LEVEL COUNTEREXAMPLE FOUND; stopping. ***", file=sys.stderr)
                break

    all_records = old_records + new_records
    write_challenge_targets(challenge_path, all_records)

    counts: Dict[str, int] = {}
    for r in all_records:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    dt = time.time() - t0
    print("\nSummary", file=sys.stderr)
    for k in sorted(counts):
        print(f"  {k:20s} {counts[k]:5d}", file=sys.stderr)
    print(f"  elapsed              {dt:.1f}s", file=sys.stderr)
    print(f"  CSV:       {csv_path}", file=sys.stderr)
    print(f"  JSONL:     {jsonl_path}", file=sys.stderr)
    print(f"  challenges:{challenge_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
