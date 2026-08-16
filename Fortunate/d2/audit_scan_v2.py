#!/usr/bin/env python3
"""
audit_scan_v2.py -- Independent audit of fortunate_scan.py scan JSON files.

This version fixes the original audit_scan.py mismatch in the "nasty" set.
The scanner does not call every prime offset m > sp(n)/2 a nasty: it first
strips any prime divisor of Q_n + m in (sp(n)/2, TRIAL_B].  Only a rough
prime offset whose shifted value is still composite after that medium-prime
stage is deferred as a nasty.

For every record in a scan_*.json file, this program independently:

  * reconstructs the n-th semiprimorial Q_n;
  * checks sp_n, q0, 2*q0, digits_N, and lnN;
  * rescans every offset 2 <= m < the recorded scan front-runner;
  * strips inherited prime factors p <= floor(sp(n)/2);
  * independently strips medium prime factors in
        (floor(sp(n)/2), TRIAL_B];
  * classifies the remaining cofactor at strong-PRP level;
  * detects any smaller resolved semiprime ("UPSET");
  * reconstructs and compares nasties_below exactly;
  * reconstructs and compares alive_composites exactly;
  * verifies that the recorded scan front-runner itself is a resolved
    semiprime at PRP level.

Important:
  This audits the SCAN LAYER.  A nasty below the scan front-runner may later
  be resolved as a semiprime by the factorization campaign, so a_prp in the
  scan JSON need not be the final exact A226525 value.

Confidence:
  gmpy2.is_prime(..., 30) is used for large cofactors.  Thus "semiprime"
  here is at the same strong probable-prime level as the campaign unless
  separate primality certificates are supplied.

Output:
  A human-readable report is written to stdout.
  A derived-state JSON is also written.  Its default filename is unique to
  the input scan file, so several audits may safely run simultaneously.

Examples:
  python3 audit_scan_v2.py scan_1_65.json
  python3 audit_scan_v2.py scan_66_92.json
  python3 audit_scan_v2.py scan_93_117.json

  # Explicit output path:
  python3 audit_scan_v2.py scan_1_65.json --derived-out derived_1_65.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import gmpy2
from gmpy2 import gcd, mpz
from sympy import factorint


TRIAL_B_DEFAULT = 10**6
PRP_ROUNDS = 30


def is_prp(x: int | mpz) -> bool:
    return bool(gmpy2.is_prime(mpz(x), PRP_ROUNDS))


def primes_upto(limit: int) -> List[int]:
    """Simple sieve, used for q0 lookup and a fallback primorial builder."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    r = int(limit**0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = (
                b"\x00" * (((limit - p * p) // p) + 1)
            )
    return [i for i in range(2, limit + 1) if sieve[i]]


def primorial_upto(limit: int, primes: List[int]) -> mpz:
    """
    Product of primes <= limit.

    gmpy2.primorial is preferred; the fallback keeps the script portable
    across older gmpy2 builds.
    """
    if limit < 2:
        return mpz(1)
    if hasattr(gmpy2, "primorial"):
        return mpz(gmpy2.primorial(limit))
    out = mpz(1)
    for p in primes:
        if p > limit:
            break
        out *= p
    return out


def omega_small_exact(m: int) -> int:
    """Exact Omega(m) for the small offsets occurring in these scans."""
    return sum(factorint(m).values())


def semiprimes_up_to_index(nmax: int) -> List[int]:
    """
    Generate the first nmax semiprimes exactly.

    The scan ranges in this project are small enough that factorint-based
    generation is inexpensive and independent of fortunate_scan.py's sieve.
    """
    sps: List[int] = []
    k = 3
    while len(sps) < nmax:
        k += 1
        if omega_small_exact(k) == 2:
            sps.append(k)
    return sps


def strip_known_prime(h: mpz, p: int) -> Tuple[mpz, int]:
    """Strip the full p-adic valuation from h."""
    count = 0
    pp = mpz(p)
    while h % pp == 0:
        h //= pp
        count += 1
    return h, count


def strip_medium(h: mpz, midprod: mpz) -> Tuple[mpz, int, List[int]]:
    """
    Strip all prime factors of h represented in the squarefree product midprod.

    midprod = product_{half < p <= TRIAL_B} p.

    gcd(h, midprod) contains each distinct medium prime at most once.  We
    factor that gcd, then strip the full multiplicity of each such prime from h.
    """
    d = gcd(h, midprod)
    if d == 1:
        return h, 0, []

    # d is composed entirely of primes <= TRIAL_B, so factorint is quick here.
    fs = sorted(int(p) for p in factorint(int(d)))
    count = 0
    for p in fs:
        h, c = strip_known_prime(h, p)
        count += c
    return h, count, fs


def structural_kind(m: int, half: int, sp: int) -> str:
    """A descriptive offset-shape label, independent of the scan outcome."""
    fm = factorint(m)
    if len(fm) == 1 and next(iter(fm.values())) == 1:
        return "inherited" if m <= half else "rough-prime"

    small = [p for p in fm if p <= half]
    if len(small) == 1:
        p = small[0]
        e = fm[p]
        r = m // (p**e)
        if e == 1 and r > half and is_prp(r):
            return "admitted"
        if e == 2 and r == 1 and sp / 3.0 < p <= half:
            return "case-ii"
    return "auto-dead"


def classify_offset(
    Q: mpz,
    m: int,
    sp: int,
    half: int,
    midprod: mpz,
) -> Dict[str, object]:
    """
    Reproduce the mathematical classification used by fortunate_scan.py,
    but from scratch.

    Returns enough internal state to audit nasties, alive composites, and
    accidental semiprime upsets.
    """
    g = Q + m
    h = mpz(g)
    stripped = 0

    fm = factorint(m)
    m_is_prime = (len(fm) == 1 and next(iter(fm.values())) == 1)
    m_is_composite = not m_is_prime
    kind = structural_kind(m, half, sp)

    # Stage 1: inherited factors p <= half.  By the inheritance lemma,
    # these are exactly the small primes that can divide Q+m.
    for p in sorted(fm):
        if p > half:
            continue
        h, c = strip_known_prime(h, int(p))
        stripped += c
        if stripped > 2:
            return {
                "status": "DEAD",
                "omega": None,
                "stripped": stripped,
                "h": h,
                "m_is_prime": m_is_prime,
                "m_is_composite": m_is_composite,
                "kind": kind,
                "medium_primes": [],
                "reached_prp_stage": False,
            }

    # Stage 2: factors in (half, TRIAL_B].
    h, cmed, med_ps = strip_medium(h, midprod)
    stripped += cmed
    if stripped > 2:
        return {
            "status": "DEAD",
            "omega": None,
            "stripped": stripped,
            "h": h,
            "m_is_prime": m_is_prime,
            "m_is_composite": m_is_composite,
            "kind": kind,
            "medium_primes": med_ps,
            "reached_prp_stage": False,
        }

    # This is the point at which fortunate_scan.py tests the remaining
    # cofactor.  Preserve the same condition used for alive_composites.
    reached_prp_stage = (stripped == 1 and h > 1)

    # Stage 3: cofactor classification.
    if h == 1:
        omega = stripped
        resolved = True
    elif is_prp(h):
        omega = stripped + 1
        resolved = True
    else:
        omega = None
        resolved = (stripped >= 1)

    if resolved:
        status = "SEMIPRIME" if omega == 2 else "DEAD"
    else:
        status = "NASTY"

    return {
        "status": status,
        "omega": omega,
        "stripped": stripped,
        "h": h,
        "m_is_prime": m_is_prime,
        "m_is_composite": m_is_composite,
        "kind": kind,
        "medium_primes": med_ps,
        "reached_prp_stage": reached_prp_stage,
    }


def default_derived_name(scan_path: Path) -> Path:
    stem = scan_path.stem
    suffix = stem[5:] if stem.startswith("scan_") else stem
    return scan_path.with_name(f"derived_state_v2_{suffix}.json")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Independent audit of fortunate_scan.py JSON output."
    )
    ap.add_argument("scan_json", help="scan_*.json file to audit")
    ap.add_argument(
        "--trial-bound",
        type=int,
        default=TRIAL_B_DEFAULT,
        help=f"medium-prime stripping bound (default {TRIAL_B_DEFAULT})",
    )
    ap.add_argument(
        "--derived-out",
        default=None,
        help="derived-state JSON path (default is unique to input file)",
    )
    args = ap.parse_args()

    scan_path = Path(args.scan_json)
    with scan_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not data:
        print("Empty scan file.", file=sys.stderr)
        return 2

    nmax = max(int(rec["n"]) for rec in data)
    sps = semiprimes_up_to_index(nmax)

    # The scanner's medium-prime stage is defined by primes up to TRIAL_B.
    print(f"sieving primes to {args.trial_bound} ...", file=sys.stderr, flush=True)
    primes = primes_upto(args.trial_bound)
    all_prime_product = primorial_upto(args.trial_bound, primes)

    # Build Q_n and log(Q_n) once, incrementally.
    Q_by_n: Dict[int, mpz] = {}
    L_by_n: Dict[int, float] = {}
    Q = mpz(1)
    L = 0.0
    for n, q in enumerate(sps, start=1):
        Q *= q
        L += math.log(q)
        Q_by_n[n] = mpz(Q)
        L_by_n[n] = L

    problems: List[str] = []
    derived: Dict[int, dict] = {}

    header = (
        f"{'n':>3} {'win':>5} {'kind':>5} {'inh':>4} {'adm':>4} "
        f"{'nasty':>5} {'json':>5} {'match':>5} "
        f"{'med-kill':>8} {'checks':>7}"
    )
    print(header)

    for rec in sorted(data, key=lambda r: int(r["n"])):
        n = int(rec["n"])
        sp = sps[n - 1]
        half = sp // 2
        y = sp / 2.0
        Q = Q_by_n[n]
        L = L_by_n[n]

        # ---- metadata ----
        if int(rec["sp_n"]) != sp:
            problems.append(f"n={n}: sp_n {rec['sp_n']} != {sp}")

        q0 = next(p for p in primes if p > half)
        if int(rec["q0"]) != q0:
            problems.append(f"n={n}: q0 {rec['q0']} != {q0}")
        if int(rec["two_q0"]) != 2 * q0:
            problems.append(
                f"n={n}: two_q0 {rec['two_q0']} != {2*q0}"
            )

        qdigits = len(str(Q))
        if int(rec["digits_N"]) != qdigits:
            problems.append(
                f"n={n}: digits_N {rec['digits_N']} != {qdigits}"
            )

        if abs(float(rec["lnN"]) - L) > 1e-6 * max(1.0, L):
            problems.append(f"n={n}: lnN {rec['lnN']} vs {L}")

        winner = rec.get("a_prp")
        if winner is None:
            problems.append(f"n={n}: no recorded scan front-runner")
            continue
        winner = int(winner)

        # Product of primes half < p <= TRIAL_B.
        small_primorial = primorial_upto(half, primes)
        midprod = all_prime_product // small_primorial

        # Sanity check behind the scanner's "rough m is prime" statement.
        if winner >= q0 * q0:
            problems.append(
                f"n={n}: winner={winner} reaches/exceeds q0^2={q0*q0}; "
                "rough composite offsets are no longer automatically impossible"
            )

        derived_nasties: List[int] = []
        derived_alive: List[Tuple[int, str]] = []
        upsets: List[Tuple[str, int]] = []

        inherited_count = 0
        admitted_count = 0
        medium_kills = 0
        checks = 0

        # ---- full independent sweep below the scan front-runner ----
        for m in range(2, winner):
            c = classify_offset(Q, m, sp, half, midprod)
            checks += 1

            if c["kind"] == "inherited":
                inherited_count += 1
            elif c["kind"] in ("admitted", "case-ii"):
                admitted_count += 1

            if c["medium_primes"] and c["status"] == "DEAD":
                medium_kills += 1

            if c["status"] == "SEMIPRIME":
                upsets.append((str(c["kind"]), m))
            elif c["status"] == "NASTY":
                # In the protected race range a nasty should be a rough prime.
                if not c["m_is_prime"]:
                    problems.append(
                        f"n={n}: rough composite m={m} classified NASTY"
                    )
                derived_nasties.append(m)

            if (
                c["m_is_composite"]
                and c["reached_prp_stage"]
                and c["status"] != "SEMIPRIME"
            ):
                derived_alive.append((m, "lost"))

        # ---- verify the scan front-runner itself ----
        wc = classify_offset(Q, winner, sp, half, midprod)
        actual_wkind = (
            "composite-m" if wc["m_is_composite"] else "prime-m"
        )
        if actual_wkind != rec["winner_kind"]:
            problems.append(
                f"n={n}: winner_kind {rec['winner_kind']} but "
                f"winner {winner} is {actual_wkind}"
            )
        if wc["status"] != "SEMIPRIME":
            problems.append(
                f"n={n}: recorded winner {winner} reclassifies as "
                f"{wc['status']} (stripped={wc['stripped']})"
            )

        # ---- compare exact nasty set ----
        json_nasties = sorted(int(x) for x in rec.get("nasties_below", []))
        match = (derived_nasties == json_nasties)
        if not match:
            json_extra = sorted(set(json_nasties) - set(derived_nasties))
            derived_extra = sorted(set(derived_nasties) - set(json_nasties))
            problems.append(
                f"n={n}: nasty mismatch, json-extra={json_extra}, "
                f"derived-extra={derived_extra}"
            )

        # ---- compare alive_composites ----
        json_alive = [(int(m), str(tag)) for m, tag in rec.get("alive_composites", [])]
        if derived_alive != json_alive:
            problems.append(
                f"n={n}: alive_composites mismatch, "
                f"json={json_alive}, derived={derived_alive}"
            )

        for kind, m in upsets:
            problems.append(
                f"n={n}: *** UPSET: {kind} offset {m} is a resolved "
                f"probable semiprime below recorded winner {winner}"
            )

        print(
            f"{n:>3} {winner:>5} {rec['winner_kind'][:5]:>5} "
            f"{inherited_count:>4} {admitted_count:>4} "
            f"{len(derived_nasties):>5} {len(json_nasties):>5} "
            f"{'yes' if match else 'NO':>5} "
            f"{medium_kills:>8} {checks:>7}"
        )

        derived[n] = {
            "sp": sp,
            "half": half,
            "y": y,
            "q0": q0,
            "L": L,
            "winner": winner,
            "winner_kind": rec["winner_kind"],
            "winner_recheck": wc["status"],
            "deferred": derived_nasties,
            "alive_composites": derived_alive,
            "Q_digits": qdigits,
            "trial_bound": args.trial_bound,
        }

    derived_out = (
        Path(args.derived_out)
        if args.derived_out
        else default_derived_name(scan_path)
    )
    with derived_out.open("w", encoding="utf-8") as fh:
        json.dump(derived, fh, indent=1)

    print(f"\nDerived state: {derived_out}")

    if problems:
        print(f"\n*** {len(problems)} problem(s):")
        for p in problems:
            print("   ", p)
        return 1

    print(
        "\nAll metadata, medium-prime filtering, scan classifications, "
        "nasty sets,\nalive-composite sets, front-runners, and PRP cofactor "
        "checks passed.\nNo resolved semiprime occurs below any recorded "
        "scan front-runner."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
