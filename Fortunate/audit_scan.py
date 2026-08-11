#!/usr/bin/env python3
"""
audit_scan.py -- Independent re-derivation of the scan layer for
scan_66_92.json.

For each index this script recomputes Q_n exactly and re-classifies every
offset 2 <= m < winner from scratch:

  auto-dead   : m has >= 2 distinct prime factors <= y, or has exactly one
                small prime in a shape excluded by Lemma 4 (Omega >= 3
                unconditionally, no test needed)
  inherited   : m prime <= y  -> Q_n + m = m*(Q_n/m + 1); BPSW the cofactor
  admitted    : m = p*r, p <= y < r both prime (Lemma 4(i)) ->
                BPSW (Q_n + m)/p
  case-ii     : m = p^2, sp/3 < p <= sp/2 -> BPSW Q_n/p + p
                (asserted absent below these winners; checked anyway)
  deferred    : m prime > y -> undecidable without a factor (the nasties)

Every inherited/admitted/case-ii offset below the winner must have a
composite cofactor, or the recorded winner is wrong (a smaller semiprime
exists). The derived deferred set is compared against the JSON's
nasties_below, the JSON's alive_composites are checked to be a subset of
the derived admitted set, and the metadata (sp_n, q0, 2q0, digits, lnN)
is recomputed. Winners <= y would be self-certifying; winners > y get a
direct BPSW confirmation only if a factor is known, so here their prp
status is re-checked only for inherited/admitted shapes (prime winners
> y are exactly the shapes whose semiprimality needs a nontrivial factor
and are handled by the resolution driver, not this audit).

Primality screening uses gmpy2.is_prime (BPSW + Miller-Rabin rounds);
'prp' below means probable in the paper's sense.

Output: audit report to stdout; derived_state.json for the resolver.
"""

import json
import math
import sys

import gmpy2
from gmpy2 import mpz

from sympy import factorint


def semiprimes_up_to_index(nmax):
    sps = []
    k = 3
    while len(sps) < nmax:
        k += 1
        if sum(factorint(k).values()) == 2:
            sps.append(k)
    return sps


def is_prp(x):
    return gmpy2.is_prime(mpz(x), 30)


def small_primes_upto(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:: i] = b"\x00" * len(range(i * i, limit + 1, i))
    return [i for i in range(limit + 1) if sieve[i]]


def main():
    data = json.load(open(sys.argv[1]))
    nmax = max(rec["n"] for rec in data)
    sps = semiprimes_up_to_index(nmax)
    primes_1k = small_primes_upto(1000)
    prime_set_1k = set(primes_1k)

    problems = []
    derived = {}

    print(f"{'n':>3} {'win':>4} {'kind':>5} {'inh':>4} {'adm':>4} "
          f"{'nasty':>5} {'json':>5} {'match':>5} {'cofactor-checks':>15}")

    for rec in sorted(data, key=lambda r: r["n"]):
        n = rec["n"]
        sp = sps[n - 1]
        y = sp / 2.0
        Q = mpz(1)
        for q in sps[:n]:
            Q *= q
        L = sum(math.log(q) for q in sps[:n])

        # ---- metadata ----
        if rec["sp_n"] != sp:
            problems.append(f"n={n}: sp_n {rec['sp_n']} != {sp}")
        q0 = next(p for p in primes_1k if p > y)
        if rec["q0"] != q0:
            problems.append(f"n={n}: q0 {rec['q0']} != {q0}")
        if rec["two_q0"] != 2 * q0:
            problems.append(f"n={n}: two_q0 mismatch")
        if rec["digits_N"] != len(str(Q)):
            problems.append(f"n={n}: digits {rec['digits_N']} != "
                            f"{len(str(Q))}")
        if abs(rec["lnN"] - L) > 1e-6 * L:
            problems.append(f"n={n}: lnN {rec['lnN']} vs {L}")

        # ---- full sweep below the winner ----
        winner = rec["a_prp"]
        inherited_dead, admitted_dead, deferred = [], [], []
        upsets = []
        for m in range(2, winner):
            if m in prime_set_1k:
                if m <= y:
                    cof = Q // m + 1
                    if is_prp(cof):
                        upsets.append(("inherited", m))
                    else:
                        inherited_dead.append(m)
                else:
                    deferred.append(m)
                continue
            f = factorint(m)
            small = [p for p in f if p <= y]
            if len(small) != 1:
                continue  # >=2 small primes: auto-dead; 0 small: impossible
            p = small[0]
            e = f[p]
            r = m // (p ** e)
            if e == 1 and r > 1 and r in prime_set_1k and r > y:
                cof = (Q + m) // p
                if is_prp(cof):
                    upsets.append(("admitted", m))
                else:
                    admitted_dead.append(m)
            elif e == 2 and r == 1 and sp / 3.0 < p <= y:
                cof = Q // p + p
                if is_prp(cof):
                    upsets.append(("case-ii", m))
                else:
                    admitted_dead.append(m)
            # every other one-small-prime shape: auto-dead by Lemma 4

        # ---- winner shape ----
        wkind = "prime-m" if winner in prime_set_1k or is_prp(winner) \
            else "composite-m"
        if wkind != rec["winner_kind"]:
            problems.append(f"n={n}: winner_kind {rec['winner_kind']} "
                            f"but winner {winner} is {wkind}")
        if wkind == "composite-m":
            f = factorint(winner)
            ok = (len(f) == 2 and all(e == 1 for e in f.values()))
            if ok:
                p, r = sorted(f)
                ok = p <= y < r
            if not ok:
                problems.append(f"n={n}: composite winner {winner} not "
                                f"Lemma 4(i) shape")
            else:
                p = min(f)
                if not is_prp((Q + winner) // p):
                    problems.append(f"n={n}: composite winner {winner} "
                                    f"cofactor fails BPSW")

        # ---- compare with JSON ----
        jn = rec["nasties_below"]
        match = (deferred == sorted(jn))
        if not match:
            extra = set(jn) - set(deferred)
            miss = set(deferred) - set(jn)
            problems.append(f"n={n}: nasty mismatch, json-extra={sorted(extra)},"
                            f" derived-extra={sorted(miss)}")
        for (mm, tag) in rec["alive_composites"]:
            if tag != "lost":
                problems.append(f"n={n}: composite {mm} tagged '{tag}'")
            if mm not in admitted_dead:
                problems.append(f"n={n}: json composite {mm} not in derived "
                                f"admitted-dead set")
        for (kind, mm) in upsets:
            problems.append(f"n={n}: *** UPSET: {kind} offset {mm} has PRP "
                            f"cofactor below recorded winner {winner}")

        checks = len(inherited_dead) + len(admitted_dead)
        print(f"{n:>3} {winner:>4} {rec['winner_kind'][:5]:>5} "
              f"{len(inherited_dead):>4} {len(admitted_dead):>4} "
              f"{len(deferred):>5} {len(jn):>5} "
              f"{'yes' if match else 'NO':>5} {checks:>15}")

        derived[n] = {
            "sp": sp, "y": y, "q0": q0, "L": L,
            "winner": winner, "winner_kind": rec["winner_kind"],
            "deferred": deferred,
            "Q_digits": len(str(Q)),
        }

    with open("derived_state.json", "w") as fh:
        json.dump(derived, fh)

    if problems:
        print(f"\n*** {len(problems)} problem(s):")
        for p in problems:
            print("   ", p)
        sys.exit(1)
    print("\nAll metadata, classifications, and cofactor checks passed: no "
          "inherited or admitted\noffset below any recorded winner is a "
          "probable semiprime, and every derived deferred\nset matches the "
          "JSON exactly.")


if __name__ == "__main__":
    main()
