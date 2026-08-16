#!/usr/bin/env python3
"""
campaign_ledger.py -- Parse, verify, and price a Fortunate-semiprimes
deferred-target worklist (next_targets.txt format).

The worklist format:
  # ---- n=NN: winner a_prp=AAA (prime-m|composite-m); K deferred ----
  # n=NN m=MMM (D digits)
  <full decimal expansion of Q_NN + MMM>

For each index this script:

  1. Recomputes Q_n exactly (product of the first n semiprimes) and checks
     every listed integer equals Q_n + m exactly.
  2. Checks structural claims: each deferred m is prime and exceeds
     y = sp(n)/2 (Lemma 3: deferred offsets are always prime, and an
     inherited prime m <= y would never be deferred); each deferred m lies
     below the winner; the stated digit count matches.
  3. Checks the winner's shape: a prime-m winner must be prime; a
     composite-m winner must have the Lemma 4(i) form p*Q with
     p <= y < Q, both prime. For composite-m winners the semiprimality
     of Q_n + a is additionally screened: the cofactor (Q_n + a)/p is
     BPSW-tested (sympy.isprime), which verifies the "prp" in a_prp.
     Prime-m winners rest on factorizations not present in this file and
     are marked unverifiable-here.
  4. Prices the race under Model 1: every deferred offset is a prime
     m > y, so Q_n + m is automatically y-rough and its semiprime hazard
     is h = (e^gamma log y / L) log(u - 1) with L = log Q_n, u = L/log y.
     Reports h, the expected number of semiprimes hiding among the
     deferred offsets (k*h), and the survival probability of the winner
     (1-h)^k -- for composite-m winners this is the model price of a
     composite term at that index.
  5. Emits two CSVs:
       deferred_ledger.csv    -- one row per deferred target
                                 (n, m, digits, status=open)
       values_provisional.csv -- one row per index
                                 (n, a, status=provisional)
     Fill in the ledger with DEAD / SEMIPRIME verdicts as factors arrive;
     an index whose every deferred m below the (possibly updated) winner
     is DEAD becomes exact, its row in the values file flips to
     proved/probable, and the file then feeds hazard_transform_test.py.

Usage:
  python3 campaign_ledger.py next_targets.txt [--outdir DIR]
"""

import argparse
import csv
import math
import os
import re
import sys

from sympy import isprime, factorint

EULER_GAMMA = 0.5772156649015328606

HEADER_RE = re.compile(
    r"#\s*----\s*n=(\d+):\s*winner\s+a_prp=(\d+)\s*\((prime-m|composite-m)\);"
    r"\s*(\d+)\s+deferred\s*----"
)
TARGET_RE = re.compile(r"#\s*n=(\d+)\s+m=(\d+)\s*\((\d+)\s+digits\)")


def semiprimes_up_to_index(nmax):
    sps = []
    k = 3
    while len(sps) < nmax:
        k += 1
        if sum(factorint(k).values()) == 2:
            sps.append(k)
    return sps


def parse_worklist(path):
    """Return {n: {'winner': a, 'kind': str, 'declared': k,
                   'targets': [(m, digits, bigint), ...]}}"""
    indices = {}
    current = None
    pending = None  # (n, m, digits) awaiting its integer line
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            mh = HEADER_RE.match(line)
            if mh:
                n = int(mh.group(1))
                indices[n] = {
                    "winner": int(mh.group(2)),
                    "kind": mh.group(3),
                    "declared": int(mh.group(4)),
                    "targets": [],
                }
                current = n
                pending = None
                continue
            mt = TARGET_RE.match(line)
            if mt:
                pending = (int(mt.group(1)), int(mt.group(2)),
                           int(mt.group(3)))
                continue
            if line.startswith("#"):
                continue
            if line.isdigit():
                if pending is None or current is None:
                    sys.exit(f"line {lineno}: integer with no target header")
                tn, tm, td = pending
                if tn != current:
                    sys.exit(f"line {lineno}: target n={tn} under "
                             f"section n={current}")
                indices[current]["targets"].append((tm, td, int(line)))
                pending = None
    return indices


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("worklist")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    indices = parse_worklist(args.worklist)
    if not indices:
        sys.exit("No indices parsed.")
    nmax = max(indices)
    sps = semiprimes_up_to_index(nmax)

    problems = []
    rows_ledger = []
    rows_values = []

    print(f"{'n':>4} {'sp(n)':>6} {'y':>7} {'winner':>7} {'kind':>11} "
          f"{'#def':>4} {'h_II':>6} {'E[semi]':>7} {'P(stand)':>8} "
          f"{'verify':>7}")

    for n in sorted(indices):
        info = indices[n]
        sp = sps[n - 1]
        y = sp / 2.0
        Q = 1
        for q in sps[:n]:
            Q *= q
        L = sum(math.log(q) for q in sps[:n])
        logy = math.log(y)
        u = L / logy
        h = math.exp(EULER_GAMMA) * logy / L * math.log(u - 1.0)

        ok = True
        targets = info["targets"]
        if len(targets) != info["declared"]:
            ok = False
            problems.append(f"n={n}: header declares {info['declared']} "
                            f"deferred, file lists {len(targets)}")

        seen_m = set()
        for (m, digits, big) in targets:
            if Q + m != big:
                ok = False
                problems.append(f"n={n} m={m}: integer != Q_n + m "
                                f"(differs by {big - Q - m})")
            if len(str(big)) != digits:
                ok = False
                problems.append(f"n={n} m={m}: stated {digits} digits, "
                                f"actual {len(str(big))}")
            if not isprime(m):
                ok = False
                problems.append(f"n={n} m={m}: deferred offset not prime "
                                f"(violates Lemma 3 bookkeeping)")
            if m <= y:
                ok = False
                problems.append(f"n={n} m={m}: deferred offset <= y={y} "
                                f"(inherited primes are never deferred)")
            if m >= info["winner"]:
                ok = False
                problems.append(f"n={n} m={m}: deferred offset not below "
                                f"winner {info['winner']}")
            if m in seen_m:
                ok = False
                problems.append(f"n={n} m={m}: duplicate target")
            seen_m.add(m)
            rows_ledger.append([n, m, digits, "open", ""])

        # Winner shape checks
        a = info["winner"]
        if info["kind"] == "prime-m":
            if not isprime(a):
                ok = False
                problems.append(f"n={n}: prime-m winner {a} is composite")
            wverify = "shape"  # semiprimality rests on external factors
        else:
            f = factorint(a)
            shape_ok = (len(f) == 2
                        and all(e == 1 for e in f.values()))
            if shape_ok:
                p, bigq = sorted(f)
                shape_ok = p <= y < bigq
            if not shape_ok:
                ok = False
                problems.append(f"n={n}: composite-m winner {a} not of "
                                f"Lemma 4(i) shape p*Q with p<=y<Q")
                wverify = "FAIL"
            else:
                p = min(f)
                num = Q + a
                if num % p != 0:
                    ok = False
                    problems.append(f"n={n}: p={p} does not divide "
                                    f"Q_n + {a}")
                    wverify = "FAIL"
                else:
                    cof = num // p
                    wverify = "prp-ok" if isprime(cof) else "PRP-FAIL"
                    if wverify == "PRP-FAIL":
                        ok = False
                        problems.append(
                            f"n={n}: cofactor (Q_n+{a})/{p} fails BPSW; "
                            f"winner is not even a probable semiprime")

        k = len(targets)
        p_stand = (1.0 - h) ** k
        rows_values.append([n, a, "provisional"])
        print(f"{n:>4} {sp:>6} {y:>7.1f} {a:>7} {info['kind']:>11} "
              f"{k:>4} {h:>6.3f} {k * h:>7.2f} {p_stand:>8.3f} "
              f"{('OK' if ok else 'ISSUE'):>4}/{wverify}")

    os.makedirs(args.outdir, exist_ok=True)
    lpath = os.path.join(args.outdir, "deferred_ledger.csv")
    with open(lpath, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "m", "digits", "status", "factor"])
        w.writerows(rows_ledger)
    vpath = os.path.join(args.outdir, "values_provisional.csv")
    with open(vpath, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "a", "status"])
        w.writerows(rows_values)

    print(f"\nWrote {lpath} ({len(rows_ledger)} targets) and {vpath} "
          f"({len(rows_values)} indices).")
    print("Ledger statuses: open -> DEAD (factor found, cofactor composite "
          "or offset otherwise excluded)\n                 open -> "
          "SEMIPRIME (factor found, cofactor prp/proved: new winner)")

    if problems:
        print(f"\n*** {len(problems)} verification problem(s):")
        for p in problems:
            print("   ", p)
        sys.exit(1)
    print("\nAll structural and arithmetic checks passed.")


if __name__ == "__main__":
    main()
