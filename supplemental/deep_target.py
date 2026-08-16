#!/usr/bin/env python3
"""
deep_target.py N M SECONDS B1 -- deep ECM on one specific Q_n + m.

Resolving a single-target index closes it outright: n=69 (only 137 open)
and n=81 (only 139 open) each become exact on one verdict, either way.
Appends any verdict to resolutions.csv in the same schema.
"""

import csv
import os
import subprocess
import sys
import time

import gmpy2
from gmpy2 import mpz
from sympy import factorint

n, m, budget, B1 = int(sys.argv[1]), int(sys.argv[2]), \
    float(sys.argv[3]), sys.argv[4]


def sps_upto(nm):
    s, k = [], 3
    while len(s) < nm:
        k += 1
        if sum(factorint(k).values()) == 2:
            s.append(k)
    return s


Q = mpz(1)
for q in sps_upto(n)[:n]:
    Q *= q
N = Q + m
print(f"n={n} m={m}: N has {N.num_digits()} digits; ECM B1={B1}, "
      f"budget {budget:.0f}s", flush=True)

deadline = time.time() + budget
curves = 0
while time.time() < deadline:
    left = deadline - time.time()
    batch = max(1, int(left // 2))
    try:
        r = subprocess.run(["ecm", "-q", "-one", "-c", str(batch), B1],
                           input=str(N), capture_output=True, text=True,
                           timeout=left)
    except subprocess.TimeoutExpired:
        break
    curves += batch
    hit = None
    for tok in r.stdout.split():
        if tok.isdigit():
            t = mpz(tok)
            if 1 < t < N and N % t == 0:
                hit = t if hit is None else min(hit, t)
    if hit is not None:
        c = N // hit
        dp, cp = gmpy2.is_prime(hit, 30), gmpy2.is_prime(c, 30)
        verdict = "SEMIPRIME" if (dp and cp) else "DEAD"
        print(f"FACTOR {hit} ({hit.num_digits()}d, prime={dp}); "
              f"cofactor {c.num_digits()}d prime={cp} -> {verdict}",
              flush=True)
        rows = list(csv.DictReader(open("resolutions.csv")))
        rows.append({"n": n, "m": m, "tier": f"deep{B1}",
                     "verdict": verdict, "factor": str(hit),
                     "cofactor_prp": dp and cp})
        with open("resolutions.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["n", "m", "tier", "verdict",
                                               "factor", "cofactor_prp"])
            w.writeheader()
            w.writerows(rows)
        sys.exit(0)
print(f"no factor after ~{curves} curves at B1={B1}", flush=True)
with open("attempts.csv", "a") as fh:
    fh.write(f"{n},{m},deep{B1}\n")
