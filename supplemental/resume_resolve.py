#!/usr/bin/env python3
"""
resume_resolve.py -- Continue the tiered resolution, reconstructing state
from derived_state.json + resolutions.csv (append-only). Runs until its
per-invocation deadline, then exits cleanly; re-invoke to continue.

Usage: python3 resume_resolve.py SECONDS TIER
  TIER in {T2, T3, T3b}:
    T2  = P-1 @ 1e7 then ECM 120 @ 5e4
    T3  = ECM 100 @ 25e4
    T3b = ECM 250 @ 1e6   (deep pass for stragglers)
"""

import csv
import json
import subprocess
import sys
import time

import gmpy2
from gmpy2 import mpz

from sympy import factorint

DEADLINE = time.time() + float(sys.argv[1])
TIER = sys.argv[2]

PAPER_WINNERS = {81: 151, 83: 233, 84: 337}
PAPER_CONFIRMED = {81, 83}


def is_prp(x):
    return gmpy2.is_prime(mpz(x), 30)


def semiprimes_up_to_index(nmax):
    sps = []
    k = 3
    while len(sps) < nmax:
        k += 1
        if sum(factorint(k).values()) == 2:
            sps.append(k)
    return sps


def rebuild_state():
    derived = json.load(open("derived_state.json"))
    try:
        res = list(csv.DictReader(open("resolutions.csv")))
    except FileNotFoundError:
        res = []
    state = {}
    for k, rec in derived.items():
        n = int(k)
        winner = PAPER_WINNERS.get(n, rec["winner"])
        wconf = "paper" if n in PAPER_CONFIRMED else None
        # replay verdicts in recorded order
        verdicts = [r for r in res if int(r["n"]) == n]
        dead = set()
        for r in verdicts:
            m = int(r["m"])
            if r["verdict"] == "DEAD":
                dead.add(m)
            elif r["verdict"] == "SEMIPRIME":
                if m < winner:
                    winner = m
                    wconf = f"{r['tier']}:{r['factor']}"
                elif m == winner:
                    wconf = f"{r['tier']}:{r['factor']}"
        opens = sorted(m for m in rec["deferred"]
                       if m < winner and m not in dead)
        state[n] = {"winner": winner, "open": opens, "wconf": wconf}
    return state, res


def try_ecm(N, args, timeout):
    try:
        r = subprocess.run(["ecm", "-q", "-one"] + args, input=str(N),
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    best = None
    for tok in r.stdout.split():
        if tok.isdigit():
            t = mpz(tok)
            if 1 < t < N and N % t == 0:
                if best is None or t < best:
                    best = t
    return best


def runner(N):
    if TIER == "T2":
        d = try_ecm(N, ["-pm1", "10000000"], 90)
        if d is not None:
            return d
        return try_ecm(N, ["-c", "120", "50000"], 300)
    if TIER == "T3":
        return try_ecm(N, ["-c", "100", "250000"], 600)
    if TIER == "T3b":
        return try_ecm(N, ["-c", "250", "1000000"], 3000)
    raise SystemExit(f"unknown tier {TIER}")


def main():
    sps = semiprimes_up_to_index(92)
    state, res = rebuild_state()
    done_keys = {(int(r["n"]), int(r["m"]), r["tier"]) for r in res}
    import os
    if os.path.exists("attempts.csv"):
        for line in open("attempts.csv"):
            a, b, c = line.strip().split(",")
            done_keys.add((int(a), int(b), c))
    Qc = {}

    def Q(n):
        if n not in Qc:
            v = mpz(1)
            for q in sps[:n]:
                v *= q
            Qc[n] = v
        return Qc[n]

    def record(n, m, verdict, factor, cp):
        res.append({"n": n, "m": m, "tier": TIER, "verdict": verdict,
                    "factor": factor, "cofactor_prp": cp})
        with open("resolutions.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["n", "m", "tier", "verdict",
                                               "factor", "cofactor_prp"])
            w.writeheader()
            w.writerows(res)

    for n in sorted(state):
        st = state[n]
        slots = [(m, "d") for m in list(st["open"])]
        if st["wconf"] is None:
            slots.append((st["winner"], "w"))
        for m, kind in slots:
            if time.time() > DEADLINE:
                print("[deadline]", flush=True)
                report(state)
                return
            if kind == "d" and (m >= st["winner"] or m not in st["open"]):
                continue
            if kind == "w" and st["wconf"] is not None:
                continue
            if (n, m, TIER) in done_keys:
                continue  # this tier already attempted on this target
            N = Q(n) + m
            t0 = time.time()
            d = runner(N)
            if d is None:
                done_keys.add((n, m, TIER))
                with open("attempts.csv", "a") as fh:
                    fh.write(f"{n},{m},{TIER}\n")
                print(f"n={n} m={m}: survived {TIER} "
                      f"({time.time()-t0:.0f}s)", flush=True)
                continue
            c = N // d
            if is_prp(d) and is_prp(c):
                verdict = "SEMIPRIME"
            else:
                verdict = "DEAD"
            record(n, m, verdict, str(d), is_prp(c) and is_prp(d))
            if kind == "w":
                st["wconf"] = f"{TIER}:{d}"
                print(f"n={n}: winner {m} confirmed ({d})", flush=True)
            elif verdict == "DEAD":
                st["open"].remove(m)
                print(f"n={n} m={m}: DEAD (factor {d})", flush=True)
            else:
                print(f"*** n={n}: DISPLACEMENT to {m} (factor {d})",
                      flush=True)
                st["winner"] = m
                st["wconf"] = f"{TIER}:{d}"
                st["open"] = [x for x in st["open"] if x < m]
    report(state)


def report(state):
    tot = sum(len(s["open"]) for s in state.values())
    unconf = [n for n, s in state.items() if s["wconf"] is None]
    print(f"[state] open: {tot} across "
          f"{[(n, s['open']) for n, s in sorted(state.items()) if s['open']]}")
    print(f"[state] winners unconfirmed: {unconf}", flush=True)


if __name__ == "__main__":
    main()
