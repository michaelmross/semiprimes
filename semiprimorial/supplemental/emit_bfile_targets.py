#!/usr/bin/env python3
"""
emit_bfile_targets.py -- build the decidable-window completion targets.

Usage:
    python3 emit_bfile_targets.py out.txt scan_66_92.json scan_93_130.json [...]

Reads fortunate_scan JSON files, selects every index n <= 117 (the decidable
window: digits(N_n) <= 250) whose race has deferred ("nasty") offsets below
the PRP winner, rebuilds N_n, and writes one decimal per deferred offset in
ecm_targets format for ecm_hunt.sh / ecm_audit.py.

Skips n = 81 (its two surviving gates are the posted challenge pair; its
resolved gates need no work). Includes prime-winner indices (needed for the
consecutive b-file extension) and the candidate indices 83, 84, 98, 103
(whose full resolution is a certifiable verdict on the conjecture).
"""
import sys, json
import gmpy2
from gmpy2 import mpz

B = 10**6
s = bytearray([1])*(B+1); s[0:2] = b"\x00\x00"
for i in range(2, int(B**0.5)+1):
    if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
PR = [i for i in range(B+1) if s[i]]
def om(x):
    c, t = 0, x
    for p in PR:
        if p*p > t: break
        while t % p == 0: t //= p; c += 1
    return c + (1 if t > 1 else 0)

def main():
    out, *jsons = sys.argv[1:]
    rows = []
    for j in jsons:
        rows += json.load(open(j))
    rows = {r["n"]: r for r in rows}          # dedupe by n, last wins
    todo = []
    for n in sorted(rows):
        r = rows[n]
        if n > 117 or n == 81:               # decidable window; challenge pair excluded
            continue
        if not r.get("a_prp") or not r.get("nasties_below"):
            continue
        todo.append((n, r))
    # semiprimes once, to the max n needed
    nmax = max(n for n, _ in todo)
    sps, x = [], 3
    while len(sps) < nmax:
        x += 1
        if om(x) == 2: sps.append(x)
    with open(out, "w") as f:
        f.write("# Decidable-window completion targets (n<=117, n!=81).\n")
        f.write("# Every DEAD advances the certified b-file; every SEMIPRIME is a\n")
        f.write("# new certified term (prime-winner rows) or a candidate verdict\n")
        f.write("# component (rows at n in {83,84,98,103}).\n\n")
        N = mpz(1); k = 0
        total = 0
        for n, r in todo:
            while k < n:
                N *= sps[k]; k += 1
            f.write(f"# ---- n={n}: winner a_prp={r['a_prp']} ({r['winner_kind']}); "
                    f"{len(r['nasties_below'])} deferred ----\n")
            for m in r["nasties_below"]:
                g = N + m
                f.write(f"# n={n} m={m} ({len(str(g))} digits)\n{g}\n")
                total += 1
        print(f"wrote {total} targets across {len(todo)} indices to {out}")

if __name__ == "__main__":
    main()
