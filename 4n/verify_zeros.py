"""
verify_zeros.py -- independent verification of the zero-census artifacts.

Usage:
    python verify_zeros.py [N] [--bpsw K]

Four checks (1-3 always cover the FULL lists; 4 covers n <= N, default 10^6):

  1. checksums    -- every file listed in CHECKSUMS.sha256 that is present
                     matches its recorded SHA-256.
  2. structure    -- header "n"; strictly increasing; all entries >= 2;
                     zeros_J.csv is exactly the intersection of the other two.
  3. soundness    -- every listed zero is confirmed to be a zero by direct
                     enumeration of the admissible x for that n against a
                     freshly built prime sieve (no vectorized slicing).
  4. completeness -- the zero sets are rebuilt from scratch for n <= N by an
                     independently written vectorized recount and must match
                     the CSVs on that range.

  --bpsw K        -- additionally re-test each record plus K random indices
                     per list using sympy's Baillie-PSW primality test, a
                     primality path sharing nothing with the sieve.

Soundness of zeros_J.csv follows from check 2 (J = R intersect U) together
with check 3 on the R and U lists; it is also covered directly by check 4 on
n <= N, which certifies zeros_R.csv and zeros_J.csv end to end whenever N
exceeds their final entries.
"""
import sys, csv, hashlib, random
from math import isqrt
import numpy as np

def sha(fn):
    h = hashlib.sha256()
    with open(fn, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()

def load(fn):
    with open(fn) as f:
        r = csv.reader(f)
        hdr = next(r)
        assert hdr == ["n"], f"{fn}: bad header {hdr}"
        ns = [int(row[0]) for row in r]
    return ns

def main():
    args = [a for a in sys.argv[1:]]
    K = 0
    if "--bpsw" in args:
        i = args.index("--bpsw"); K = int(args[i+1]); del args[i:i+2]
    N = int(float(args[0])) if args else 10**6

    ok = True
    # ---- 1. checksums ----
    try:
        recorded = {}
        for line in open("CHECKSUMS.sha256"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue                       # comment / blank lines (sha256sum -c style)
            parts = line.split()
            if len(parts) >= 2:
                fn = parts[-1].lstrip("*")
                if fn == "CHECKSUMS.sha256":
                    continue                   # a self-entry can never validate
                recorded[fn] = parts[0]
        checked = mismatched = 0
        import os
        for fn, h in recorded.items():
            if os.path.exists(fn):
                checked += 1
                if sha(fn) != h:
                    mismatched += 1; ok = False
                    print(f"[1 checksums] MISMATCH: {fn}")
        print(f"[1 checksums] {checked} files checked against CHECKSUMS.sha256: "
              f"{'OK' if mismatched == 0 else f'{mismatched} FAIL'}")
    except FileNotFoundError:
        print("[1 checksums] CHECKSUMS.sha256 not found: SKIPPED")

    # ---- 2. structure ----
    R, U, J = load("zeros_R.csv"), load("zeros_U.csv"), load("zeros_J.csv")
    s2 = all(ns == sorted(set(ns)) and (not ns or ns[0] >= 2) for ns in (R, U, J))
    s2 &= (set(J) == set(R) & set(U))
    ok &= s2
    print(f"[2 structure ] monotone, >=2, J = R \u2229 U "
          f"({len(R)}+{len(U)} rows, |J| = {len(J)}): {'OK' if s2 else 'FAIL'}")

    # ---- sieve for checks 3-4 ----
    maxn = max(R[-1], U[-1], N)
    LIM = 2 * maxn + 2 * isqrt(4 * maxn) + 10
    sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
    for i in range(2, isqrt(LIM) + 1):
        if sieve[i]: sieve[i*i::i] = False
    pr = sieve  # local alias

    def lower_empty(n):
        m = 2 * n
        for x in range(1, isqrt(n) + 1, 2):
            if pr[m - x] and pr[m + x]:
                return False
        return True

    def upper_empty(n):
        m = 2 * n + 1
        lo = isqrt(3 * n)
        if lo * lo < 3 * n + 1: lo += 1
        x = lo + (lo & 1)
        hi = isqrt(4 * n)
        while x <= hi:
            if pr[m - x] and pr[m + x]:
                return False
            x += 2
        return True

    # ---- 3. soundness (full lists, direct enumeration) ----
    badR = [n for n in R if not lower_empty(n)]
    badU = [n for n in U if not upper_empty(n)]
    s3 = not badR and not badU
    ok &= s3
    print(f"[3 soundness ] all {len(R)} R-zeros and {len(U)} U-zeros confirmed by "
          f"direct enumeration (J follows from 2+3): "
          f"{'OK' if s3 else f'FAIL R:{badR[:3]} U:{badU[:3]}'}")

    # ---- 4. completeness on n <= N (independent vectorized recount) ----
    Rc = np.zeros(N + 1, dtype=np.int32)
    for x in range(1, isqrt(N) + 1, 2):
        idx = np.arange(max(x * x, 2), N + 1)
        Rc[idx[0]:] += (pr[2 * idx - x] & pr[2 * idx + x]).astype(np.int32)
    Uc = np.zeros(N + 1, dtype=np.int32)
    for x in range(2, 2 * isqrt(N) + 2, 2):
        n1 = max((x * x + 3) // 4, 2); n2 = min((x * x - 1) // 3, N)
        if n1 > N: break
        if n1 > n2: continue
        idx = np.arange(n1, n2 + 1)
        Uc[n1:n2 + 1] += (pr[2 * idx + 1 - x] & pr[2 * idx + 1 + x]).astype(np.int32)
    rebuiltR = set((np.nonzero(Rc[2:N+1] == 0)[0] + 2).tolist())
    rebuiltU = set((np.nonzero(Uc[2:N+1] == 0)[0] + 2).tolist())
    s4 = (rebuiltR == {n for n in R if n <= N}) and (rebuiltU == {n for n in U if n <= N})
    ok &= s4
    certR = "end to end" if R[-1] <= N else f"up to {N}"
    certJ = "end to end" if J[-1] <= N else f"up to {N}"
    print(f"[4 complete  ] recount on n <= {N} matches CSVs: {'OK' if s4 else 'FAIL'} "
          f"(certifies zeros_R.csv {certR}, zeros_J.csv {certJ})")

    # ---- optional BPSW cross-machinery spot check ----
    if K:
        from sympy import isprime
        def le_b(n):
            m = 2 * n
            return not any(isprime(m - x) and isprime(m + x) for x in range(1, isqrt(n) + 1, 2))
        def ue_b(n):
            m = 2 * n + 1
            lo = isqrt(3 * n)
            if lo * lo < 3 * n + 1: lo += 1
            xs = range(lo + (lo & 1), isqrt(4 * n) + 1, 2)
            return not any(isprime(m - x) and isprime(m + x) for x in xs)
        rng = random.Random(20260802)
        pickR = {R[-1]} | set(rng.sample(R, min(K, len(R))))
        pickU = {U[-1]} | set(rng.sample(U, min(K, len(U))))
        s5 = all(le_b(n) for n in pickR) and all(ue_b(n) for n in pickU)
        ok &= s5
        print(f"[5 bpsw      ] records + {K}/list random via Baillie-PSW: {'OK' if s5 else 'FAIL'}")

    print("ALL CHECKS PASSED" if ok else "*** FAILURES ABOVE ***")

if __name__ == "__main__":
    main()
