"""
census_final.py -- one-step-Fermat semiprimes in J_n = [4n^2-n, 4n^2+n]: full census.

Usage:  python census_final.py [X]        (default X = 1_000_000)

Counts, for 1 <= n <= X:
  R(n) = #{x odd,  1 <= x <= sqrt(n)          : 2n-x,   2n+x   both prime}   (lower half)
  U(n) = #{x even, 3n+1 <= x^2 <= 4n          : 2n+1-x, 2n+1+x both prime}   (upper half)

Reports (all bands scale with X):
  - exact-identity regression checks (direct interval scan, n <= 300; lower-half sum at 1e5)
  - HL band ratios for both halves; upper/lower aggregate vs (2 - sqrt 3)
  - zero censuses for R, U, R+U with per-band observed vs Poisson-predicted counts
  - singular-series stratification over (X/2, X], correctly labelled
  - dumps zero lists to CSV with SHA256 checksums (zeros_R.csv, zeros_U.csv, zeros_J.csv)

Memory at X = 1e8: ~400MB prime sieve + ~800MB SPF sieve + ~400MB counts. Fits in 4GB.
Corrects two defects of item3_census.py: hardcoded band limits (bands beyond 1e6 were
silently unreported at larger X) and a mislabelled stratification header (the range is
(X/2, X], which the header now states).
"""
import sys, hashlib
import numpy as np
from math import isqrt, log, sqrt

X = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
LIM = 4 * X + 2 * isqrt(4 * X) + 10

# ---------------- sieves ----------------
sieve = np.ones(LIM + 1, dtype=bool); sieve[:2] = False
for i in range(2, isqrt(LIM) + 1):
    if sieve[i]:
        sieve[i * i::i] = False

M = 2 * X + 1                       # SPF to 2X+1 (factor n and 2n+1)
spf = np.arange(M + 1, dtype=np.int32)
for i in range(2, isqrt(M) + 1):
    mask = spf[i * i::i] == np.arange(i * i, M + 1, i, dtype=np.int32)
    spf[i * i::i] = np.where(mask, i, spf[i * i::i])

C2 = 1.0
for p in range(3, min(10**6, LIM)):
    if sieve[p]:
        C2 *= 1 - 1 / ((p - 1) * (p - 1))

def S_of(n):
    """Goldbach singular series S(2m) computed from the odd primes of the argument n."""
    v = 2 * C2
    while n > 1:
        p = int(spf[n])
        if p > 2:
            v *= (p - 1) / (p - 2)
        while n % p == 0:
            n //= p
    return v

def omega_odd(n):
    c = 0
    while n > 1:
        p = int(spf[n]); c += (p > 2)
        while n % p == 0:
            n //= p
    return c

# ---------------- counts ----------------
R = np.zeros(X + 1, dtype=np.int16)
for x in range(1, isqrt(X) + 1, 2):
    n0 = max(x * x, 1)
    if n0 > X:
        break
    R[n0:X + 1] += (sieve[2 * n0 - x:2 * X - x + 1:2] & sieve[2 * n0 + x:2 * X + x + 1:2])

U = np.zeros(X + 1, dtype=np.int16)
for x in range(2, 2 * isqrt(X) + 2, 2):
    n1 = (x * x + 3) // 4; n2 = min((x * x - 1) // 3, X)
    if n1 > X:
        break
    if n1 > n2:
        continue
    U[n1:n2 + 1] += (sieve[2 * n1 + 1 - x:2 * n2 + 1 - x + 1:2] & sieve[2 * n1 + 1 + x:2 * n2 + 1 + x + 1:2])

# ---------------- regression checks ----------------
if X >= 100_000:
    s = int(R[:100_001].sum())
    print(f"regression: sum R(n), n<=1e5 = {s} (expect 319013) {'OK' if s == 319013 else 'FAIL'}")
NN = 4 * 300 * 300 + 300 + 2
spf2 = np.arange(NN, dtype=np.int64)
for i in range(2, isqrt(NN) + 1):
    mask = spf2[i * i::i] == np.arange(i * i, NN, i, dtype=np.int64)
    spf2[i * i::i] = np.where(mask, i, spf2[i * i::i])
def is_sp(N):
    p = int(spf2[N]); q = N // p
    return q > 1 and int(spf2[q]) == q
ok = True
for n in range(2, 301):
    lo = sum(1 for N in range(4*n*n - n, 4*n*n + 1)
             if is_sp(N) and isqrt((isqrt(N-1)+1)**2 - N)**2 == (isqrt(N-1)+1)**2 - N)
    up = sum(1 for N in range(4*n*n + 1, 4*n*n + n + 1)
             if is_sp(N) and isqrt((isqrt(N-1)+1)**2 - N)**2 == (isqrt(N-1)+1)**2 - N)
    ok &= (lo == R[n] and up == U[n])
print(f"regression: direct interval scan n<=300 exact match: {'OK' if ok else 'FAIL'}")

# ---------------- per-n HL predictions (chunked; no full float arrays) ----------------
def band_edges(X):
    e = [10**3]
    while e[-1] * 10 <= X:
        e.append(e[-1] * 10)
    if e[-1] != X:
        e.append(X)
    return [(e[i], e[i + 1]) for i in range(len(e) - 1)]

CH = 10**6
def band_stats(lo, hi):
    aR = int(R[lo + 1:hi + 1].sum()); aU = int(U[lo + 1:hi + 1].sum())
    pR = pU = eZR = eZJ = 0.0
    oZR = oZJ = 0
    for a in range(lo + 1, hi + 1, CH):
        b = min(a + CH - 1, hi)
        ns = np.arange(a, b + 1, dtype=float)
        L2 = np.log(2 * ns) ** 2
        S1 = np.array([S_of(n) for n in range(a, b + 1)])
        S2 = np.array([S_of(2 * n + 1) for n in range(a, b + 1)])
        lamR = S1 * np.sqrt(ns) / L2
        lamU = S2 * np.maximum(0, np.sqrt(4 * ns) - np.sqrt(3 * ns + 1)) / L2
        pR += lamR.sum(); pU += lamU.sum()
        eZR += np.exp(-lamR).sum(); eZJ += np.exp(-(lamR + lamU)).sum()
        oZR += int((R[a:b + 1] == 0).sum())
        oZJ += int(((R[a:b + 1] + U[a:b + 1]) == 0).sum())
    return aR, pR, aU, pU, oZR, eZR, oZJ, eZJ

print(f"\nX = {X}")
print("band                 lower act/HL        upper act/HL        R-zeros obs/Poisson   J-zeros obs/Poisson")
for lo, hi in band_edges(X):
    aR, pR, aU, pU, oZR, eZR, oZJ, eZJ = band_stats(lo, hi)
    print(f"({lo:.0e},{hi:.0e}]   {aR}/{pR:.0f} = {aR/pR:.4f}   {aU}/{pU:.0f} = {aU/pU:.4f}   {oZR}/{eZR:.1f}   {oZJ}/{eZJ:.1f}")
hl = max(X // 10, 1000)
ur = U[hl + 1:].sum() / max(R[hl + 1:].sum(), 1)
print(f"upper/lower aggregate on ({hl:.0e},{X:.0e}]: {ur:.4f}  (2-sqrt3 = {2 - sqrt(3):.4f})")

# ---------------- zero censuses + dumps ----------------
zR = np.nonzero(R[2:X + 1] == 0)[0] + 2
zU = np.nonzero(U[2:X + 1] == 0)[0] + 2
zJ = np.nonzero((R[2:X + 1] + U[2:X + 1]) == 0)[0] + 2
print(f"\nR-zeros: {len(zR)}, last {zR[-1] if len(zR) else '-'}")
print(f"U-zeros: {len(zU)}, last {zU[-1] if len(zU) else '-'}")
print(f"J-zeros (R+U=0): {len(zJ)}, last {zJ[-1] if len(zJ) else '-'}")
for name, arr in [("zeros_R.csv", zR), ("zeros_U.csv", zU), ("zeros_J.csv", zJ)]:
    data = "n\n" + "\n".join(map(str, arr.tolist())) + "\n"
    with open(name, "w", newline="\n") as f:
        f.write(data)
    h = hashlib.sha256(data.encode()).hexdigest()
    print(f"  wrote {name} ({len(arr)} rows), sha256 {h[:16]}...")

# ---------------- stratification over (X/2, X], correctly labelled ----------------
print(f"\nstratification over ({X // 2:.0e},{X:.0e}] by omega = #distinct odd prime divisors of n:")
print("omega     #n           mean r        mean S(4n)    ratio")
acc = {}
for a in range(X // 2 + 1, X + 1, CH):
    b = min(a + CH - 1, X)
    ns = np.arange(a, b + 1, dtype=float)
    rn = R[a:b + 1] * np.log(2 * ns) ** 2 / np.sqrt(ns)
    Sv = np.array([S_of(n) for n in range(a, b + 1)])
    om = np.array([omega_odd(n) for n in range(a, b + 1)], dtype=np.int8)
    pr = sieve[a:b + 1]
    m105 = (np.arange(a, b + 1) % 105 == 0)
    for key, msk in [*((w, om == w) for w in range(0, 8)), ("prime", pr), ("0mod105", m105)]:
        if key not in acc:
            acc[key] = [0, 0.0, 0.0]
        acc[key][0] += int(msk.sum()); acc[key][1] += float(rn[msk].sum()); acc[key][2] += float(Sv[msk].sum())
for key in list(range(0, 8)) + ["prime", "0mod105"]:
    c, sr, ss = acc.get(key, [0, 0, 0])
    if c < 50:
        continue
    print(f"  {str(key):8s} {c:10d}    {sr/c:8.4f}      {ss/c:8.4f}      {(sr/c)/(ss/c):.4f}")
