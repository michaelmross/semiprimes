"""
Tier-1 verification for the averaged one-step-Fermat / near-midpoint-Goldbach count.

R(n) = #{ 1 <= x <= sqrt(n) : 2n - x and 2n + x both prime }
T(X) = sum_{n <= X} R(n)

Claimed exact reduction:
T(X) = #{ (p,q) primes, p < q, q - p ≡ 2 (mod 4), (q-p)^2 <= p+q, p+q <= 4X }

(The congruence q-p ≡ 2 mod 4 is *equivalent* to p+q ≡ 0 mod 4 for odd p,q,
so the mod-4 condition lives entirely on the shift.)

Also checks:
  - avg of singular series:  (1/X) sum_{n<=X} S(4n)  ->  2
  - avg over shift class:    (4/(2 sqrt(X)... )) i.e. mean of S(d) over d ≡ 2 mod 4 -> 2
  - T(X) against HL main term sum_{n<=X} S(4n) sqrt(n)/log(2n)^2
"""
from math import isqrt, log

X = 100_000
LIMIT = 4 * X + 2 * isqrt(X) + 10

# ---- sieve ----
sieve = bytearray([1]) * (LIMIT + 1)
sieve[0] = sieve[1] = 0
for i in range(2, isqrt(LIMIT) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))

# smallest prime factor sieve up to X (for S(4n))
spf = list(range(X + 1))
for i in range(2, isqrt(X) + 1):
    if spf[i] == i:
        for j in range(i*i, X + 1, i):
            if spf[j] == j:
                spf[j] = i

def odd_prime_divisors(n):
    s = set()
    while n > 1:
        p = spf[n]
        if p > 2:
            s.add(p)
        while n % p == 0:
            n //= p
    return s

# twin-prime constant C2 (over primes below 10^6; tail correction negligible here)
C2 = 1.0
for p in range(3, 10**6):
    if p < len(sieve) and sieve[p]:
        C2 *= (1.0 - 1.0/((p-1)*(p-1)))
print(f"C2 ~ {C2:.10f}  (true 0.6601618158...)")

def S4n(n):
    """Goldbach singular series S(4n) = 2*C2 * prod_{p | n, p odd} (p-1)/(p-2)."""
    v = 2.0 * C2
    for p in odd_prime_divisors(n):
        v *= (p - 1.0) / (p - 2.0)
    return v

# ---- direct count T(X) ----
T_direct = 0
Rn = [0] * (X + 1)
for n in range(1, X + 1):
    r = 0
    m = 2 * n
    for x in range(1, isqrt(n) + 1):
        if sieve[m - x] and sieve[m + x]:
            r += 1
    Rn[n] = r
    T_direct += r

# ---- pair-count formulation ----
primes = [i for i in range(2, 2 * X + isqrt(X) + 5) if sieve[i]]
T_pairs = 0
d = 2
Dmax = 2 * isqrt(X) + 2
import bisect
while d <= Dmax:
    # p >= ceil((d^2 - d)/2), 2p + d <= 4X  =>  p <= (4X - d)/2
    plo = (d * d - d) // 2      # (d^2 - d) is even, so this is exact
    phi = (4 * X - d) // 2
    if plo < 2:
        plo = 2
    i0 = bisect.bisect_left(primes, plo)
    i1 = bisect.bisect_right(primes, phi)
    for p in primes[i0:i1]:
        if sieve[p + d]:
            # double check exact conditions
            # d^2 <= p + (p+d)  and  p + (p+d) <= 4X  -- guaranteed by range
            T_pairs += 1
    d += 4  # d ≡ 2 (mod 4)

print(f"X = {X}")
print(f"T_direct (sum of R(n))      = {T_direct}")
print(f"T_pairs  (shift-d count)    = {T_pairs}")
print(f"EXACT MATCH: {T_direct == T_pairs}")

# ---- singular series averages ----
sumS = 0.0
main_term = 0.0
for n in range(2, X + 1):
    s = S4n(n)
    sumS += s
    main_term += s * (n ** 0.5) / (log(2 * n) ** 2)

print(f"(1/X) sum S(4n)             = {sumS / X:.6f}   (should -> 2)")

# mean of S over shift class d ≡ 2 mod 4  (S(d) with odd part of d = d/2 odd)
sumSd, cnt = 0.0, 0
d = 2
while d <= Dmax:
    v = 2.0 * C2
    m = d
    while m % 2 == 0:
        m //= 2
    mm = m
    # factor m (small)
    f = 3
    while f * f <= mm:
        if mm % f == 0:
            v *= (f - 1.0) / (f - 2.0)
            while mm % f == 0:
                mm //= f
        f += 2
    if mm > 1:
        v *= (mm - 1.0) / (mm - 2.0)
    sumSd += v
    cnt += 1
    d += 4
print(f"mean S(d), d≡2(4), d<=2√X   = {sumSd / cnt:.6f}   (should -> 2)")

# ---- HL comparison ----
print(f"HL main term sum S(4n)√n/log²(2n) = {main_term:.1f}")
print(f"ratio T(X)/main               = {T_direct / main_term:.5f}")

# cumulative bands
for lo, hi in [(1000, 10000), (10000, 50000), (50000, 100000)]:
    a = sum(Rn[lo+1:hi+1])
    pr = sum(S4n(n) * (n ** 0.5) / (log(2 * n) ** 2) for n in range(lo + 1, hi + 1))
    print(f"band ({lo},{hi}]: actual {a}, HL {pr:.1f}, ratio {a/pr:.5f}")
