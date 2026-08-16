#!/usr/bin/env python3
"""
Scanner for A226525 (Fortunate semiprimes): a(n) = least m > 1 such that
m + N_n is a semiprime (Omega = 2), where N_n = A112141(n) = product of
the first n semiprimes.

Goal: hunt for a COMPOSITE term a(n), falsifying the OEIS comment
"we conjecture that all terms are prime".

Structural facts used (all elementary):
  (S1) A prime p divides N_n  iff  2p <= sp(n) (i.e. p <= sp(n)/2), since
       2p is then one of the first n semiprimes; small p divide N_n to
       high multiplicity.
  (S2) For p <= sp(n)/2:  p | m + N_n  <=>  p | m.  So the small-prime
       structure of the shifted number is inherited from the offset m.
  (S3) Consequently, in the race range m = O(log^2 N) << (sp(n)/2)^2, any
       m with NO prime factor <= sp(n)/2 is 1 or prime ("rough m").
  (S4) A candidate g = m + N_n with a stripped small-prime count >= 1
       whose remaining cofactor is composite has Omega >= 3: resolved
       WITHOUT factoring.  The only unresolvable-at-PRP-level candidates
       ("nasty") are rough PRIME m whose g is composite with no factor
       below the trial bound.  Nasty m are always prime, so an unresolved
       nasty can never hide a composite counterexample -- it can only
       relocate a(n) to a smaller PRIME value.  Certification of a
       composite hit therefore requires ECM only on the nasties BELOW it.

Per candidate m the classification is:
  g = m + N_n
  1. strip: for each prime p | m with p <= sp(n)/2, divide out p^{v_p(g)}
     exactly (v_p(g) computed directly; by (S2) no other p <= sp(n)/2
     divides g).
  2. strip further: gcd with a primorial-style product of the primes in
     (sp(n)/2, TRIAL_B], factor the gcd, divide out.
  3. cofactor h: h == 1        -> Omega = stripped
                 h PRP         -> Omega = stripped + 1
                 h composite   -> Omega >= stripped + 2
                                  (resolved non-semiprime if stripped>=1,
                                   else NASTY: needs ECM, m is prime)
Semiprime <=> Omega == 2 (with h-part at PRP confidence).

Output per n: a_prp(n); whether the winner m is composite (HIT);
nasty m below the winner (certification burden); the threshold 2*q0
(q0 = least prime not dividing N_n) marking where composite candidates
enter the race; and all composite candidates that were alive in the race.
"""

import sys, time, json
import gmpy2
from gmpy2 import mpz, is_prime, gcd

TRIAL_B = 10**6          # trial-division bound for medium primes
MR_EXTRA = 24            # extra Miller-Rabin rounds on top of gmpy2 BPSW-style test

# ---------------------------------------------------------------- sieve
def primes_upto(B):
    sieve = bytearray([1]) * (B + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(B**0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(B + 1) if sieve[i]]

print("sieving primes to", TRIAL_B, "...", file=sys.stderr)
PRIMES = primes_upto(TRIAL_B)
PRIME_SET = set(PRIMES)

def is_prp(x):
    return is_prime(x, MR_EXTRA)      # BPSW + extra MR rounds

# ------------------------------------------------- semiprimes / semiprimorials
def omega_small(x):
    """Exact Omega for small x (trial division; x < TRIAL_B^2 assumed)."""
    cnt, t = 0, x
    for p in PRIMES:
        if p * p > t:
            break
        while t % p == 0:
            t //= p
            cnt += 1
    if t > 1:
        cnt += 1
    return cnt

def gen_semiprimes(limit_count):
    out, x = [], 3
    while len(out) < limit_count:
        x += 1
        if omega_small(x) == 2:
            out.append(x)
    return out

# ---------------------------------------------------------------- per-n scan
def scan_n(n, sps, prod_cache, report_all_composites=True, race_cap_mult=8):
    """
    Scan index n.  Returns a dict of results at PRP level.
    """
    sp_n = sps[n - 1]
    N = prod_cache[n]
    half = sp_n // 2                      # primes <= half divide N  (S1)
    # least prime NOT dividing N:
    q0 = next(p for p in PRIMES if p > half)
    lnN = gmpy2.log(N)
    race_cap = int(race_cap_mult * lnN) + 50   # generous stopping guard

    # product of primes in (half, TRIAL_B]  (for gcd stripping); cache by 'half'
    key = ("midprod", half)
    if key not in prod_cache:
        P = mpz(1)
        for p in PRIMES:
            if p > half:
                P *= p
        prod_cache[key] = P
    MIDP = prod_cache[key]

    nasties = []            # unresolved prime m (need ECM to certify)
    alive_composites = []   # composite m whose g survived to the PRP stage
    winner = None           # (m, kind) kind in {"prime-m","composite-m"}

    m = 1
    while m < race_cap:
        m += 1
        g = N + m
        stripped = 0
        h = g
        # --- step 1: strip primes p | m, p <= half (exact, by (S2))
        t = m
        small_ps = []
        for p in PRIMES:
            if p > half or p * p > t:
                break
            if t % p == 0:
                small_ps.append(p)
                while t % p == 0:
                    t //= p
        if 1 < t <= half:
            small_ps.append(t)            # t is a prime factor of m, <= half
            t = 1
        # (any remaining t is a prime factor of m exceeding half: does NOT divide g)
        for p in small_ps:
            while h % p == 0:
                h //= p
                stripped += 1
                if stripped > 2:
                    break
            if stripped > 2:
                break
        if stripped > 2:
            continue                      # Omega >= 3, resolved
        # --- step 2: strip medium primes in (half, TRIAL_B]
        if h > 1:
            d = gcd(h, MIDP)
            while d > 1:
                # factor d over the medium primes (d is smallish & squarefree-ish)
                dd = int(d)
                for p in PRIMES:
                    if p <= half:
                        continue
                    if p * p > dd and dd > 1:
                        # dd itself is a medium prime
                        p = dd
                    if dd % p == 0:
                        while h % p == 0:
                            h //= p
                            stripped += 1
                        while dd % p == 0:
                            dd //= p
                    if dd == 1 or stripped > 2:
                        break
                if stripped > 2:
                    break
                d = gcd(h, MIDP)
        if stripped > 2:
            continue
        # --- step 3: classify cofactor
        if h == 1:
            om = stripped
            resolved = True
        elif is_prp(h):
            om = stripped + 1
            resolved = True
        else:
            om = None
            resolved = (stripped >= 1)    # composite h: Omega >= stripped+2 >= 3
        m_is_composite = omega_small(m) >= 2

        if resolved:
            if om == 2:
                winner = (m, "composite-m" if m_is_composite else "prime-m")
                break
            # else non-semiprime, continue
        else:
            # nasty: rough m (necessarily prime by (S3)), g composite & rough
            nasties.append(m)
        if m_is_composite and not resolved:
            pass  # cannot happen: composite m always resolves (stripped>=1)
        if m_is_composite and resolved and om is not None and om != 2:
            if om == 3 and report_all_composites:
                pass
        # track composite candidates that reached the PRP stage (stripped==1, h tested)
        if m_is_composite and stripped == 1 and h > 1:
            alive_composites.append((m, "won" if (winner and winner[0] == m) else "lost"))

    return {
        "n": n, "sp_n": sp_n, "digits_N": len(str(N)), "lnN": float(lnN),
        "q0": q0, "two_q0": 2 * q0,
        "a_prp": winner[0] if winner else None,
        "winner_kind": winner[1] if winner else "NONE-IN-CAP",
        "nasties_below": [x for x in nasties if winner and x < winner[0]],
        "alive_composites": alive_composites,
    }

# ---------------------------------------------------------------- driver
def main():
    n_lo = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n_hi = int(sys.argv[2]) if len(sys.argv) > 2 else 65

    sps = gen_semiprimes(n_hi + 5)
    prod_cache = {0: mpz(1)}
    for i in range(1, n_hi + 1):
        prod_cache[i] = prod_cache[i - 1] * sps[i - 1]

    OEIS = [2,2,2,5,2,3,2,7,3,19,11,3,23,5,61,29,31,3,29,31,13,19,5,7,23,47,
            3,53,47,19,13,7,41,53,2,43,7,103,2,61,59,71,17,59,79,43,167,71,
            97,7,151,37,103,83,127,103,11,53,29,7,67,83,151,107,37]

    results = []
    t0 = time.time()
    for n in range(n_lo, n_hi + 1):
        r = scan_n(n, sps, prod_cache)
        results.append(r)
        ok = ""
        if n <= 65:
            ok = "OK" if r["a_prp"] == OEIS[n - 1] else f"MISMATCH(oeis={OEIS[n-1]})"
        comp_flag = "  <<< COMPOSITE WINNER" if r["winner_kind"] == "composite-m" else ""
        print(f"n={n:3d} sp={r['sp_n']:4d} dig(N)={r['digits_N']:3d} "
              f"2q0={r['two_q0']:4d} a_prp={r['a_prp']} [{r['winner_kind']}] "
              f"nasties<{r['a_prp']}: {len(r['nasties_below'])} "
              f"aliveC={len(r['alive_composites'])} {ok}{comp_flag}", flush=True)
        if r["alive_composites"]:
            print(f"      composite candidates in race: {r['alive_composites']}", flush=True)
        if r["nasties_below"]:
            print(f"      nasty (unresolved, prime m): {r['nasties_below']}", flush=True)
    print(f"[{time.time()-t0:.1f}s total]", file=sys.stderr)
    with open(f"scan_{n_lo}_{n_hi}.json", "w") as f:
        json.dump(results, f, indent=1)

if __name__ == "__main__":
    main()
