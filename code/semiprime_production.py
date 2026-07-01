#!/usr/bin/env python3
"""
Semiprimes in short intervals: production statistics.

Streams a segmented Omega-sieve over [2, HI). For every segment it computes
Omega(n) exactly via the "remaining-cofactor" method (divide out all primes
p <= sqrt(B); whatever residual > 1 remains is a single large prime). From the
one Omega array per segment we get, for free:
    primes   = (Omega == 1)
    E2       = (Omega == 2)   # exact semiprimes, p^2 included
    E3       = (Omega == 3)

Experiments produced:
  [A] global density: counts at 10^k vs the Landau leading term
  [B] worst-case gap: G_E2(x) at decade checkpoints + location  (Conjecture 1)
  [C] square intervals: E2(n), E3(n), pi(n) per n; detrended sparsity; odd/even
  [D] E3/E2 density ratio (Landau hierarchy check)
  [E] fixed-width windows: mean, var/mean (dispersion)
  [F] random sqrt-scale windows: typical vs worst, fraction empty
  [G] normalized-gap distribution per decade (self-similarity / "fractal")

Resumable: checkpoint pickle every CKPT_EVERY segments. Re-running resumes.

Tune HI / SEG at the top. Defaults are a ~30-60 min job; set HI=10**11 for the
full run (multi-hour). Memory ~ 9 * SEG bytes (int64 cofactor + uint8 omega).
"""

import numpy as np
import math, time, pickle, os, sys

# ----------------------------- configuration -----------------------------
HI         = 10**10          # upper bound (exclusive). Raise to 10**11 for full run.
SEG        = 5 * 10**7       # segment length. Lower if RAM-constrained (rem is int64).
CKPT_EVERY = 5               # checkpoint every this many segments
CKPT_FILE  = "semi_ckpt.pkl"
OUT_FILE   = "semi_results.npz"
MAXG       = 4096            # gap-histogram cap (gaps below 10^11 are << this)
RESUME     = True            # resume from checkpoint if present
# ------------------------------------------------------------------------

# Optional CLI overrides:  python semiprime_production.py [HI] [SEG]
# accepts forms like  1e11  10**11  100000000000
def _parse(tok):
    return int(float(eval(tok, {"__builtins__": {}}))) if tok else None
if len(sys.argv) > 1 and _parse(sys.argv[1]): HI  = _parse(sys.argv[1])
if len(sys.argv) > 2 and _parse(sys.argv[2]): SEG = _parse(sys.argv[2])

A0 = 2
NUM_DEC = int(math.log10(HI)) + 1     # decades 0..NUM_DEC-1
NMAX = math.isqrt(HI) - 1             # last fully-contained square interval index

def base_prime_sieve(limit):
    """All primes <= limit via Eratosthenes (limit = sqrt(HI), small)."""
    s = np.ones(limit + 1, dtype=bool); s[:2] = False
    for i in range(2, int(limit**0.5) + 1):
        if s[i]:
            s[i*i::i] = False
    return np.nonzero(s)[0].astype(np.int64)

def omega_segment(A, B, primes):
    """Exact Omega(n) for n in [A, B). Returns uint8 array of length B-A."""
    L = B - A
    om  = np.zeros(L, dtype=np.uint8)
    rem = np.arange(A, B, dtype=np.int64)
    for p in primes:
        p = int(p)
        if p * p > B:          # all remaining primes exceed sqrt(B): stop
            break
        m0 = ((A + p - 1) // p) * p     # first multiple of p in [A, B)
        if m0 >= B:
            continue
        s = m0 - A
        sub  = om[s::p]                 # strided VIEW -> setitem writes through
        rsub = rem[s::p]
        sub += 1                        # first power: every entry divisible by p
        rsub //= p
        mask = (rsub % p == 0)          # which have p^2 | n
        while mask.any():
            sub[mask] += 1
            rsub[mask] //= p
            mask = (rsub % p == 0)
    om += (rem > 1).astype(np.uint8)    # leftover large prime contributes 1
    return om

def sqrt_index(v):
    """n with n^2 < v <= (n+1)^2, i.e. n = isqrt(v-1), integer-safe (vectorized)."""
    n = np.sqrt((v - 1).astype(np.float64)).astype(np.int64)
    n += ((n + 1) * (n + 1) <= (v - 1))   # nudge up if float undershot
    n -= (n * n > (v - 1))                 # nudge down if float overshot
    return n

# ----------------------------- state -----------------------------
def fresh_state():
    return dict(
        next_A      = A0,
        seg_count   = 0,
        last_sv     = -1,                                   # last semiprime value seen
        cum_semi    = 0, cum_tri = 0, cum_prime = 0,
        dec_semi    = np.zeros(NUM_DEC, dtype=np.int64),    # count <= 10^k
        dec_tri     = np.zeros(NUM_DEC, dtype=np.int64),
        dec_prime   = np.zeros(NUM_DEC, dtype=np.int64),
        gmax        = 0, gmax_loc = 0,                      # running worst gap
        dec_gmax    = np.zeros(NUM_DEC, dtype=np.int64),    # worst gap with end <= 10^k
        E2arr       = np.zeros(NMAX + 2, dtype=np.int64),   # per square interval
        E3arr       = np.zeros(NMAX + 2, dtype=np.int64),
        P1arr       = np.zeros(NMAX + 2, dtype=np.int64),
        gap_hist    = np.zeros((NUM_DEC, MAXG + 1), dtype=np.int64),
    )

def save_ckpt(st):
    tmp = CKPT_FILE + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(st, f, protocol=4)
    os.replace(tmp, CKPT_FILE)

# ----------------------------- main pass -----------------------------
def run():
    t0 = time.time()
    primes = base_prime_sieve(math.isqrt(HI) + 1)
    print(f"base primes <= {math.isqrt(HI)+1}: {len(primes)}   ({time.time()-t0:.1f}s)")

    st = fresh_state()
    if RESUME and os.path.exists(CKPT_FILE):
        with open(CKPT_FILE, "rb") as f:
            st = pickle.load(f)
        print(f"resumed from checkpoint at A={st['next_A']}  ({st['seg_count']} segs done)")

    decade_bounds = [10**k for k in range(NUM_DEC)]

    while st["next_A"] < HI:
        A = st["next_A"]; B = min(A + SEG, HI)
        om = omega_segment(A, B, primes)
        v  = np.arange(A, B, dtype=np.int64)

        prime_mask = (om == 1); semi_mask = (om == 2); tri_mask = (om == 3)
        sv = v[semi_mask]

        # --- [B][G] gap walk across the segment boundary ---
        if st["last_sv"] >= 0 and sv.size:
            full = np.concatenate(([st["last_sv"]], sv))
        else:
            full = sv
        if full.size >= 2:
            gaps = np.diff(full)
            ends = full[1:]
            gi = gaps.argmax()
            if gaps[gi] > st["gmax"]:
                st["gmax"] = int(gaps[gi]); st["gmax_loc"] = int(ends[gi])
            # per-decade histogram + per-decade running max, keyed by end value
            dk = np.floor(np.log10(ends)).astype(np.int64)
            gclip = np.minimum(gaps, MAXG)
            for k in np.unique(dk):
                m = (dk == k)
                st["gap_hist"][k] += np.bincount(gclip[m], minlength=MAXG + 1)
                gmk = int(gaps[m].max())
                if gmk > st["dec_gmax"][k]:
                    st["dec_gmax"][k] = gmk
        if sv.size:
            st["last_sv"] = int(sv[-1])

        # --- [C] square intervals (segment-safe via unique n attribution) ---
        for mask, arr in ((semi_mask, "E2arr"), (tri_mask, "E3arr"), (prime_mask, "P1arr")):
            vv = v[mask]
            if vv.size:
                n = sqrt_index(vv)
                n = n[n <= NMAX]
                if n.size:
                    st[arr] += np.bincount(n, minlength=NMAX + 2)

        # --- [A] decade density snapshots (exact, handles straddling) ---
        for k, D in enumerate(decade_bounds):
            if A <= D < B:
                st["dec_semi"][k]  = st["cum_semi"]  + int((sv <= D).sum())
                st["dec_tri"][k]   = st["cum_tri"]   + int((v[tri_mask]   <= D).sum())
                st["dec_prime"][k] = st["cum_prime"] + int((v[prime_mask] <= D).sum())
        st["cum_semi"]  += int(semi_mask.sum())
        st["cum_tri"]   += int(tri_mask.sum())
        st["cum_prime"] += int(prime_mask.sum())

        st["next_A"] = B; st["seg_count"] += 1
        if st["seg_count"] % CKPT_EVERY == 0:
            save_ckpt(st)
        frac = (B - A0) / (HI - A0)
        el = time.time() - t0
        print(f"  [{B:>13d}] {100*frac:5.1f}%  semi={st['cum_semi']:>12d}  "
              f"gmax={st['gmax']}@{st['gmax_loc']}  {el:6.0f}s  eta {el/frac-el:6.0f}s",
              flush=True)

    save_ckpt(st)
    print(f"main pass done in {time.time()-t0:.0f}s")
    return st, primes

# ----------------------------- reporting -----------------------------
def report(st, primes):
    print("\n========== [A] global density vs Landau x*lglg/lg ==========")
    for k in range(4, NUM_DEC):
        x = 10**k; c = int(st["dec_semi"][k])
        if c == 0: continue
        pred = x*math.log(math.log(x))/math.log(x)
        print(f"  <=1e{k}: count={c:>13d}  Landau={pred:>15.0f}  ratio={c/pred:.4f}")

    print("\n========== [B] worst-case gap G_E2(x) (Conjecture 1) ==========")
    # dec_gmax[j] is the worst gap whose END value lies in [10^j, 10^{j+1}).
    # Cumulative G(10^k) = worst gap with end <= 10^k = max over decades 0..k-1.
    Gcum = np.maximum.accumulate(np.concatenate(([0], st["dec_gmax"])))  # Gcum[k]=G(10^k)
    for k in range(4, NUM_DEC):
        x = 10**k; G = int(Gcum[k])
        if G == 0: continue
        print(f"  x=1e{k}: G={G:>4d}  sqrt={math.sqrt(x):>12.0f}  (lgx)^2={math.log(x)**2:7.1f}"
              f"  G/sqrt={G/math.sqrt(x):.3e}")
    print(f"  overall worst gap below {HI}: G={st['gmax']} at {st['gmax_loc']}")
    st["Gcum"] = Gcum

    print("\n========== [C] square intervals: E2(n), detrended, odd/even ==========")
    E2 = st["E2arr"]; P1 = st["P1arr"]; E3 = st["E3arr"]
    n = np.arange(NMAX + 2, dtype=np.float64)
    x = n*(n+1)                           # representative height of (n^2,(n+1)^2]
    with np.errstate(divide="ignore", invalid="ignore"):
        dens = np.log(np.log(x))/np.log(x)
        Ehat = (2*n + 1)*dens             # corrected leading prediction
        r = np.where(Ehat > 0, E2/Ehat, np.nan)   # detrended ratio
    for k in range(1, int(math.log10(NMAX)) + 1):
        lo, hi = 10**k, min(10**(k+1), NMAX + 1)
        if lo >= hi: break
        idx = np.arange(lo, hi)
        e2 = E2[idx]; rr = r[idx]; pr = P1[idx]
        odd = idx % 2 == 1
        jmin = idx[np.nanargmin(rr)]
        print(f"  n in [{lo},{hi}): E2 min={e2.min()} mean={e2.mean():.1f} | "
              f"r mean={np.nanmean(rr):.3f} (odd {np.nanmean(rr[odd]):.3f}, "
              f"even {np.nanmean(rr[~odd]):.3f}) | most-deficient n={jmin} "
              f"({'odd' if jmin%2 else 'even'}) r={r[jmin]:.3f} E2={E2[jmin]} | "
              f"pi mean odd {pr[odd].mean():.1f} even {pr[~odd].mean():.1f}")
    print(f"  global E2(n) min over n in [2,{NMAX}] = {E2[2:NMAX+1].min()} "
          f"(never 0 supports Conjecture 2)")

    print("\n========== [D] E3/E2 density ratio (Landau hierarchy) ==========")
    for k in range(4, NUM_DEC):
        if st["dec_semi"][k] == 0: continue
        x = 10**k
        print(f"  x=1e{k}: E3/E2={st['dec_tri'][k]/st['dec_semi'][k]:.3f}  "
              f"(predict lglg/2={math.log(math.log(x))/2:.3f})")

    print("\n========== [E] fixed-width windows: dispersion ==========")
    H = 4000
    for x0 in sorted({10**6, 10**8, HI//2}):
        if x0 + 300*H >= HI: continue
        om = omega_segment(x0, x0 + 300*H, primes)
        semi = (om == 2)
        cnts = np.array([int(semi[i*H:(i+1)*H].sum()) for i in range(300)])
        pred = H*math.log(math.log(x0))/math.log(x0)
        print(f"  x0={x0:>13d}: mean={cnts.mean():.2f} pred={pred:.2f} "
              f"var/mean={cnts.var()/cnts.mean():.2f} min={cnts.min()}")

    print("\n========== [F] random sqrt-scale windows: typical vs worst ==========")
    x0 = HI // 2; W = int(2*math.sqrt(x0))
    BAND = min(10**8, HI - x0 - W - 1)        # bounded band: caps memory & time
    om = omega_segment(x0, x0 + BAND, primes); semi = (om == 2)
    cs = np.concatenate(([0], np.cumsum(semi, dtype=np.int64)))
    rng = np.random.default_rng(0)
    offs = rng.integers(0, BAND - W, size=2000)
    cnts = cs[offs + W] - cs[offs]
    print(f"  x0={x0} W={W} band={BAND}: mean={cnts.mean():.1f} min={int(cnts.min())} "
          f"max={int(cnts.max())} frac_empty={(cnts == 0).mean():.5f}")

    print("\n========== [G] normalized-gap distribution per decade ==========")
    bins = np.arange(MAXG + 1)
    for k in range(2, NUM_DEC):
        h = st["gap_hist"][k]; tot = h.sum()
        if tot == 0: continue
        mean = (bins*h).sum()/tot
        var  = ((bins - mean)**2 * h).sum()/tot
        gmax = int(np.nonzero(h)[0].max())
        print(f"  1e{k}-1e{k+1}: mean_gap={mean:.3f}  CV={math.sqrt(var)/mean:.3f}  "
              f"max_gap={gmax}  max/mean={gmax/mean:.2f}")

    # ---- validation (self-contained at any height) ----
    print("\n========== [CHECK] validation ==========")
    # (i) reference counts proven by two independent methods in the pilot
    REF = {4: 2625, 5: 23378, 6: 210035, 7: 1904324, 8: 17427258}
    for k, c in REF.items():
        if k < NUM_DEC and int(st["dec_semi"][k]) > 0:
            got = int(st["dec_semi"][k])
            print(f"  count<=1e{k}: {got:>13d}  ref {c:>13d}  match={got == c}")
    # (ii) brute-force Omega spot-check near the top of the range
    try:
        from sympy import factorint
        A = (HI // 2) | 1; B = A + 2000
        om = omega_segment(A, B, primes)
        bad = sum(int(om[i]) != sum(factorint(A + i).values()) for i in range(B - A))
        print(f"  brute-force Omega on [{A},{B}): {'PASS' if bad == 0 else f'FAIL ({bad})'}")
    except ImportError:
        print("  (sympy not installed; skipping brute-force spot-check)")

    np.savez(OUT_FILE, E2arr=st["E2arr"], E3arr=st["E3arr"], P1arr=st["P1arr"],
             gap_hist=st["gap_hist"], dec_semi=st["dec_semi"], dec_tri=st["dec_tri"],
             dec_prime=st["dec_prime"], dec_gmax=st["dec_gmax"], Gcum=st.get("Gcum"),
             gmax=st["gmax"], gmax_loc=st["gmax_loc"], HI=HI)
    print(f"\nsaved arrays to {OUT_FILE}")

if __name__ == "__main__":
    state, prms = run()
    report(state, prms)
