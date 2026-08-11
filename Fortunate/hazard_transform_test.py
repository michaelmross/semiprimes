#!/usr/bin/env python3
"""
hazard_transform_test.py -- Exp(1) hazard-transform test for Model 1 of
"The Fortunate Semiprimes" (A226525).

Idea. Model 1 defines, for each index n, a discrete hazard h_n(m) on offsets
m >= 2 and predicts survival P(a(n) > t) ~ exp(-H_n(t)) with
H_n(t) = sum_{2 <= m <= t} h_n(m). If the model is correct, the transformed
values Z_n = H_n(a(n)) are (approximately) i.i.d. Exp(1) across n. This
script computes Z_n for every exactly determined term, tests the sample
against Exp(1), estimates the global hazard deflation as 1/mean(Z), and
localizes any miscalibration by channel (inherited vs rough) and by whether
the race ended below or above y.

Hazards (Model 1, Section 5 of the paper), with y = sp(n)/2,
L = log Q_n, u = L/log y, rho = e^gamma * log y:

  channel I   : m prime, m <= y            h = rho/L
  channel II  : m prime, y < m < q0^2      h = (rho/L) * log(u-1)
  channel III : m = p*Q, p <= y < Q both prime   h = rho/L
                m = p^2, sp(n)/3 < p <= sp(n)/2  h = rho/L
  otherwise   : h = 0   (>= 2 small prime factors, or excluded shapes of
                         Lemma 4's proof)

Discreteness. The transform Z = H_n(a(n)) is exactly Exp(1) only for a
continuous hazard. With discrete per-offset hazards h the winning offset
contributes its full h, biasing Z upward by ~h/2. We therefore report both
Z_full (hazard through a(n) inclusive) and the midpoint-corrected
Z_mid = H_n(a(n)-) + h_n(a(n))/2, and use Z_mid for the tests. At these
heights h is 0.02-0.10 per offset, so the correction is small but free.

Input. A CSV with columns n,a,status (header optional). Only rows whose
status is one of {proved, probable, exact} are used; anything else
(provisional, frontrunner, open, ...) is skipped, which enforces the
censoring discipline of Remark 5: provisional front-runners are upper
bounds resolved by preferential ECM and must not enter the calibration.

Usage:
  python3 hazard_transform_test.py values.csv [--nmin 66] [--nmax 93]
                                             [--boot 20000] [--seed 1]

The --nmin/--nmax filters reproduce the paper's calibration window; omit
them to use every exact term in the file. --boot sets the parametric
bootstrap replicates for the estimated-scale (Lilliefors-type) KS p-value.

Dependencies: sympy, scipy, numpy.
"""

import argparse
import csv
import math
import sys

import numpy as np
from scipy import stats
from sympy import isprime, factorint, nextprime

EULER_GAMMA = 0.5772156649015328606


# ----------------------------------------------------------------------
# Semiprimes and per-index model parameters
# ----------------------------------------------------------------------

def semiprimes_up_to_index(nmax):
    """Return the list [q_1, ..., q_nmax] of the first nmax semiprimes."""
    sps = []
    k = 3
    while len(sps) < nmax:
        k += 1
        f = factorint(k)
        if sum(f.values()) == 2:
            sps.append(k)
    return sps


class IndexModel:
    """Model-1 parameters and hazards for a single index n."""

    def __init__(self, n, sps):
        self.n = n
        self.sp = sps[n - 1]
        self.y = self.sp / 2.0
        self.L = sum(math.log(q) for q in sps[:n])  # log Q_n
        self.logy = math.log(self.y)
        self.u = self.L / self.logy
        self.rho = math.exp(EULER_GAMMA) * self.logy
        self.h_base = self.rho / self.L                    # channels I, III
        self.h_rough = self.h_base * math.log(self.u - 1)  # channel II
        self.q0 = int(nextprime(math.floor(self.y)))
        self.q0sq = self.q0 * self.q0

    def hazard(self, m):
        """Return (h, channel) for offset m; channel in {'I','II','III',None}."""
        if m < 2:
            return 0.0, None
        if isprime(m):
            if m <= self.y:
                return self.h_base, "I"
            if m < self.q0sq:
                return self.h_rough, "II"
            return 0.0, None  # beyond the modeled window
        # composite: admitted shapes of Lemma 4 only
        f = factorint(m)
        if len(f) == 2:
            (p1, e1), (p2, e2) = sorted(f.items())
            if e1 == 1 and e2 == 1 and p1 <= self.y < p2 and m < self.q0sq:
                return self.h_base, "III"
        elif len(f) == 1:
            (p, e), = f.items()
            if e == 2 and self.sp / 3.0 < p <= self.sp / 2.0:
                return self.h_base, "III"
        return 0.0, None

    def transform(self, a):
        """Cumulative hazard at a(n) = a, with channel decomposition.

        Returns dict with Z_full, Z_mid, per-channel totals, H at y,
        and the winning offset's own hazard.
        """
        H = 0.0
        by_channel = {"I": 0.0, "II": 0.0, "III": 0.0}
        H_below_y = 0.0
        h_last = 0.0
        for m in range(2, a + 1):
            h, ch = self.hazard(m)
            if h == 0.0:
                continue
            H += h
            by_channel[ch] += h
            if m <= self.y:
                H_below_y += h
            if m == a:
                h_last = h
        if h_last == 0.0:
            raise ValueError(
                f"n={self.n}: reported term a={a} has zero model hazard "
                f"(not an admitted offset). Check the value."
            )
        return {
            "Z_full": H,
            "Z_mid": H - 0.5 * h_last,
            "h_last": h_last,
            "H_I": by_channel["I"],
            "H_II": by_channel["II"],
            "H_III": by_channel["III"],
            "H_below_y": H_below_y,
        }


# ----------------------------------------------------------------------
# Input
# ----------------------------------------------------------------------

EXACT_STATUSES = {"proved", "probable", "exact"}


def read_values(path):
    rows = []
    with open(path, newline="") as fh:
        for raw in csv.reader(fh):
            if not raw or raw[0].strip().startswith("#"):
                continue
            if raw[0].strip().lower() in {"n", "index"}:
                continue  # header
            n = int(raw[0])
            a = int(raw[1])
            status = raw[2].strip().lower() if len(raw) > 2 else "exact"
            rows.append((n, a, status))
    return rows


# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

def ks_exp1(z):
    """KS test against Exp(1), fully specified null."""
    return stats.ks_1samp(z, stats.expon.cdf)


def ks_exp_fitted_bootstrap(z, nboot, rng):
    """KS against Exp(theta) with theta = mean(z) estimated from the data.

    The null distribution of the statistic is not the standard KS law when
    the scale is estimated (Lilliefors), so the p-value is obtained by
    parametric bootstrap: simulate Exp(1) samples of the same size,
    re-estimate the scale each time, and compare statistics.
    """
    z = np.asarray(z, dtype=float)
    theta = z.mean()
    d_obs = stats.ks_1samp(z / theta, stats.expon.cdf).statistic
    n = len(z)
    count = 0
    for _ in range(nboot):
        sim = rng.exponential(1.0, size=n)
        d = stats.ks_1samp(sim / sim.mean(), stats.expon.cdf).statistic
        if d >= d_obs:
            count += 1
    return theta, d_obs, (count + 1) / (nboot + 1)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("values", help="CSV of n,a,status")
    ap.add_argument("--nmin", type=int, default=None)
    ap.add_argument("--nmax", type=int, default=None)
    ap.add_argument("--boot", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rows = read_values(args.values)
    rows = [(n, a, s) for (n, a, s) in rows if s in EXACT_STATUSES]
    if args.nmin is not None:
        rows = [r for r in rows if r[0] >= args.nmin]
    if args.nmax is not None:
        rows = [r for r in rows if r[0] <= args.nmax]
    if not rows:
        sys.exit("No exact values in range; nothing to test.")
    rows.sort()

    nmax = max(n for n, _, _ in rows)
    sps = semiprimes_up_to_index(nmax)

    print(f"{'n':>4} {'a(n)':>6} {'y':>7} {'won':>6} "
          f"{'Z_mid':>7} {'Z_full':>7} {'H_I':>6} {'H_II':>6} {'H_III':>6}")
    z_mid, z_full = [], []
    below_y_obs = 0
    z_below, z_above = [], []
    for n, a, status in rows:
        im = IndexModel(n, sps)
        t = im.transform(a)
        won = "I" if a <= im.y else "II/III"
        if a <= im.y:
            below_y_obs += 1
            z_below.append(t["Z_mid"])
        else:
            z_above.append(t["Z_mid"])
        z_mid.append(t["Z_mid"])
        z_full.append(t["Z_full"])
        print(f"{n:>4} {a:>6} {im.y:>7.1f} {won:>6} "
              f"{t['Z_mid']:>7.3f} {t['Z_full']:>7.3f} "
              f"{t['H_I']:>6.3f} {t['H_II']:>6.3f} {t['H_III']:>6.3f}")

    z = np.array(z_mid)
    N = len(z)
    print(f"\nSample size: {N} exact terms")
    print(f"mean(Z_mid)  = {z.mean():.3f}   (Exp(1) predicts 1.000; "
          f"SE at N={N} is {1/math.sqrt(N):.3f})")
    print(f"implied global hazard deflation 1/mean(Z) = {1/z.mean():.3f}")
    print(f"median(Z_mid) = {np.median(z):.3f}   (Exp(1) predicts "
          f"{math.log(2):.3f})")

    ks = ks_exp1(z)
    print(f"\nKS vs Exp(1) [fully specified null]: "
          f"D = {ks.statistic:.3f}, p = {ks.pvalue:.4f}")

    rng = np.random.default_rng(args.seed)
    theta, d_obs, p_boot = ks_exp_fitted_bootstrap(z, args.boot, rng)
    print(f"KS vs Exp(theta), theta = mean(Z) = {theta:.3f} "
          f"[scale estimated, bootstrap p]: D = {d_obs:.3f}, "
          f"p = {p_boot:.4f}  ({args.boot} replicates)")
    print("  -> If Exp(1) is rejected but Exp(theta) is not, the model's "
          "shape is right\n     and the misfit is a uniform hazard "
          "inflation by 1/theta.")

    # Localization: below-y vs above-y races
    exp_below = sum(1.0 - math.exp(-self_hy(IndexModel(n, sps)))
                    for n, _, _ in rows)
    print(f"\nRaces won below y: observed {below_y_obs}/{N}, "
          f"model expects {exp_below:.1f}")
    if z_below and z_above:
        print(f"mean Z_mid | won below y : {np.mean(z_below):.3f}  "
              f"(N={len(z_below)})")
        print(f"mean Z_mid | won above y : {np.mean(z_above):.3f}  "
              f"(N={len(z_above)})")
        print("  -> Inflation concentrated in the below-y group indicts "
              "channel I;\n     in the above-y group, channel II's "
              "log(u-1) factor.")


def self_hy(im):
    """Cumulative model hazard at t = y for index model im (channel I only,
    since no channel II/III offset lies at or below y)."""
    H = 0.0
    m = 2
    while m <= im.y:
        if isprime(m):
            H += im.h_base
        m += 1
    return H


if __name__ == "__main__":
    main()
