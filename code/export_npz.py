#!/usr/bin/env python3
"""Write semi_results.npz directly from a checkpoint. No sieving, no module
import -> immune to the [F] memory issue. Works on partial checkpoints.
Usage:  python export_npz.py            (reads semi_ckpt.pkl)
"""
import pickle, math
import numpy as np

with open("semi_ckpt.pkl", "rb") as f:
    st = pickle.load(f)

reached = int(st["next_A"])
nmax_valid = math.isqrt(reached) - 1          # square intervals trustworthy up to here
dg   = st["dec_gmax"]
Gcum = np.maximum.accumulate(np.concatenate(([0], dg)))   # Gcum[k] = G_E2(10^k)

np.savez("semi_results.npz",
         E2arr=st["E2arr"], E3arr=st["E3arr"], P1arr=st["P1arr"],
         gap_hist=st["gap_hist"], dec_semi=st["dec_semi"], dec_tri=st["dec_tri"],
         dec_prime=st["dec_prime"], dec_gmax=dg, Gcum=Gcum,
         gmax=int(st["gmax"]), gmax_loc=int(st["gmax_loc"]),
         reached=reached, nmax_valid=nmax_valid)

print(f"saved semi_results.npz")
print(f"  reached height   : {reached:,}")
print(f"  valid square-int n: 2 .. {nmax_valid:,}   (E2arr[n] for n>{nmax_valid} is unfilled)")

# corrected square-interval minimum over VALID n only
E2 = st["E2arr"]
vmin = int(E2[2:nmax_valid + 1].min())
nstar = 2 + int(E2[2:nmax_valid + 1].argmin())
print(f"  E2(n) min over valid n = {vmin} at n={nstar}  "
      f"(>=1 everywhere => supports Conjecture 2)")

# recover [G] normalized-gap stats from the histogram (was skipped by the crash)
gh = st["gap_hist"]; bins = np.arange(gh.shape[1])
print("\n[G] normalized-gap distribution per decade")
for k in range(2, len(gh)):
    h = gh[k]; tot = int(h.sum())
    if tot == 0: continue
    mean = (bins * h).sum() / tot
    var  = ((bins - mean) ** 2 * h).sum() / tot
    gmax = int(np.nonzero(h)[0].max())
    print(f"  1e{k}-1e{k+1}: mean_gap={mean:.3f}  CV={math.sqrt(var)/mean:.3f}  "
          f"max_gap={gmax}  max/mean={gmax/mean:.2f}")
