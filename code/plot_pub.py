import numpy as np, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "figure.dpi": 200})

import os
_cands = ["semi_results.npz", "data/semi_results.npz",
          os.path.join(os.path.dirname(__file__), "..", "data", "semi_results.npz")]
_npz = next((p for p in _cands if os.path.exists(p)), _cands[0])
d = np.load(_npz)
reached = int(d["reached"]); nmax = int(d["nmax_valid"])
Gcum = d["Gcum"]; gap_hist = d["gap_hist"]
E2 = d["E2arr"]; dec_semi = d["dec_semi"]

# ---------- Fig 1: worst-case gap ----------
ks = np.arange(4, 12); xs = 10.0**ks; G = Gcum[ks].astype(float)
lx = np.log(xs)
c = np.mean(G / lx**2)                       # fitted constant for c (log x)^2
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
a1.loglog(xs, G, "o-", lw=2, ms=7, color="#1f4e79", label=r"$G_{E_2}(x)$ observed")
a1.loglog(xs, np.sqrt(xs), "--", color="#d9730d", label=r"$\sqrt{x}$ (threshold)")
a1.loglog(xs, lx**2, ":", color="#2e8b57", label=r"$(\log x)^2$")
a1.loglog(xs, c*lx**2, "-.", color="#999999", lw=1.2, label=fr"$%.2f(\log x)^2$ fit" % c)
a1.set_xlabel("x"); a1.set_ylabel("gap"); a1.grid(True, which="both", alpha=.25)
a1.set_title("Worst-case exact-semiprime gap"); a1.legend(fontsize=8.5)
a2.semilogx(xs, G/np.sqrt(xs), "s-", color="crimson", lw=2, ms=6)
a2.set_xlabel("x"); a2.set_ylabel(r"$G_{E_2}(x)/\sqrt{x}$")
a2.set_title(r"Gap relative to $\sqrt{x}$ (collapsing)"); a2.grid(True, which="both", alpha=.25)
fig.tight_layout(); fig.savefig("fig1_gap_vs_sqrt.png"); plt.close(fig)

# ---------- Fig 2: detrended square intervals ----------
n = np.arange(2, nmax + 1); e2 = E2[2:nmax + 1].astype(float)
x = n*(n + 1.0); Ehat = (2*n + 1)*np.log(np.log(x))/np.log(x); r = e2/Ehat
odd = n % 2 == 1
edges = np.unique(np.round(np.geomspace(2, nmax + 1, 48)).astype(int))
cen, ro, re_, rmin = [], [], [], []
for i in range(len(edges)-1):
    m = (n >= edges[i]) & (n < edges[i+1])
    if m.sum() < 5: continue
    cen.append(math.sqrt(edges[i]*edges[i+1]))
    ro.append(r[m & odd].mean()); re_.append(r[m & ~odd].mean()); rmin.append(r[m].min())
fig, ax = plt.subplots(figsize=(8, 4.6))
ax.semilogx(cen, ro, "o-", ms=4, color="#1f4e79", label="odd n (mean)")
ax.semilogx(cen, re_, "s--", ms=4, color="#d9730d", label="even n (mean)")
ax.semilogx(cen, rmin, ".:", color="gray", label="min in bin")
ax.axhline(1.0, color="k", lw=.7)
ax.set_xlabel("n"); ax.set_ylabel(r"$E_2(n)\,/\,\widehat{E_2}(n)$")
ax.set_title(r"Detrended count in $(n^2,(n+1)^2]$ to $n=316{,}226$ — odd vs even")
ax.grid(True, which="both", alpha=.25); ax.legend(fontsize=9)
fig.tight_layout(); fig.savefig("fig2_square_detrend.png"); plt.close(fig)

# ---------- Fig 3: gap self-similarity (honest title) ----------
bins = np.arange(gap_hist.shape[1])
fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4.2))
cvs, kk = [], []
for k in range(2, len(gap_hist)):
    h = gap_hist[k].astype(float); tot = h.sum()
    if tot == 0 or 10**(k+1) > reached: continue
    mean = (bins*h).sum()/tot; var = ((bins-mean)**2*h).sum()/tot
    cvs.append(math.sqrt(var)/mean); kk.append(k)
    sel = h > 0
    b1.plot(bins[sel]/mean, h[sel]/tot*mean, drawstyle="steps-mid", alpha=.7,
            label=fr"$10^{{{k}}}$–$10^{{{k+1}}}$")
b1.set_xlabel("gap / mean gap"); b1.set_ylabel("density (scaled)")
b1.set_title("Normalized gap distribution"); b1.set_xlim(0, 8)
b1.grid(alpha=.25); b1.legend(fontsize=7, ncol=1)
b2.plot(kk, cvs, "o-", lw=2, color="#1f4e79")
b2.set_xlabel(r"decade $k$  ($x\in[10^k,10^{k+1})$)"); b2.set_ylabel("coefficient of variation")
b2.set_title(r"CV drifts $0.72\!\to\!0.83$ (shape stable, tail broadens)")
b2.set_ylim(0, 1); b2.grid(alpha=.25)
fig.tight_layout(); fig.savefig("fig3_gap_selfsim.png"); plt.close(fig)

# print exact numbers for the writeup
print("fit c (log x)^2:", round(c, 3))
for k in range(4, 11):
    x = 10**k; lan = x*math.log(math.log(x))/math.log(x)
    print(f"1e{k}: ratio {int(dec_semi[k])/lan:.4f}")
tot = 13959990342  # exact count < 1e11 from run log
lan11 = 1e11*math.log(math.log(1e11))/math.log(1e11)
print(f"1e11 total {tot}  ratio {tot/lan11:.4f}")
print("done")
