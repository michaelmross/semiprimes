# Semiprimorial Triad

This directory contains the computational material supporting the paper

> M. M. Ross, *A Semiprimorial Triad: Prime First Hits, Semiprime First Hits, and Semiprime Recurrence* (2026).
> [doi.org/10.5281/zenodo.21969573](https://doi.org/10.5281/zenodo.21969573)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21970222.svg)](https://doi.org/10.5281/zenodo.21970222)

Let

$$
Q_n=\prod_{j=1}^n s_j,
$$

where $$s_j$$ is the $$j$$-th semiprime, and define

$$
d_1(n)=\min\{m>1:\Omega(Q_n+m)=1\},
\qquad
d_2(n)=\min\{m>1:\Omega(Q_n+m)=2\}.
$$

The paper studies three related phenomena:

1. $$d_1(n)$$ appears always to be prime;
2. $$d_2(n)$$ appears always to be prime;
3. $$Q_n+1$$ appears to be semiprime infinitely often.

The files here reproduce the finite results reported in the paper. In particular:

- all $$100$$ values $$d_1(n)$$, $$1\le n\le100$$, are determined exactly and proved prime;
- the **prime-property** of $$d_2(n)$$ is rigorously verified for every $$1\le n\le100$$, without requiring every exact minimum to be known;
- thirteen values of $$Q_n+1$$ with $$n\le100$$ are rigorously certified semiprime, while six cases remain deliberately unclassified.

## Repository layout

```text
semiprimorial/
├── d1/             exact d1 verification through n = 100
├── d2/             d2 prime-property verification through n = 100
├── Qn+1/         fixed-shift Q_n + 1 certificates and open cases
├── supplemental/   campaign and factor-discovery Python programs
└── README.md
```

The first three folders contain the material needed for the finite claims in the paper. The `supplemental/` folder is research provenance: it records useful search and factor-discovery code but is not required to verify the stated finite theorems.

---

## 1. `d1/`: exact prime first hits through $$n=100$$

The paper defines

$$
d_1(n)=\min\{m>1:Q_n+m\text{ is prime}\}.
$$

Every $$d_1(n)$$ through $$n=100$$ is determined exactly.

Files:

- `finite100_d1_witness_manifest.csv` — the $$100$$ exact values, together with $$s_n$$, $$q_0$$, the $$q_0^2$$ protection boundary, and related audit data;
- `finite100_d1_certify.gp` — independently reconstructs the semiprimorials, checks every smaller offset, proves $$Q_n+d_1(n)$$ prime, proves $$d_1(n)$$ prime, and generates ECPP certificates where needed;
- `finite100_d1_verify.gp` — independently reconstructs the data and validates the saved certificates;
- `finite100_d1_certify.log` — certification transcript;
- `finite100_d1_verify.log` — independent verification transcript.

To rerun:

```bash
cd d1

gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  finite100_d1_certify.gp 2>&1 | tee finite100_d1_certify.log

gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  finite100_d1_verify.gp 2>&1 | tee finite100_d1_verify.log
```

The expected final verifier message is:

```text
ALL 100 EXACT d1 VALUES AND SAVED CERTIFICATES VALID
```

---

## 2. `d2/`: prime-property verification through $$n=100$$

Here

$$
d_2(n)=\min\{m>1:\Omega(Q_n+m)=2\}.
$$

This computation is intentionally different from the $$d_1$$ computation. Large semiprime candidates can be difficult to factor completely, so the paper verifies the **property** that the true minimum has a prime offset rather than insisting on determining every exact value.

The method has two layers:

1. an independent audit exhausts the offsets below each recorded scan front-runner and checks that no composite-offset semiprime can overturn the prime-property;
2. a positive certificate supplies a rigorously proved semiprime at a prime offset for every $$n\le100$$.

Files:

- `finite100_witness_manifest.csv` — one certified prime-offset semiprime witness for every $$1\le n\le100$$;
- `finite100_certify_v2.gp` — proves the witness factorizations and primality of their factors;
- `finite100_verify_v2.gp` — independently rechecks the factorizations and saved certificates;
- `finite100_certify_v2.log` — certification transcript;
- `finite100_verify_v2.log` — verification transcript;
- `audit_scan_v2.py` — authoritative independent audit of the scan layer;
- `scan_1_65.json`;
- `scan_66_92.json`;
- `scan_93_117.json`;
- `audit_v2_1_65.log`;
- `audit_v2_66_92.log`;
- `audit_v2_93_117.log`.

The older `audit_scan.py` is superseded; `audit_scan_v2.py` is the version used for the paper.

### Rerun the independent scan audit

```bash
cd d2

python3 audit_scan_v2.py scan_1_65.json   2>&1 | tee audit_v2_1_65.log
python3 audit_scan_v2.py scan_66_92.json  2>&1 | tee audit_v2_66_92.log
python3 audit_scan_v2.py scan_93_117.json 2>&1 | tee audit_v2_93_117.log
```

Each audit should finish by reporting that the metadata, medium-prime filtering, classifications, deferred sets, alive-composite sets, front-runners, and PRP cofactor checks passed, with no resolved semiprime below a recorded scan front-runner.

### Certify and verify the positive witnesses

```bash
gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  finite100_certify_v2.gp 2>&1 | tee finite100_certify_v2.log

gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  finite100_verify_v2.gp 2>&1 | tee finite100_verify_v2.log
```

The expected final verifier message is:

```text
ALL 100 WITNESS FACTORIZATIONS AND SAVED CERTIFICATES VALID
```

### What this proves

For every $$1\le n\le100$$,

$$
d_2(n)\ \text{is prime}.
$$

This is a property verification. Some exact $$d_2(n)$$ values beyond the fully resolved initial range may still depend on unresolved factorizations, but any such unresolved competitor in the audited race lies at a prime offset and therefore cannot invalidate the theorem.

---

## 3. `Qn+1/`: the fixed shift $$Q_n+1$$

The third part of the paper asks whether

$$
\Omega(Q_n+1)=2
$$

for infinitely many $$n$$.

Through $$n=100$$, thirteen semiprime occurrences have been rigorously established:

```text
2, 3, 8, 9, 15, 16, 19, 21, 23, 27, 29, 43, 65
```

The first four are small and are factored explicitly in the paper. The nine later cases are certified here.

Files:

- `qplus1_semiprime9_certify_v2.gp` — reconstructs $$Q_n+1$$, checks the exact two-factor decompositions, and proves both factors prime;
- `qplus1_semiprime9_verify_v2.gp` — independently validates the saved certificates;
- `qplus1_semiprime9_certify_v2.log` — certification transcript;
- `qplus1_semiprime9_verify_v2.log` — verification transcript;
- `qplus1_open_after_resolve.txt` — the six cases still unclassified in the $$n\le100$$ census.

The unresolved indices are

```text
52, 72, 74, 77, 80, 86
```

No conclusion in the paper depends on how these six cases eventually factor. The finite statement is therefore the lower bound

$$
\#\{n\le100:\Omega(Q_n+1)=2\}\ge13,
$$

not an assertion that $$13$$ is the exact count.

To rerun the nine later certifications:

```bash
cd Qn+1

gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  qplus1_semiprime9_certify_v2.gp 2>&1 | tee qplus1_semiprime9_certify_v2.log

gp -fq -D parisize=64000000 -D parisizemax=2000000000 \
  qplus1_semiprime9_verify_v2.gp 2>&1 | tee qplus1_semiprime9_verify_v2.log
```

The expected final verifier message is:

```text
ALL 9 Q_n+1 SEMIPRIME FACTORIZATIONS AND SAVED CERTIFICATES VALID
```

---

## 4. `supplemental/`: campaign and factor-discovery code

This folder contains Python programs used during the broader computational investigation. They are retained as provenance and as tools for extending the search, but they are **not dependencies of the finite verification above**.

Examples include programs for:

- producing baseline semiprimorial scans;
- targeted factor discovery;
- resumable ECM campaigns;
- campaign bookkeeping;
- the $$Q_n+1$$ FactorDB/local-resolution campaign;
- model and hazard diagnostics;
- b-file target generation and value finalization.

Only the Python source programs are retained in this supplemental folder; campaign journals, provisional CSVs, deeper scan outputs, and other transient research artifacts are not needed for the paper's verification package.

GMP-ECM is useful for some of these supplemental factor-discovery programs, but it is **not required** to verify the finite claims in `d1/`, `d2/`, or `Qn+1/`.

---

## Software

The core verification uses:

- Python 3;
- `gmpy2`;
- `sympy`;
- PARI/GP with `isprime`, `primecert`, `primecertisvalid`, and `primecertexport`.

For example, the Python dependencies can be installed with:

```bash
python3 -m pip install gmpy2 sympy
```

The PARI/GP certification scripts create `.cert.gp` and `.primo` files for large primes. These are **generated outputs**: the certifier creates them and the verifier reads them back and validates them. They therefore do not need to be stored in the Git repository in order to reproduce the verification.

Likewise, the derived-state JSON files written by `audit_scan_v2.py` are reproducible outputs and need not be committed.

---

## Reproducibility scope

The repository distinguishes three kinds of computational statement:

- **exact determination** — used for $$d_1(n)$$ through $$100$$;
- **property verification** — used for $$d_2(n)$$ through $$100$$;
- **finite recurrence evidence** — used for $$Q_n+1$$, where thirteen hits are proved and six cases remain open.

The first two are finite theorems. The third is evidence for an infinitude conjecture: no finite computation can prove or disprove that infinitely many $$Q_n+1$$ are semiprime.

No verification claim here depends on a live FactorDB entry.

---

## Citation

If you use this census, please cite the repository archive on Zenodo — see the
badge in the top-level `README.md` — and the accompanying paper.

## License

Code MIT, data CC BY 4.0, as for the rest of the repository.
