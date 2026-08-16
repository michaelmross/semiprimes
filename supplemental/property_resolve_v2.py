#!/usr/bin/env python3
"""
property_resolve.py -- breadth-first resolver for d2(n)-prime challenges.

Purpose
-------
This is the second-stage companion to property_scan.py.

For a semiprimorial

    Q_n = product of the first n semiprimes,

let

    d2(n) = min { m > 1 : Q_n + m is semiprime }.

property_scan.py marks an index n as CHALLENGE when it has already found a
resolved semiprime at a COMPOSITE offset c, but one or more smaller unresolved
offsets p < c remain.  Those unresolved offsets are prime.

For the Boolean conjecture "d2(n) is prime", we do NOT need the exact d2(n):

    * If ANY unresolved prime offset p gives a semiprime Q_n+p, n is PASS.
      Stop work on that n immediately.
    * Only if EVERY unresolved prime offset is proved non-semiprime (DEAD)
      does the known composite witness c become a counterexample.

This program therefore runs factorization attempts breadth-first across all
open challenge indices, stopping an index on its first semiprime hit.

Input
-----
The challenge file written by property_scan.py, e.g.

    d2_challenges_251_999.txt

with rows

    n,m,offset_kind,witness_m,witness_kind

The expected challenge rows have offset_kind=prime and
witness_kind=composite-m.

Dependencies
------------
  * Python 3
  * gmpy2
  * GMP-ECM executable named "ecm" on PATH

Typical use
-----------
Cheap first pass:

    python3 property_resolve.py d2_challenges_251_999.txt \
        --seconds 3600 --tiers T0

Escalate only the T0 NOFACTOR survivors with a micro-tier:

    python3 property_resolve.py d2_challenges_251_999.txt \
        --seconds 28800 --tiers T1a --escalate-from T0

Continue later (journal makes runs resumable):

    python3 property_resolve.py d2_challenges_251_999.txt \
        --seconds 7200 --tiers T2 --escalate-from T1a

Deeper passes only on survivors:

    python3 property_resolve.py d2_challenges_251_999.txt \
        --seconds 21600 --tiers T3

    python3 property_resolve.py d2_challenges_251_999.txt \
        --seconds 43200 --tiers T3b

Tier definitions
----------------
T0   P-1 B1=1e6, then 6 ECM curves at B1=1e4
T1a  4 ECM curves at B1=5e4 (micro-tier for T0 NOFACTOR survivors)
T1   24 ECM curves at B1=5e4
T2   Existing-repo tier: P-1 B1=1e7, then 120 ECM curves at B1=5e4
T3   Existing-repo tier: 100 ECM curves at B1=2.5e5
T3b  Existing-repo tier: 250 ECM curves at B1=1e6

Notes on confidence
-------------------
A factor hit is classified as SEMIPRIME only when both factor and cofactor pass
gmpy2.is_prime(..., 30), matching the PRP-level logic of the existing Fortunate
scripts.  Thus PASS-PRP and COUNTEREXAMPLE-PRP are computational statements,
not formal primality certificates.

A COUNTEREXAMPLE-PRP status also relies on the challenge file's already-resolved
composite-offset semiprime witness.  This program does not re-prove that witness.

Output
------
By default, next to the challenge file:

  d2_property_resolve_journal.csv   append-only attempt journal
  d2_property_resolve_status.csv    current per-index state (rewritten)
  d2_property_resolve_open.txt      still-open n,m targets

The journal is the resume state.  Keep it.

Filtered escalation is safe: --escalate-from changes only which targets are
attempted at the higher tier. It does not remove untried targets from the
original challenge group. COUNTEREXAMPLE-PRP is still emitted only when every
original saving prime offset for that n has independently been marked DEAD.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

try:
    import gmpy2
    from gmpy2 import mpz
except ImportError as e:
    raise SystemExit(
        "property_resolve.py requires gmpy2. "
        "Install it in the same environment used for the Fortunate scripts."
    ) from e


PRP_ROUNDS = 30

# A tier is a sequence of GMP-ECM invocations.
# argv pieces are inserted after: ecm -q -one
# timeout is per invocation, in seconds.
TIER_STEPS = {
    "T0": [
        (["-pm1", "1000000"], 45),
        (["-c", "6", "10000"], 120),
    ],
    "T1a": [
        # Micro-tier for targets that already survived T0.
        # Four B1=50k curves keeps the breadth-first sweep cheap.
        (["-c", "4", "50000"], 120),
    ],
    "T1": [
        (["-c", "24", "50000"], 300),
    ],
    # Preserve the existing repository's T2/T3/T3b definitions.
    "T2": [
        (["-pm1", "10000000"], 90),
        (["-c", "120", "50000"], 300),
    ],
    "T3": [
        (["-c", "100", "250000"], 600),
    ],
    "T3b": [
        (["-c", "250", "1000000"], 3000),
    ],
}


@dataclass(frozen=True)
class ChallengeGroup:
    n: int
    witness_m: int
    witness_kind: str
    candidates: Tuple[int, ...]


@dataclass
class AttemptResult:
    outcome: str                 # SEMIPRIME, DEAD, NOFACTOR, INCOMPLETE
    factor: str = ""
    factor_digits: str = ""
    cofactor_prp: str = ""
    seconds: float = 0.0
    note: str = ""


# ---------------------------------------------------------------------------
# Small semiprime machinery: enough to reconstruct Q_n up to max challenge n.

def primes_upto(B: int) -> List[int]:
    if B < 2:
        return []
    sieve = bytearray(b"\x01") * (B + 1)
    sieve[0:2] = b"\x00\x00"
    r = int(B ** 0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            sieve[p * p : B + 1 : p] = b"\x00" * (((B - p * p) // p) + 1)
    return [i for i in range(2, B + 1) if sieve[i]]


def first_semiprimes(nmax: int) -> List[int]:
    """Return the first nmax semiprimes, counting Omega with multiplicity."""
    if nmax <= 0:
        return []
    B = max(64, 4 * nmax)
    while True:
        omega = bytearray(B + 1)
        for p in primes_upto(B):
            pk = p
            while pk <= B:
                for x in range(pk, B + 1, pk):
                    omega[x] += 1
                if pk > B // p:
                    break
                pk *= p
        sps = [x for x in range(4, B + 1) if omega[x] == 2]
        if len(sps) >= nmax:
            return sps[:nmax]
        B *= 2


def build_Q_cache(needed_n: Iterable[int]) -> Dict[int, mpz]:
    needed = sorted(set(needed_n))
    if not needed:
        return {}
    nmax = needed[-1]
    sps = first_semiprimes(nmax)
    need = set(needed)
    out: Dict[int, mpz] = {}
    Q = mpz(1)
    for i, s in enumerate(sps, start=1):
        Q *= s
        if i in need:
            out[i] = mpz(Q)
    return out


# ---------------------------------------------------------------------------
# Input parsing

def parse_challenges(path: Path) -> Dict[int, ChallengeGroup]:
    raw: Dict[int, Dict[str, object]] = {}

    with path.open("r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [x.strip() for x in line.split(",")]
            if len(parts) != 5:
                raise SystemExit(
                    f"{path}:{lineno}: expected 5 comma-separated fields, got {len(parts)}"
                )
            ns, ms, offset_kind, ws, witness_kind = parts
            n, m, witness_m = int(ns), int(ms), int(ws)

            if offset_kind != "prime":
                raise SystemExit(
                    f"{path}:{lineno}: expected offset_kind=prime, got {offset_kind!r}"
                )
            if witness_kind != "composite-m":
                raise SystemExit(
                    f"{path}:{lineno}: expected witness_kind=composite-m, got {witness_kind!r}"
                )
            if m >= witness_m:
                raise SystemExit(
                    f"{path}:{lineno}: challenge offset m={m} is not below witness {witness_m}"
                )

            if n not in raw:
                raw[n] = {
                    "witness_m": witness_m,
                    "witness_kind": witness_kind,
                    "candidates": [],
                }
            else:
                if int(raw[n]["witness_m"]) != witness_m:
                    raise SystemExit(
                        f"{path}:{lineno}: inconsistent witness_m for n={n}"
                    )
                if str(raw[n]["witness_kind"]) != witness_kind:
                    raise SystemExit(
                        f"{path}:{lineno}: inconsistent witness_kind for n={n}"
                    )
            raw[n]["candidates"].append(m)  # type: ignore[index]

    groups: Dict[int, ChallengeGroup] = {}
    for n, rec in raw.items():
        candidates = tuple(sorted(set(int(x) for x in rec["candidates"])))  # type: ignore[index]
        if not candidates:
            continue
        groups[n] = ChallengeGroup(
            n=n,
            witness_m=int(rec["witness_m"]),
            witness_kind=str(rec["witness_kind"]),
            candidates=candidates,
        )
    return groups


# ---------------------------------------------------------------------------
# Journal / state

JOURNAL_FIELDS = [
    "timestamp_utc",
    "n",
    "m",
    "tier",
    "outcome",
    "factor",
    "factor_digits",
    "cofactor_prp",
    "seconds",
    "digits_N",
    "note",
]


def load_journal(path: Path) -> List[dict]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def append_journal(path: Path, row: dict) -> None:
    exists = path.exists() and path.stat().st_size > 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=JOURNAL_FIELDS)
        if not exists:
            w.writeheader()
        w.writerow(row)
        fh.flush()
        os.fsync(fh.fileno())


def reconstruct_state(
    groups: Dict[int, ChallengeGroup],
    journal_rows: Sequence[dict],
) -> Tuple[
    Dict[int, Set[int]],
    Dict[int, Optional[int]],
    Set[Tuple[int, int, str]],
    Dict[Tuple[int, int, str], str],
]:
    """
    Return:
      dead[n]              offsets proved non-semiprime
      semiprime_hit[n]     first recorded saving prime offset, or None
      tier_done            (n,m,tier) for completed attempts
      tier_outcome         latest conclusive outcome for (n,m,tier)

    tier_outcome supports --escalate-from, so a higher tier can be restricted
    to targets whose lower-tier result was specifically NOFACTOR.
    """
    dead: Dict[int, Set[int]] = {n: set() for n in groups}
    semiprime_hit: Dict[int, Optional[int]] = {n: None for n in groups}
    tier_done: Set[Tuple[int, int, str]] = set()
    tier_outcome: Dict[Tuple[int, int, str], str] = {}

    valid_candidates = {n: set(g.candidates) for n, g in groups.items()}

    for r in journal_rows:
        try:
            n, m = int(r["n"]), int(r["m"])
            tier = r["tier"]
            outcome = r["outcome"]
        except (KeyError, ValueError):
            continue

        if n not in groups or m not in valid_candidates[n]:
            continue

        if outcome in {"SEMIPRIME", "DEAD", "NOFACTOR"}:
            tier_done.add((n, m, tier))
            tier_outcome[(n, m, tier)] = outcome
        # INCOMPLETE deliberately does not mark the tier complete and does not
        # overwrite an earlier conclusive outcome.

        if outcome == "DEAD":
            dead[n].add(m)
        elif outcome == "SEMIPRIME":
            old = semiprime_hit[n]
            semiprime_hit[n] = m if old is None else min(old, m)

    return dead, semiprime_hit, tier_done, tier_outcome


def group_status(
    g: ChallengeGroup,
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
) -> str:
    if semiprime_hit[g.n] is not None:
        return "PASS-PRP"
    if all(m in dead[g.n] for m in g.candidates):
        return "COUNTEREXAMPLE-PRP"
    return "OPEN"


def summarize(
    groups: Dict[int, ChallengeGroup],
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    for g in groups.values():
        counts[group_status(g, dead, semiprime_hit)] += 1
    return dict(counts)


def write_status(
    path: Path,
    groups: Dict[int, ChallengeGroup],
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
) -> None:
    fields = [
        "n", "status", "witness_m", "witness_kind",
        "total_targets", "dead_targets", "open_targets", "semiprime_m",
    ]
    rows = []
    for n in sorted(groups):
        g = groups[n]
        status = group_status(g, dead, semiprime_hit)
        open_targets = [
            m for m in g.candidates
            if m not in dead[n] and semiprime_hit[n] is None
        ]
        rows.append({
            "n": n,
            "status": status,
            "witness_m": g.witness_m,
            "witness_kind": g.witness_kind,
            "total_targets": len(g.candidates),
            "dead_targets": len(dead[n] & set(g.candidates)),
            "open_targets": len(open_targets),
            "semiprime_m": "" if semiprime_hit[n] is None else semiprime_hit[n],
        })

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


def write_open_targets(
    path: Path,
    groups: Dict[int, ChallengeGroup],
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write("# Still-open property-only d2 targets.\n")
        fh.write("# n,m,witness_m\n")
        for n in sorted(groups):
            g = groups[n]
            if group_status(g, dead, semiprime_hit) != "OPEN":
                continue
            for m in g.candidates:
                if m not in dead[n]:
                    fh.write(f"{n},{m},{g.witness_m}\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# GMP-ECM

def extract_factor(text: str, N: mpz) -> Optional[mpz]:
    best: Optional[mpz] = None
    for tok in text.replace(":", " ").replace("=", " ").split():
        if tok.isdigit():
            t = mpz(tok)
            if 1 < t < N and N % t == 0:
                if best is None or t < best:
                    best = t
    return best


def run_ecm_step(
    N: mpz,
    args: Sequence[str],
    timeout: float,
) -> Tuple[Optional[mpz], bool, str]:
    """
    Return (factor, completed, note).
    completed=False means the subprocess timed out.
    """
    try:
        r = subprocess.run(
            ["ecm", "-q", "-one"] + list(args),
            input=str(N),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as e:
        text = ""
        if e.stdout:
            text += e.stdout if isinstance(e.stdout, str) else e.stdout.decode(errors="ignore")
        if e.stderr:
            text += " " + (e.stderr if isinstance(e.stderr, str) else e.stderr.decode(errors="ignore"))
        d = extract_factor(text, N)
        if d is not None:
            return d, False, "factor recovered from timed-out ECM output"
        return None, False, f"timeout after {timeout:.0f}s"

    text = (r.stdout or "") + " " + (r.stderr or "")
    d = extract_factor(text, N)
    if d is not None:
        return d, True, ""
    return None, True, f"ecm exit={r.returncode}" if r.returncode != 0 else ""


def run_tier(
    N: mpz,
    tier: str,
    global_deadline: float,
    timeout_scale: float,
) -> AttemptResult:
    t0 = time.time()

    for args, base_timeout in TIER_STEPS[tier]:
        left = global_deadline - time.time()
        if left <= 1.0:
            return AttemptResult(
                outcome="INCOMPLETE",
                seconds=time.time() - t0,
                note="global deadline reached before tier completed",
            )

        timeout = min(base_timeout * timeout_scale, max(1.0, left))
        d, completed, note = run_ecm_step(N, args, timeout)

        if d is not None:
            c = N // d
            dp = bool(gmpy2.is_prime(d, PRP_ROUNDS))
            cp = bool(gmpy2.is_prime(c, PRP_ROUNDS))
            verdict = "SEMIPRIME" if (dp and cp) else "DEAD"
            return AttemptResult(
                outcome=verdict,
                factor=str(d),
                factor_digits=str(d.num_digits()),
                cofactor_prp=str(bool(dp and cp)),
                seconds=time.time() - t0,
                note=note,
            )

        if not completed:
            return AttemptResult(
                outcome="INCOMPLETE",
                seconds=time.time() - t0,
                note=note,
            )

    return AttemptResult(
        outcome="NOFACTOR",
        seconds=time.time() - t0,
    )


# ---------------------------------------------------------------------------
# Breadth-first scheduling

def make_round_robin_schedule(
    groups: Dict[int, ChallengeGroup],
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
    tier_done: Set[Tuple[int, int, str]],
    tier_outcome: Dict[Tuple[int, int, str], str],
    tier: str,
    escalate_from: Optional[str] = None,
) -> List[Tuple[int, int]]:
    """
    Breadth-first across n:
      first pending target from every open n,
      then second pending target from every open n, etc.

    Within an n, offsets stay in ascending order for reproducibility/auditability.
    """
    pending: Dict[int, List[int]] = {}
    for n, g in groups.items():
        if group_status(g, dead, semiprime_hit) != "OPEN":
            continue
        xs = []
        for m in g.candidates:
            if m in dead[n] or (n, m, tier) in tier_done:
                continue
            if escalate_from is not None:
                # Escalate only targets that completed the named lower tier
                # with NOFACTOR. Untried/INCOMPLETE lower-tier targets are
                # excluded, but nothing is inferred about them.
                if tier_outcome.get((n, m, escalate_from)) != "NOFACTOR":
                    continue
            xs.append(m)
        if xs:
            pending[n] = xs

    if not pending:
        return []

    # Put smaller open burdens first inside each breadth layer.
    ns = sorted(pending, key=lambda n: (len(pending[n]), n))
    max_len = max(len(xs) for xs in pending.values())
    schedule: List[Tuple[int, int]] = []
    for rank in range(max_len):
        for n in ns:
            xs = pending[n]
            if rank < len(xs):
                schedule.append((n, xs[rank]))
    return schedule


def print_summary(
    groups: Dict[int, ChallengeGroup],
    dead: Dict[int, Set[int]],
    semiprime_hit: Dict[int, Optional[int]],
) -> None:
    counts = summarize(groups, dead, semiprime_hit)
    open_targets = 0
    dead_targets = 0
    for n, g in groups.items():
        dead_targets += len(dead[n] & set(g.candidates))
        if group_status(g, dead, semiprime_hit) == "OPEN":
            open_targets += sum(1 for m in g.candidates if m not in dead[n])

    print("\nSummary", flush=True)
    for k in ("PASS-PRP", "OPEN", "COUNTEREXAMPLE-PRP"):
        print(f"  {k:22s} {counts.get(k, 0)}", flush=True)
    print(f"  DEAD targets           {dead_targets}", flush=True)
    print(f"  OPEN targets           {open_targets}", flush=True)


# ---------------------------------------------------------------------------
# Driver

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Breadth-first property-only resolver for d2(n)-prime challenges."
    )
    ap.add_argument("challenge_file", type=Path)
    ap.add_argument(
        "--seconds", type=float, default=3600.0,
        help="wall-clock budget for this invocation (default: 3600)"
    )
    ap.add_argument(
        "--tiers", default="T0,T1",
        help="comma-separated tiers to run, in order (default: T0,T1)"
    )
    ap.add_argument(
        "--escalate-from", default=None,
        help=(
            "only run selected tier(s) on targets whose named lower tier "
            "finished with NOFACTOR; e.g. --tiers T1a --escalate-from T0"
        )
    )
    ap.add_argument(
        "--timeout-scale", type=float, default=1.0,
        help="multiply per-ECM-call timeouts by this factor (default: 1.0)"
    )
    ap.add_argument(
        "--journal", type=Path, default=None,
        help="append-only journal CSV (default: beside challenge file)"
    )
    ap.add_argument(
        "--status-out", type=Path, default=None,
        help="current per-n status CSV (default: beside challenge file)"
    )
    ap.add_argument(
        "--open-out", type=Path, default=None,
        help="current still-open target list (default: beside challenge file)"
    )
    ap.add_argument(
        "--stop-on-counterexample", action="store_true",
        help="stop immediately if all saving prime offsets for an n are killed"
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="parse/reconstruct/report state but do not call ECM"
    )
    args = ap.parse_args()

    challenge_file = args.challenge_file.resolve()
    if not challenge_file.exists():
        raise SystemExit(f"challenge file not found: {challenge_file}")

    base = challenge_file.parent
    journal = (args.journal or (base / "d2_property_resolve_journal.csv")).resolve()
    status_out = (args.status_out or (base / "d2_property_resolve_status.csv")).resolve()
    open_out = (args.open_out or (base / "d2_property_resolve_open.txt")).resolve()

    tiers = [x.strip() for x in args.tiers.split(",") if x.strip()]
    if not tiers:
        raise SystemExit("no tiers selected")
    bad = [t for t in tiers if t not in TIER_STEPS]
    if bad:
        raise SystemExit(
            f"unknown tier(s): {', '.join(bad)}; choose from {', '.join(TIER_STEPS)}"
        )
    if args.escalate_from is not None and args.escalate_from not in TIER_STEPS:
        raise SystemExit(
            f"unknown --escalate-from tier {args.escalate_from!r}; "
            f"choose from {', '.join(TIER_STEPS)}"
        )
    if args.seconds <= 0:
        raise SystemExit("--seconds must be > 0")
    if args.timeout_scale <= 0:
        raise SystemExit("--timeout-scale must be > 0")

    groups = parse_challenges(challenge_file)
    if not groups:
        raise SystemExit("no challenge rows found")

    journal_rows = load_journal(journal)
    dead, semiprime_hit, tier_done, tier_outcome = reconstruct_state(groups, journal_rows)

    print(
        f"Loaded {len(groups)} challenge indices and "
        f"{sum(len(g.candidates) for g in groups.values())} prime-offset targets.",
        flush=True,
    )
    print(f"Journal: {journal}", flush=True)
    print_summary(groups, dead, semiprime_hit)
    write_status(status_out, groups, dead, semiprime_hit)
    write_open_targets(open_out, groups, dead, semiprime_hit)

    if args.dry_run:
        print(f"\nDry run only. Status: {status_out}\nOpen:   {open_out}", flush=True)
        return

    if shutil.which("ecm") is None:
        raise SystemExit(
            "GMP-ECM executable 'ecm' was not found on PATH. "
            "Install GMP-ECM or activate the environment used by your existing scripts."
        )

    # Build only the Q_n values that can still matter.
    open_ns = [
        n for n, g in groups.items()
        if group_status(g, dead, semiprime_hit) == "OPEN"
    ]
    if not open_ns:
        print("\nNo open challenge indices remain.", flush=True)
        return

    print(f"\nBuilding semiprimorials through n={max(open_ns)} ...", flush=True)
    Qcache = build_Q_cache(open_ns)

    deadline = time.time() + args.seconds
    attempted_this_run: Set[Tuple[int, int, str]] = set()
    total_attempts = 0

    for tier in tiers:
        if time.time() >= deadline:
            break

        schedule = make_round_robin_schedule(
            groups,
            dead,
            semiprime_hit,
            tier_done,
            tier_outcome,
            tier,
            args.escalate_from,
        )
        filter_note = (
            f"; only {args.escalate_from}=NOFACTOR survivors"
            if args.escalate_from is not None
            else ""
        )
        print(
            f"\n[{tier}] {len(schedule)} currently eligible target-tier attempts "
            f"(breadth-first{filter_note}).",
            flush=True,
        )

        for n, m in schedule:
            if time.time() >= deadline:
                print("[global deadline]", flush=True)
                break

            g = groups[n]
            if group_status(g, dead, semiprime_hit) != "OPEN":
                continue
            if m in dead[n]:
                continue
            if (n, m, tier) in tier_done:
                continue
            if (n, m, tier) in attempted_this_run:
                continue

            attempted_this_run.add((n, m, tier))
            N = Qcache[n] + m
            digits_N = int(N.num_digits())

            print(
                f"{tier} n={n} m={m}  ({digits_N} digits; "
                f"witness={g.witness_m}) ...",
                end=" ",
                flush=True,
            )

            result = run_tier(N, tier, deadline, args.timeout_scale)
            total_attempts += 1

            row = {
                "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "n": n,
                "m": m,
                "tier": tier,
                "outcome": result.outcome,
                "factor": result.factor,
                "factor_digits": result.factor_digits,
                "cofactor_prp": result.cofactor_prp,
                "seconds": f"{result.seconds:.3f}",
                "digits_N": digits_N,
                "note": result.note,
            }
            append_journal(journal, row)

            if result.outcome in {"SEMIPRIME", "DEAD", "NOFACTOR"}:
                tier_done.add((n, m, tier))
                tier_outcome[(n, m, tier)] = result.outcome

            if result.outcome == "SEMIPRIME":
                semiprime_hit[n] = m
                print(
                    f"SEMIPRIME -> PASS-PRP "
                    f"(factor {result.factor_digits} digits, {result.seconds:.1f}s)",
                    flush=True,
                )
                # This is the whole point: n is closed immediately.
            elif result.outcome == "DEAD":
                dead[n].add(m)
                print(
                    f"DEAD (factor {result.factor_digits} digits, {result.seconds:.1f}s)",
                    flush=True,
                )
                if group_status(g, dead, semiprime_hit) == "COUNTEREXAMPLE-PRP":
                    print(
                        f"\n*** COUNTEREXAMPLE-PRP at n={n}: "
                        f"all {len(g.candidates)} saving prime offsets are DEAD; "
                        f"composite witness m={g.witness_m} remains first. ***\n",
                        flush=True,
                    )
                    write_status(status_out, groups, dead, semiprime_hit)
                    write_open_targets(open_out, groups, dead, semiprime_hit)
                    if args.stop_on_counterexample:
                        print_summary(groups, dead, semiprime_hit)
                        return
            elif result.outcome == "NOFACTOR":
                print(f"no factor ({result.seconds:.1f}s)", flush=True)
            else:
                print(f"INCOMPLETE ({result.seconds:.1f}s): {result.note}", flush=True)

            # Persist human-readable current state after every conclusive target.
            if result.outcome in {"SEMIPRIME", "DEAD"}:
                write_status(status_out, groups, dead, semiprime_hit)
                write_open_targets(open_out, groups, dead, semiprime_hit)

        if time.time() >= deadline:
            break

    write_status(status_out, groups, dead, semiprime_hit)
    write_open_targets(open_out, groups, dead, semiprime_hit)

    print(f"\nAttempts this run: {total_attempts}", flush=True)
    print_summary(groups, dead, semiprime_hit)
    print(f"\nStatus:  {status_out}", flush=True)
    print(f"Open:    {open_out}", flush=True)
    print(f"Journal: {journal}", flush=True)


if __name__ == "__main__":
    main()
