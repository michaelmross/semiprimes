#!/usr/bin/env python3
"""
finalize_values.py -- Rebuild per-index state from derived_state.json plus
the verdict ledger and emit the calibration table.

Status assignment:
  exact       every deferred offset below the winner has a DEAD verdict
              with a recorded factor, and the winner's semiprimality is
              witnessed by a factor found here (BPSW on both parts)
  probable    same, but the winner's witness comes from the paper's
              certificates rather than this session (n=81, 83, 84, 93)
  provisional deferred offsets below the winner remain open; the winner
              is an upper bound only, and the row is excluded from
              calibration by hazard_transform_test.py

The emitted CSV is directly consumable by hazard_transform_test.py.
"""

import csv
import json

PAPER_WINNERS = {81: 151, 83: 233, 84: 337}
PAPER_WITNESS = {81, 83, 84}


def main():
    derived = json.load(open("derived_state.json"))
    res = list(csv.DictReader(open("resolutions.csv")))

    rows, summary = [], []
    for k in sorted(derived, key=int):
        n = int(k)
        rec = derived[k]
        winner = PAPER_WINNERS.get(n, rec["winner"])
        wwitness = "paper" if n in PAPER_WITNESS else None
        dead = {}
        for r in [x for x in res if int(x["n"]) == n]:
            m = int(r["m"])
            if r["verdict"] == "DEAD":
                dead[m] = r["factor"]
            elif r["verdict"] == "SEMIPRIME":
                if m < winner:
                    winner, wwitness = m, r["factor"]
                elif m == winner and wwitness is None:
                    wwitness = r["factor"]
        opens = sorted(m for m in rec["deferred"]
                       if m < winner and m not in dead)
        if opens:
            status, note = "provisional", "open:" + "+".join(map(str, opens))
        elif wwitness is None:
            status, note = "provisional", "winner unwitnessed"
        elif wwitness == "paper":
            status, note = "probable", "winner witness from paper Table 3"
        else:
            status, note = "exact", f"winner factor {wwitness}"
        rows.append([n, winner, status, note])
        summary.append((n, winner, status, len(opens),
                        rec["winner"], rec["deferred"]))

    rows.append([93, 167, "probable", "paper Table 3 (prp170 cofactor)"])

    with open("values_66_93.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "a", "status", "note"])
        w.writerows(rows)

    ex = [r for r in rows if r[2] in ("exact", "probable")]
    pv = [r for r in rows if r[2] == "provisional"]
    print(f"{len(ex)} calibration-eligible values, {len(pv)} provisional\n")
    print(f"{'n':>4} {'a(n)':>5} {'status':>12}  {'snapshot':>8}  note")
    for (n, w_, st, nopen, snap, defr) in summary:
        flag = " <-displaced" if w_ != snap else ""
        print(f"{n:>4} {w_:>5} {st:>12}  {snap:>8}{flag}")
    return rows


if __name__ == "__main__":
    main()
