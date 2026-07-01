#!/usr/bin/env python3
"""Produce the report + semi_results.npz from a checkpoint, without re-running.
Usage:  python finalize.py 1e11        # HI MUST match the run you started
Works on a partial checkpoint too (gives results up to the point reached)."""
import sys, pickle, importlib.util
HI_ARG = sys.argv[1] if len(sys.argv) > 1 else "1e10"
sys.argv = ["sp", HI_ARG]                       # feed HI to the module's CLI parse
spec = importlib.util.spec_from_file_location("sp", "semiprime_production.py")
sp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sp)                     # runs top-level (sets HI/NUM_DEC/NMAX), NOT __main__
with open("semi_ckpt.pkl", "rb") as f:
    st = pickle.load(f)
print(f"loaded checkpoint: reached A={st['next_A']} ({st['seg_count']} segments)")
sp.report(st, sp.base_prime_sieve(sp.math.isqrt(sp.HI) + 1))
