# `gt_scan.py`

This script accompanies Proposition 1 of *Exceptional Sets for Semiprimes in Quadratic Intervals*. It checks arithmetic used in the claimed bound for the exceptional-set exponent at the square-root scale. Run it with `python3 gt_scan.py`; it uses only Python's standard library.

## What it prints

- The four lines labelled `sampled expression` evaluate the Gafni–Tao formula on a finite grid of values of σ for θ = 0.5, 0.499, 0.495 and 0.49. They help locate the apparent maximum near σ = 7/10. **A sampled maximum is not a proved upper bound**: the grid can miss a larger value between nodes or at an isolated admissible point.
- The `A(7/10)` and `A*(7/10)` line computes the tabulated bounds at that specific point using exact fractions. At θ = 1/2, the script obtains μ₂ = 97/130 and μ₄ = 183/260. It then displays 2(183/260) − 1 = 53/130.
- The final `GT point check` evaluates μ₄ exactly at θ = 17/30, σ = 7/10. The result, 7/12, matches Gafni–Tao's worked example. This is a check of the formula at one point, not a scan of all σ.

## Relationship to the proof

The paper establishes its upper bound by checking every applicable piecewise branch of the published bounds for A(σ) and A*(σ). The script is an independent diagnostic for transcription and arithmetic, not a substitute for that branch-by-branch argument. The displayed θ values below 1/2 illustrate the behavior used in the paper's limiting argument; their six-digit numerical output is not a certificate of an inequality for all sufficiently small perturbations.

## Why there is no grid line for θ = 17/30

At that θ, the relevant tabulated admissibility threshold is attained at σ = 7/10. A generic floating-point grid need not hit that point and can return the floor 1 − θ = 13/30 ≈ 0.433333, which misleadingly looks stronger than the actual calculation. The script now omits this scan and uses exact fractions for the worked-example check.
