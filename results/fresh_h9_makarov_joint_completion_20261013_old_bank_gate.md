# The previous cap bank cannot complete the retained 27-wall chamber

This is an exact conditional obstruction, not a claim about arbitrary completion. Every one of the **25,101 literal recipes** in the two previous numerical banks either misses a common strict spacelike escape or has a forbidden pairing with one of the sixteen retained caps. Therefore no acute-compatible selection from those banks can produce finite volume while all 27 facets are retained.

The literal count includes duplicates; the earlier 23,841 count was numerical deduplication, which is not used in this proof.

- [Certificate](fresh_h9_makarov_joint_completion_20261013_old_bank_gate.json).
- [Exact interval reconstruction and replay](../scripts/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.py).
- [Retained 27-facet chamber](fresh_h9_makarov_joint_20261013.md).

| Exact gate | Literal recipes |
|---|---:|
| Strictly negative escape pairing | 24,152 |
| Positive pairing with a retained cap | 553 |
| Forbidden finite-angle interval with a retained cap | 396 |
| Total | 25,101 |

The point is

    (153240957355, 189417557874, 130887763098, 72360464506,
     33886289902, 10763857490, -240266086393, -296985201824,
     -205211687438, -113438150973).

Its square is positive; its center pairing and all 27 wall pairings are strictly negative. The exact base-field values and interval signs are saved in the certificate.

## What the interval proof establishes

Each archived proposal is interpreted as an exact normal defined by nine prescribed pairings and the norm-two equation. The nine-by-ten system is solved over the base field; the prescribed finite cosines and any normalization square roots have rigorous rational intervals. The quadratic has a strictly positive discriminant. Precisely one branch lies inside the rational coordinate box of radius \(10^{-5}\) around the archived proposal, and every other branch is separated from that box in at least one coordinate.

All subsequent enclosures use integer arithmetic with outward rounding at scale \(10^{70}\). Norm two follows from the quadratic identity, not from a small numerical residual. The numerical proposal only selects the branch.

For a rejected pair, both normals have norm two. Finite Coxeter pairings must equal \(-2\cos(\pi/m)\). Rational Machin-formula and Taylor bounds enclose all labels from 2 through 120. The entire remaining interval down to -2 is conservatively retained, so a finite-angle rejection excludes **every** finite label without a field-degree assumption. Dotted and parallel relations are not incorrectly rejected.

If an entry has a positive pairing with a retained facet, it cannot be added as a distinct acute-compatible wall. Positive-proportional copies of existing facets cause no loophole: they already have negative pairing with the saved escape and fall under that first gate.

For every genuinely compatible entry, the escape pairing is strictly negative. Thus any finite choice of such entries leaves the point strict and spacelike. Replacements of retained facets, other roots, other branches outside the certified boxes, and new recipes are not excluded.

## Reproduction

    python3 scripts/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.py --verify

The replay reconstructs all 25,101 exact definitions and all rejection witnesses. It binds both banks, the exact retained chamber, and the reused interval and trigonometric implementations by hash. The next constructive search uses profiles at the 27-wall chamber's new finite vertices, rather than extending this closed bank.
