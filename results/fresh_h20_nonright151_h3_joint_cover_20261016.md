# Two bounded wall banks fail to close the adjacent H3 face

The [exact158-wall chamber](fresh_h20_compact_quad_cap_20261016.md) has a
whole compact quadrilateral, but its adjacent three-face `S16+{151}` remains
open. Two numerical searches tested whether compatible collections of
new walls could make this entire three-face compact. Both finite banks
failed the complete compatible-cover test described below.

This is a **bounded numerical obstruction**, not a proof that the chamber
or the face cannot be completed. In particular, the banks impose three
finite old-side contacts on each new wall. Continuous families with only
two finite old-side contacts remain a constructive direction.

## Candidate families and complete face checks

Every candidate is perpendicular to S16. Its four prescribed contacts
are151 and three distinct walls from

```
31,95,132,142,145,150,152,156,157.
```

The three side labels range over `{2,3,4,5,6,8,10}`. The first bank requires
a right contact with151; the second allows labels2,3,4,5,6 there. For each
contact frame the source solves the norm-two equation in its remaining
one-dimensional kernel and tests both real branches. All158 old-wall
Coxeter angle conditions and the common-center condition are tested.
Duplicates and walls cutting none of the old exterior face vertices are
discarded.

| Bank | Prescribed patterns | Real branches | Retained roots |
|---|---:|---:|---:|
| Right151 | 28,812 | 57,218 | 1,730 |
| Labels2–6 at151 | 144,060 | 272,702 | 8,009 |

No individual candidate removes all old exterior face vertices. The128
strongest candidates in each bank were additionally checked by a full
three-dimensional halfspace intersection; none closes the face.

The candidate sources and data are
[right151 source](../scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py),
[right151 data](fresh_h20_quad_adjacent_h3_closure_20261016.json),
[nonright source](../scripts/fresh_h20_nonright151_h3_closure_20261016.py), and
[nonright data](fresh_h20_nonright151_h3_closure_20261016.json).
Both generation and fresh saved replay pass.

## Compatible collections also leave replacement escapes

The joint searches compute every mutual root relation. A finite clique
backtracking search then finds collections that cut the required sample
points and pass every mutual Coxeter angle test. Whenever such a
collection exists, the source computes its **whole three-face** again.
If an exterior vertex survives, a segment from an interior finite point
to that vertex supplies an actual feasible ideal-boundary sample. Any
compact completion must remove that point, so it is added to the next
cover problem.

This prevents removal of the old exterior vertices from being mistaken
for closure: replacement escape directions are explicitly found and
returned to the search.

| Bank | Numerically compatible pairs | Final samples | Final exhaustive search |
|---|---:|---:|---|
| Right151 | 871,021 | 20 | No covering clique;562 search nodes |
| Labels2–6 at151 | 21,740,931 | 21 | No covering clique;14,212 search nodes |

Both searches terminate by exhausting the finite graph, rather than by
their node or iteration guards. They start with eight numerical exterior
vertices; the other samples are feasible ideal points introduced after
complete face reconstruction. All compatible intermediate covers still
have exterior vertices. Generation and fresh saved replay pass.

The joint sources and data are
[right151 source](../scripts/fresh_h20_adjacent_h3_joint_cover_20261016.py),
[right151 data](fresh_h20_adjacent_h3_joint_cover_20261016.json),
[nonright source](../scripts/fresh_h20_nonright151_h3_joint_cover_20261016.py), and
[nonright data](fresh_h20_nonright151_h3_joint_cover_20261016.json).

The calculations use floating-point angle tolerances, a `1e-6` sphere
classification band, and a strict-cut threshold `1e-7`. Therefore the
finite graph obstruction is not promoted to an exact geometric
nonexistence certificate. No new walls are promoted by these searches,
and the exact158-wall chamber remains the constructive state.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_adjacent_h3_joint_cover_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_nonright151_h3_closure_20261016.py --verify
OPENBLAS_NUM_THREADS=1 python3 scripts/fresh_h20_nonright151_h3_joint_cover_20261016.py --verify
```
