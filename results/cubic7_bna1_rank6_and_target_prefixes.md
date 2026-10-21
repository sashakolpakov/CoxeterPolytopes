# Cubic-seven non-unimodular rank-six and target-dimension prefixes

This record freezes an exact discovery search over

\[
 L_{d,t}=\operatorname{diag}(-\theta,1^{d-1},t),\qquad
 \theta=2\cos(2\pi/7),\quad
 t\in\{2,2-\theta\}.
\]

The positive result is an exact, dimension-parametric Vinberg-prefix engine;
the search has **not** found a new finite-volume Coxeter polytope in this
family at the archived bound.  A hyperideal ray in a prefix is not a proof of
nonreflectivity, because a later mirror can cut it off.

## Exhaustive arithmetic input

For a primitive vector `r`, reflection integrality is
`q(r) | 2 I(r)`, with
`I(r)=(-theta*r0,r1,...,r_(d-1),t*r_d)`.  The local squarefree argument gives
the complete primitive mirror-norm classes, modulo totally-positive unit
squares:

- `t=2`: `1,2`;
- `t=2-theta`: `1,2,2-theta,2(2-theta)`.

Height zero is exhausted exactly.  Modulo signed permutations of the
`d-1` ordinary coordinates it has three patterns: a coordinate root, a
signed-pair root, and the final-coordinate root.  Their expanded count is
`2(d-1)^2+2`, precisely the root system `B_(d-1) x A1`: 34 roots in H5,
884 in H22, and 1570 in H29.

Positive-height ordinary coordinates are solved as ordered square partitions,
never as all tuples or coordinate permutations.  At each exact-height level,
all candidates are tested only against lower levels; the surviving level is
then proved mutually acute before it is appended.  Thus no arbitrary
`(norm,vector)` tie-break chooses the chamber.

## Exact height-five results

| dimension | seed | pool | accepted | walls | finite rays | cusps | hyperideal rays |
|---:|:---:|---:|---:|---:|---:|---:|---:|
| 5 | L2 | 4 | 1 | 6 | 3 | 0 | 3 |
| 5 | L7 | 5 | 2 | 7 | 5 | 0 | 5 |
| 22 | L2 | 4 | 1 | 23 | 3 | 0 | 20 |
| 22 | L7 | 5 | 2 | 24 | 5 | 0 | 39 |
| 29 | L2 | 4 | 1 | 30 | 3 | 0 | 27 |
| 29 | L7 | 5 | 2 | 31 | 5 | 0 | 53 |

All current cones are irredundant and their nonspacelike rays lie on a
consistent time-cone sheet, but every one still has hyperideal rays.  No
finite-volume stopping event occurs.

For every dimension 22 through 29 the full archived data give the same root
pool sizes and accepted-positive-wall counts.  The hyperideal counts are
`d-2` for L2 and `2d-5` for L7 at this height.  This is a certified bounded
stabilization pattern, not an unbounded theorem.

## Adaptive H5 continuation

The rank-six/H5 engine was also run without geometry stopping at larger
heights.  It did not stabilize:

| height | L2 pool / accepted / walls | L7 pool / accepted / walls |
|---:|:---:|:---:|
| 50 | 29 / 9 / 14 | 36 / 7 / 12 |
| 100 | 80 / 18 / 23 | 94 / 15 / 20 |
| 150 | 134 / 30 / 35 | 156 / 23 / 28 |

This sustained wall growth is evidence worth pursuing, but it is not itself
a nonreflectivity certificate.  The obvious smallest Pell isometries do not
supply one: exact folding places them back in the reflection group (the L2
fold word in the rank-five control is `[7,0,1,2,1,0,7]`; the L7 candidate
likewise folds to the identity).

## Reproduction and scope

- Engine: `src/coxeter_polytopes/search/cubic7_bna1.py`
- Generator: `scripts/search_cubic7_bna1_prefixes.py`
- Tests: `tests/test_cubic7_bna1_search.py`
- Full H22--H29 archive: `data/enumerations/cubic7_bna1_h22_h29_height5.json`
- Compact H5/H22/H29 archive: `data/enumerations/cubic7_bna1_h5_h22_h29_height5.json`

Run:

```console
PYTHONPATH=src python3 scripts/search_cubic7_bna1_prefixes.py \
  --dimensions 22-29 --height 5 --seeds both \
  --output data/enumerations/cubic7_bna1_h22_h29_height5.json
PYTHONPATH=src pytest -q tests/test_cubic7_bna1_search.py
```

The archived claim is only: complete primitive norm classes, complete root
enumeration through height five, exact acute-wall selection, and exhaustive
extreme-ray geometry of those prefixes.  It does not claim termination,
finite covolume, compactness, or nonreflectivity.

