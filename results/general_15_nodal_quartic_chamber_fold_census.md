# General 15-nodal quartic: exact chamber-fold census

## Outcome

There is a complete negative result for this algebraic stock above dimension
9:

> **Theorem (computer-assisted, exact).** No fixed-subspace fold of the
> Dolgachev--Shimada induced chamber `D15` by a nontrivial subgroup of
> `O(S15,D15) ~= S6` is a hyperbolic Coxeter polytope in dimension greater
> than 9.

The only possible subgroup in that range is a transposition subgroup.  Its
fixed space is `H^11`; the restricted chamber has 128 genuine irredundant
facets, but 82 facet pairs have forbidden normalized squared products.  Thus
it fails the Coxeter-angle gate before compactness or finite volume is at
issue.

This is a stock-specific no-go theorem, not a general non-existence theorem
for hyperbolic Coxeter polytopes above dimension 9.

## Exact algebraic input

The official data of Dolgachev and Shimada give:

- the rank-16 lattice `S15`, with determinant `-128`, inertia `(1,15)`, and
  ample class of square `17/2`;
- the complete list of 314 walls of the induced chamber `D15`; and
- all 720 matrices in `O(S15,D15)`, with their natural `S6` labels.

The verifier replays every published wall norm and ample pairing.  It then
checks that every one of the 720 matrices preserves the Gram matrix, fixes
the ample class, permutes all 314 walls, and obeys the published `S6` group
law on the faithful wall action.  No surface-automorphism realization is
used or required.

Sources: [paper](https://arxiv.org/abs/1908.05390),
[official exact data](https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/15NodalQuart/15nodal/compdata/15nodalcompdata.txt),
and [official data dictionary](https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdatatex/15nodalCompData.pdf).

## Why the high-dimensional census is exhaustive

For a subgroup `H <= S6` and any `g in H`,

```text
Fix(H) is contained in Fix(g).
```

The exact fixed ranks of the nonidentity `S6` conjugacy types in this
rank-16 representation are:

| cycle type | fixed rank | hyperbolic dimension |
|---|---:|---:|
| `2 1^4` | 12 | 11 |
| `2^2 1^2` | 10 | 9 |
| `2^3` | 10 | 9 |
| `3 1^3` | 8 | 7 |
| `3^2` | 6 | 5 |
| `3 2 1` | 6 | 5 |
| `4 1^2` | 6 | 5 |
| `4 2` | 6 | 5 |
| `5 1` | 4 | 3 |
| `6` | 4 | 3 |

Hence a subgroup with fixed hyperbolic dimension greater than 9 can contain
only the identity and transpositions.  Two distinct transpositions multiply
to either a 3-cycle or a double transposition, both of fixed rank at most 10.
Such a subgroup therefore contains only one nonidentity element and is
exactly a transposition `C2`.  All transpositions are conjugate, so one exact
screen settles every nontrivial subgroup in the target range.

## Fold census

Each row uses exact orbit-sum normals, exact conic redundancy certificates,
strict timelike relative-interior witnesses for every retained facet, and
exact normalized products.

| generator cycle type | dimension | wall orbits | facets | forbidden pairs | status |
|---|---:|---:|---:|---:|---|
| `2 1^4` | 11 | 178 | 128 | 82 | not Coxeter |
| `2^2 1^2` | 9 | 176 | 68 | 143 | not Coxeter |
| `2^3` | 9 | 172 | 66 | 0 | Coxeter; volume not yet audited |
| `3 1^3` | 7 | 114 | 31 | 15 | not Coxeter |
| `3^2` | 5 | 112 | 21 | 0 | Coxeter; volume not yet audited |
| `5 1` | 3 | 66 | 9 | 0 | finite volume, noncompact |

For the decisive `H^11` row, the 82 forbidden ratios are `1/3` (44 pairs),
`2/3` (32 pairs), and `1/9` (6 pairs).  None is a squared cosine of an
allowed rational Coxeter angle.

The `H^3` row is a literal Vinberg-criterion control: it has 19 elliptic
edges, every edge has exactly two finite/ideal endpoints, 10 finite vertices,
and 2 ideal vertices.  It is therefore an exact finite-volume noncompact
Coxeter polyhedron.  No novelty claim is made for this low-dimensional
example.

The `H^9` and `H^5` rows pass only the Coxeter-angle gate.  Their Vinberg
edge/end-point and full parabolic-neighbourhood audits are explicitly
deferred, so the artifact makes **no finite-volume claim** for either one.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/search_15_nodal_quartic_chamber_folds.py
PYTHONPATH=src:. python3 scripts/search_15_nodal_quartic_chamber_folds.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_15_nodal_quartic_chamber_folds.py
```

Frozen outputs:

- `data/known_examples/general_15_nodal_quartic_chamber_stock.json`
- `data/enumerations/general_15_nodal_quartic_prime_order_folds.json`
