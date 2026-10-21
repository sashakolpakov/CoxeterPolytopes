# Exact bounded incidence surgery on three noncompact H9 sources

No compact H9 candidate survives this batch. Two exact tests exclude **3,131,392 discrete finite-label assignments** and **151 right/dotted incidence-splitting patterns** on three previously unused 15/16-facet sources. Both tests leave the designated dotted distances completely unrestricted. These are bounded mutation-family results, not a nonexistence theorem for compact H9 polytopes or for every modification of these sources.

The sources are the archived finite-volume noncompact chambers in [the rank-ten source archive](../data/enumerations/alexeev_parked_rank10_h9_examples.json). All indices below and in the new artifacts refer to zero-based positions in each source's `positive_wall_gram`, not ambient root indices. This is a global Gram/incidence calculation; it does not extend the earlier vertex-centered cap banks or repeat the T1/T2 four-angle search.

| Source | Facets | Old cusps | Distinct affine components | Old parallel pairs | Necessary additional changed pairs | All minimum supports using old nonright finite pairs |
|---|---:|---:|---:|---:|---:|---:|
| A14_8_1_C02 | 16 | 23 | 55 | 5 | 8 | 42 |
| A14_10_1_C02 | 15 | 15 | 35 | 3 | 7 | 21 |
| A16_14_1_C01 | 15 | 9 | 23 | 6 | 5 | 34 |

Every original parallel pair is allowed to change without counting toward the “additional changed pairs” column. The same minima 8, 7 and 5 hold when the additional changes may also use originally right pairs. All minimum supports are enumerated only for the old nonright finite pairs; the larger all-pair calculation certifies the minimum and one witness, not an exhaustive support list.

The [affine-circuit certificate](fresh_h9_incidence_surgery_20261013.json) reconstructs the connected components of every saved ideal vertex and checks a strictly positive kernel by rational elimination. A compact acute Coxeter chamber retaining the facets cannot preserve such a component. Thus every remaining affine component must contain a changed pair. Once components containing an already freed parallel pair are removed, there are respectively 50, 32 and 17 constraints. Exact bitset set-cover recursion proves the minima. Replacing a dominated covering column by a superset preserves cardinality and coverage; the exhaustive support enumeration subsequently uses the original columns.

## Finite-label mutations: exact rank obstruction

For each minimum support consisting of old nonright finite pairs, replace every selected old label independently by a **different** member of `{2,3,4,5,6}`. Each selected old label is 3 or 4, hence there are four choices. Keep every other old finite/right entry fixed and free every old parallel entry. The complete tested counts are:

| Source | Assignments | Fixed 11-node principal subsets available | Certified rank at least 11 | Unresolved |
|---|---:|---:|---:|---:|
| A14_8_1_C02 | 42 × 4⁸ = 2,752,512 | 496 | 2,752,512 | 0 |
| A14_10_1_C02 | 21 × 4⁷ = 344,064 | 210 | 344,064 | 0 |
| A16_14_1_C01 | 34 × 4⁵ = 34,816 | 36 | 34,816 | 0 |

The [modular certificate](fresh_h9_incidence_surgery_20261013_modular.json) excludes every assignment by a nonzero 11×11 principal determinant **modulo 241**. Each tested principal subset avoids every freed pair, so its determinant is independent of all free distances. H9 normals lie in a ten-dimensional Lorentz space and therefore have Gram rank at most ten. A nonzero determinant of order eleven already rules out a realization, before compactness or detailed vertex combinatorics need to be tested.

After normalizing each wall to squared norm two, the finite entries are `0, −1, −√2, −(1+√5)/2, −√3`. They belong to `Z[1/2,√2,√3,√5]`. The prime 241 admits square roots of 2, 3 and 5, so reducing these generators defines a ring homomorphism to the finite field. An exactly zero determinant would reduce to zero. Thus a nonzero modular determinant proves nonvanishing in characteristic zero. The implementation also supports prime 1009, but no assignment needed it. All native integer products and elimination updates have an explicit bound below `2p² < 2³¹`.

The [wrapper](../scripts/fresh_h9_incidence_surgery_20261013_modular.py) binds its own hash, the [native exact enumerator](../scripts/fresh_h9_incidence_surgery_20261013_modular.cpp), the exact cover artifact, and the generated per-source input. Replay recompiles and rechecks every assignment; it is not a floating-point determinant test. The earlier [numerical Schur scout](fresh_h9_incidence_surgery_20261013_search.json) found no survivors, including its alternate-basis checks, but is not needed for this exact conclusion.

## Explicit incidence splitting: 151 all-distance exclusions

A separate branch replaces each selected old nonright pair by either a right pair or a free dotted pair, and again frees all old parallel pairs. This changes incidences rather than only shortening old angles. For each minimum finite-edge support, the selection rule chooses all smallest sets of new dotted pairs that meet every old cusp having no original parallel pair. Merely weakening the finite edges to right angles in such a cusp would otherwise leave too many positive directions.

There are 54, 25 and 72 patterns, respectively, for a total of **151**. This is the stated per-support minimum-dotted rule, not all possible right/dotted assignments or a classification of incidence changes.

The [separate Fraction checker](../scripts/fresh_h9_incidence_surgery_20261013_dotted_exact.py) reconstructs these 151 patterns independently of the floating-point proposal arithmetic. For each, it checks a fixed ten-node principal Gram matrix avoiding **all** free old and new dotted entries. Exact Schur elimination yields ten positive rational pivots in every case. Its positive index is therefore at least ten, whereas an H9 Gram matrix has positive index at most nine. The [saved exact certificates](fresh_h9_incidence_surgery_20261013_dotted_exact.json) rule out every real choice of the freed entries. No nonlinear distance solver was reached.

## What remains open

The finite-label theorem is restricted to the enumerated minimum supports of originally nonright finite pairs and labels 2 through 6. The right/dotted theorem is restricted to its 151 specified patterns. Larger changed supports, changes involving old right pairs, other labels, other dotted patterns, facet replacement, and other source chambers remain outside these results. The affine hitting constraints remain useful necessary conditions for such changed constructions. No numerical feasible Gram or new polytope is claimed.

Reproduce the exact results from the repository root:

```sh
python scripts/fresh_h9_incidence_surgery_20261013.py --verify
python scripts/fresh_h9_incidence_surgery_20261013_modular.py --verify
python scripts/fresh_h9_incidence_surgery_20261013_dotted_exact.py --verify
```

The first replay uses NetworkX; the modular replay needs a C++17 compiler. The final 151-pattern replay uses only the Python standard library. Each saved replay passed on the frozen files.

The [root review](fresh_h9_incidence_surgery_20261013_review.md) records
the mathematical/source audit and fresh modular and Fraction replays.
