# Five finite-volume `H9` polytopes from the parked Alexeev rank-ten stock

## Outcome

The archived ten-stock denominator contains exactly six proper subgroup
classes with fixed rank ten.  One, `(14,8,1)-C01`, is the already frozen
`T1 H9/14` chamber.  Every one of the other five gives a finite-volume
noncompact hyperbolic Coxeter polytope:

| ID | facets | finite + ideal + hyperideal vertices | elliptic edges | affine components | relation histogram |
|---|---:|---:|---:|---:|---|
| `A14_8_0_C01` | 12 | `44+5+0` | 260 | 10 | `54m2+8m3+4m4` |
| `A14_8_1_C02` | 16 | `165+23+0` | 1,112 | 55 | `96m2+9m3+10m4+5 parallel` |
| `A14_10_1_C01` | 13 | `70+7+0` | 393 | 14 | `63m2+9m3+6m4` |
| `A14_10_1_C02` | 15 | `123+15+0` | 778 | 35 | `85m2+8m3+9m4+3 parallel` |
| `A16_14_1_C01` | 15 | `100+9+0` | 588 | 23 | `85m2+8m3+6m4+6 parallel` |

The five labelled relation diagrams are pairwise nonisomorphic and are also
distinct from the three `H9/14` factory types and all nine doubling-cube
types already in the project ledger.  The 12-facet row is exactly Im Hof's
classical `9-imhof85_03` chamber; its positive contribution is the new
Alexeev/K3 divisor-lattice extraction.  The other four are unmatched in the
pinned CoxIter lanes and remain bounded literal-diagram priority candidates.
No worldwide novelty statement is part of this certificate.

The frozen replay is
[`certify_alexeev_parked_rank10_h9_examples.py`](../scripts/certify_alexeev_parked_rank10_h9_examples.py),
with data in
[`alexeev_parked_rank10_h9_examples.json`](../data/enumerations/alexeev_parked_rank10_h9_examples.json).

## Exact denominator

The input is the previously certified complete all-subgroup denominator for
ten Alexeev wall stocks below the old fixed-rank-eleven cutoff.  Its complete
fixed-rank-ten list is

```text
(14,8,0)-C01
(14,8,1)-C01, C02
(14,10,1)-C01, C02
(16,14,1)-C01.
```

Thus the present batch plus the earlier `T1=(14,8,1)-C01` closes this exact
parked denominator positively: six rows, six finite-volume noncompact `H9`
sections.  This is not an exhaustive census of all rank-ten K3 actions or all
hyperbolic `H9` polytopes.

The `(16,14,1)` row exposes a useful implementation trap.  The vector pairing
one with each leading basis root is spacelike for this particular source
basis.  That does not make the chamber section empty.  The sum of all 29
source roots has square `552`, pairs with every source wall by at least `2`,
and its subgroup average is a strict fixed timelike interior.  The certifier
detects this basis issue exactly and uses the root sum.

## Why finite volume is proved

For every row, exact fixed-lattice reconstruction and conic reduction give a
full-dimensional pointed cone.  Every retained primitive mirror has an
integral reflection on the saturated rank-ten fixed lattice, and every
normalized Gram product is Coxeter or parallel.

GMP rational double description enumerates every extreme ray.  Every ray is
on the chosen time sheet, every Lorentz square is nonnegative, and every row
has isotropic but no hyperideal rays.  This proves finite volume and
noncompactness.

The literal Vinberg replay is independent of the ray-sign conclusion: every
elliptic rank-eight edge has exactly two finite/ideal endpoints, and every
connected affine component has a complete cusp zero-neighbourhood.  Both bad
counts are zero in all five rows.

## Divisor-lattice meaning and scope

Each fixed lattice is a saturated primitive Lorentzian sublattice of an
Alexeev K3-admissible Picard lattice.  A generic marked K3 period therefore
gives a weak full-Picard divisor carrier for the oriented mirror system.  The
source/folded norm data are:

| ID | fixed determinant / Smith | active source norms | folded primitive norms |
|---|---|---|---|
| `A14_8_0_C01` | `-1024`, `1^2 2^6 4^2` | `5*(-2)+7*(-4)` | `3*(-2)+6*(-4)+3*(-8)` |
| `A14_8_1_C02` | `-1024`, `1^2 2^6 4^2` | `9*(-2)+7*(-4)` | `6*(-2)+7*(-4)+3*(-8)` |
| `A14_10_1_C01` | `-4096`, `2^8 4^2` | `3*(-2)+10*(-4)` | `10*(-4)+3*(-8)` |
| `A14_10_1_C02` | `-4096`, `2^8 4^2` | `6*(-2)+9*(-4)` | `4*(-2)+7*(-4)+4*(-8)` |
| `A16_14_1_C01` | `-16384`, `2^6 4^4` | `4*(-2)+11*(-4)` | `4*(-2)+4*(-4)+7*(-8)` |

These are all weak divisor carriers.  A later exact Picard completion extends
the first four source actions to genuine symplectic Nikulin K3
automorphisms, preserving the archived fixed sections.  That promotion still
does not make the walls effective or any chamber a full nef or quotient-nef
cone; see the
[`parked K3-action theorem`](alexeev_parked_h9_k3_action_carriers.md).

## Lineage boundary

At pinned CoxIter commit `b800be48`, `A14_8_0_C01` has one exact
label-preserving graph match:

```text
CoxIter-master/graphs/doubly truncated orthoscheme/9-imhof85_03.coxiter
candidate -> Im Hof wall permutation: 6,11,10,9,8,1,2,3,4,5,7,12.
```

There are 21 pinned `H9/12` records; three share the global relation
histogram, but only `_03` is isomorphic.  The bounded same-header comparison
for the remaining rows gives eight records at `H9/13`, eight at `H9/15`, and
none at `H9/16`.  Both 15-facet rows and the 13-facet row have zero relation-
histogram collision and zero exact graph match.  The 16-facet row has no
same-header record.

This proves a classical identification for the first row and a bounded
corpus exclusion for the other four.  It does not prove worldwide priority:
the relevant noncompact nonsimple few-facet classifications are incomplete,
and the CoxIter bundle is not a publication-complete database.

A separate hostile replay reconstructs all five ray sets and Vinberg audits
without calling this batch builder.  It also exhausts every
reflection-stable even discriminant overclosure.  The unique maximal indices
are respectively `4,8,4,8,32`, giving the known 2-elementary classes
`(10,6,0)`, `(10,4,0)`, `(10,8,0)`, `(10,6,0)`, and `(10,4,0)`.
After the recorded sign/scaled-dual normalizations they lie in
Scharlau--Kirschmer `List9` entries 23 or 7.  Thus there are no new maximal
arithmetic reflection classes here.  See the
[`H9` hostile lineage audit](alexeev_parked_rank10_h9_examples_hostile_lineage.md).

The primary-source companion audit identifies those three maximal carrier
types integrally as
`U+D4(-1)^2`, `U(2)+D4(-1)^2`, and `U+E8(-2)`.  If used as full invariant
lattices of nonsymplectic K3 involutions, they are the classical
ten-dimensional types with fixed loci respectively
`C4+3P1`, `C3+2P1`, and two elliptic curves.  This is an algebro-geometric
interpretation of the arithmetic envelope, not a promotion of the five
original chamber actions.  Independently, four original actions are now
promoted by the separate `E8(-2)` Nikulin completion.  See the
[`H9` overclosure/K3 lineage theorem](alexeev_parked_rank10_h9_lineage_overclosures.md).

Four rows have further doubling walls:

```text
A14_8_1_C02:  7, 15, 16
A14_10_1_C01: 3, 5, 11
A14_10_1_C02: 1, 3, 6, 15
A16_14_1_C01: 1, 8, 14, 15.
```

That queue has now closed positively.  The three orthogonal walls of
`A14_10_1_C01` give three further finite diagram types.  The two rows with one
parallel pair give infinite path-times-cube families, and the two parallel
pairs of `A16_14_1_C01` give a rectangular two-parameter family with
`mn+m+n+12` facets and volume multiplier `mn`.  Six representatives pass
fresh GMP/time-sheet/Vinberg replays.  See the
[`H9` product-redoubling theorem](alexeev_parked_h9_doubling_products.md).

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_alexeev_parked_rank10_h9_examples.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_alexeev_parked_rank10_h9_examples.py
```
