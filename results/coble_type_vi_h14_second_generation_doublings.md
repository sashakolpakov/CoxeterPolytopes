# Seventeen second-generation `H14` divisor-lattice Coxeter chambers

**Status (2026-08-09): exact positive census theorem; worldwide priority
open.**

The two frozen Type-VI Coble doubles `H14/100` and `H14/101` each have 68
doubling walls.  Their intrinsic, undecorated Coxeter-diagram automorphism
groups have orders 16 and 12 and split those walls into respectively 7 and
10 orbits.  Doubling one representative from every orbit produces exactly
17 marked second-generation types.  Exact normalized-Gram refinement
separates all 17, so they are pairwise nonisometric as abstract Coxeter
chambers.

Every chamber is a finite-volume, noncompact polytope in `H14`, has index 4
and volume multiplier 4 relative to the original `H14/60` nef chamber, and
is redoublable again.  The facet census is

```text
3*H14/180 + 4*H14/181 + 4*H14/182 + 6*H14/183.
```

These are 17 new-to-project positive ledger entries.  They are not claimed
worldwide-new: Allcock already supplies infinite families in dimension 14,
and no publication-level isometry census has been completed here.

## Exhausting the second doubling choices

The normalized Coxeter graph deliberately forgets divisor square and
divisibility.  Thus its automorphisms are intrinsic to the abstract chamber
and cannot create a false distinction from the chosen surface realization.
The doubling-wall orbit representatives and resulting facet counts are:

| parent | orbit representatives, one-based | child facet counts |
|---|---|---|
| `H14/100` | `12,20,21,22,26,31,44` | `181,180,181,181,180,180,181` |
| `H14/101` | `12,16,17,20,21,25,31,32,39,45` | `183,183,183,182,183,182,182,182,183,183` |

The orbit sizes sum to all 68 walls in each parent.  An intrinsic symmetry
carrying one selected wall to another carries the doubled union with it, so
one representative per orbit is exhaustive for this generation.

For separation, the replay performs ten rounds of color refinement on the
complete edge-colored graph whose edge color is the exact normalized square

```text
(r.s)^2 / (r^2 s^2).
```

All facets start with one color and no root-length decoration is used.  The
17 final color-multiset hashes are distinct.  Isometric normalized Gram
systems would have the same refinement multiset, so distinct outputs prove
pairwise nonisometry within the 17-row census.  This does not compare the
rows with every chamber in the literature.

## Divisor-lattice construction

For a parent chamber `C` and selected doubling wall `h`, the child is

```text
C union s_h(C).
```

The construction is repeated in both the arithmetic lattice and the
physical Néron--Severi carrier `NS = I_(1,14)`.  In every row the selected
physical divisor has square `-1`; its reflection is an integral involutive
isometry.  All 180--183 child normals are primitive integral divisor
classes, span `NS` with Smith diagonal `1^15`, and have entrywise the same
normalized Gram as their arithmetic counterparts.  All wall reflections
are integral.

Exact transported facet witnesses prove irredundancy.  Each witness is zero
on its intended wall and has positive pairing with every other child wall;
the minimum positive slack in all 17 rows is 2.  Therefore the facet counts
are not merely counts of proposed inequalities.

The four repeated arithmetic profiles are:

| facets | primitive types `(-2,1);(-2,2);(-4,2)` | relations `m2;m3;m4;parallel;ultra` | doubling walls |
|---:|---:|---:|---:|
| 180 | `24;134;22` | `2968;89;244;1213;11596` | 134 |
| 181 | `26;134;21` | `2980;91;253;1211;11755` | 134 |
| 182 | `28;134;20` | `2992;93;262;1209;11915` | 134 |
| 183 | `30;134;19` | `3004;95;271;1207;12076` | 134 |

Every off-diagonal pairing has the Coxeter orientation.  Every finite-angle
normalized ratio is `0`, `1/4`, or `1/2`, corresponding to `m=2,3,4`; ratio
one is parallel and larger ratios are ultraparallel.

## Complete ray-transport lemma

The finite-volume proof depends on a complete extreme-ray transport, not on
a naive vertex count.  The required lemma is as follows.

**Union-ray lemma.** Let `D = C union s_h(C)` be a convex Coxeter double.
Every extreme ray of `D` off the glue is one of the two reflected copies of
a parent extreme ray.  A parent extreme ray on `h` survives once in `D` if
some other active wall is nonorthogonal to `h`; it is flattened if all its
other active walls are orthogonal to `h`.  There are no other extreme rays.

Indeed, away from `h` a neighborhood of the union lies in one parent copy.
On `h`, any candidate lies in the glue face, and an extreme ray of a face is
an extreme ray of the parent.  Delete the active normal `h`.  If every other
active normal is orthogonal to `h`, reflection fixes their span, which stays
one direction short.  If some active normal `a` is nonorthogonal, then

```text
s_h(a) - a
```

is a nonzero multiple of `h`; the two reflected active sets recover the
missing direction and hence full rank.  Finally, reflection is Lorentzian
and a glue ray is fixed, so retained rays keep their finite/ideal/hyperideal
type while flattening only removes rays.  This proves completeness and
causal-type preservation of the transported lists.

The four exact ray profiles are:

| facets | finite | ideal | hyperideal | parent on glue `(F,I)` | retained | flattened |
|---:|---:|---:|---:|---:|---:|---:|
| 180 | 107,264 | 3,925 | **0** | `(1052,61)` | `(732,61)` | `(320,0)` |
| 181 | 108,076 | 3,975 | **0** | `(280,11)` | `(0,11)` | `(280,0)` |
| 182 | 108,888 | 4,025 | **0** | `(1052,61)` | `(732,61)` | `(320,0)` |
| 183 | 109,700 | 4,075 | **0** | `(280,11)` | `(0,11)` | `(280,0)` |

## Literal Vinberg finite-volume test

The certificate does not infer finite volume merely from the presence of
cusps, nor merely from zero hyperideal rays.  On each complete transported
ray list it literally enumerates every rank-13 elliptic edge and checks that
it has exactly two finite-or-ideal endpoints.  At every ideal endpoint it
also checks that the active diagram is a union of affine components of total
parabolic rank 13 and is the full zero-neighborhood of each component.

| facets | rank-13 edges | finite-finite | finite-ideal | ideal-ideal | invalid cusps | bad edges |
|---:|---:|---:|---:|---:|---:|---:|
| 180 | 940,516 | 598,314 | 305,068 | 37,134 | **0** | **0** |
| 181 | 948,264 | 602,426 | 308,212 | 37,626 | **0** | **0** |
| 182 | 956,012 | 606,538 | 311,356 | 38,118 | **0** | **0** |
| 183 | 963,760 | 610,650 | 314,500 | 38,610 | **0** | **0** |

Across the 17 chambers this replays 1,847,848 finite vertices, 68,225 ideal
vertices and 16,221,212 elliptic rank-13 edges, with zero bad endpoint
counts, zero invalid affine neighborhoods, zero incomplete zero-neighborhoods
and zero hyperideal rays.  Hence all 17 chambers are finite-volume and their
ideal vertices make all 17 noncompact.

This endpoint audit depends explicitly on the independently complete seed
ray archive and the proved union-ray transport lemma.  It is not advertised
as a no-ray-import enumeration directly from the child diagrams.

## Redoubling frontier

Every child has 134 doubling walls and an explicitly archived parallel pair
among them.  Hence Allcock's redoubling mechanism applies again to every
row, providing further finite-volume `H14` chambers at higher reflection
index.  The present theorem classifies only the 17 symmetry-inequivalent
second-generation choices; it does not claim a collision-free census of all
later doubling words.

## Algebraic and novelty boundaries

What is proved is a positive `SURFACE -> POLYTOPE` divisor-lattice factory:

```text
complete Type-VI Coble nef chamber
  -> integral NS reflections
  -> 17 complete index-4 Coxeter chambers.
```

After crossing a nef wall, the reflected divisor classes are not asserted
effective on the original surface.  The child chambers are not asserted to
be nef cones or surface-automorphism fundamental domains.  This positive
result therefore does not supply the broader inverse
`POLYTOPE -> SURFACE -> algebraic classification -> theorem on P` leverage
sought on the theory track.

Likewise, `H14` was already populated infinitely before this census.  The
allowed novelty wording is **17 new-to-project, pairwise nonisometric within
the enumerated second generation; worldwide priority open**.

## Reproducibility

- [certifier](../scripts/certify_coble_type_vi_h14_second_generation_doublings.py)
- [certificate](../data/enumerations/coble_type_vi_h14_second_generation_doublings.json)
- [test](../tests/test_coble_type_vi_h14_second_generation_doublings.py)
- [first generation](coble_type_vi_h14_single_doublings.md)

Replay with

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_coble_type_vi_h14_second_generation_doublings.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_coble_type_vi_h14_second_generation_doublings.py
```

