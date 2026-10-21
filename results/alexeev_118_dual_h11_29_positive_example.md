# A finite-volume 29-wall Coxeter polytope in H11

## Theorem

Alexeev's complete 118-wall chambers for the even hyperbolic 2-elementary
lattices `(16,6,1)` and `(16,12,1)` each have exact weighted wall-diagram
automorphism group `S6`.  Exhausting all 56 subgroup conjugacy classes in
each stock leaves exactly one proper class with fixed rank at least 11: the
class of transpositions.  Its fixed section is a strict integral Coxeter
polytope in `H11` with 29 facets.

The two stocks produce one polytope, not two.  Their complete normalized
signed 29-wall Grams agree under an explicit facet permutation.  More
strongly, after taking the unique index-8 reflection-stable overclosure on
the `(16,6,1)` side and globally dividing the `(16,12,1)` fixed Gram by `2`,
an explicit determinant-one integral matrix identifies the two lattices and
conjugates all 29 reflections wall by wall.  Both normalized models are
integrally isometric to Scharlau--Kirschmer `List11/Gram11` entry 5,
`I_(11,1)(2^2)`.

The polytope has finite volume and is noncompact.  The exact Vinberg audit
finds

| invariant | value |
|---|---:|
| facets | 29 |
| elliptic rank-10 edges | 28,798 |
| finite vertices | 3,876 |
| maximal rank-10 parabolic cusps | 253 |
| bad elliptic edges | 0 |
| bad affine zero-neighborhoods | 0 |

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
affine component has a full zero-neighborhood of the required parabolic
rank.  Thus the conclusion uses Vinberg's finite-volume criterion, not the
mere existence of a finite wall list.

## What is new, and what is not

This is a newly extracted literal 29-wall chamber under the frozen
comparison used in this project.  No `H11/29` chamber occurs in the pinned
ordinary or doubling CoxIter corpora, in Alexeev's 18 explicit *full*
chamber matrices, in the earlier positive ledger, or among the 27
Scharlau--Kirschmer `List11` *full* chambers.

That bounded absence is not a worldwide priority theorem.  Allcock proved
that there are infinitely many finite-covolume Coxeter groups in every
dimension through 19, so this is neither a first example in `H11` nor a new
dimension.  See [Allcock's paper](https://arxiv.org/abs/0903.0138).

It is also not a new arithmetic commensurability class.  The common
normalized lattice is exactly Scharlau--Kirschmer `List11` entry 5.  The
[Scharlau--Kirschmer census page](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
states that each listed strongly squarefree genus has class number one; the
entry-5 full reflection chamber has 14 walls, 56 finite vertices, and 5
cusps.  The 29-wall group is a finite-index reflection subgroup of that full
group.  Its numerical index remains open because Coxeter Euler
characteristics vanish in odd hyperbolic dimension.

The defensible novelty statement is therefore:

> one apparently previously unrecorded literal finite-volume noncompact
> `H11/29` Coxeter chamber, extracted from Alexeev's divisor-lattice-inspired
> stock and now placed exactly inside a known maximal reflective lattice.

## Two algebraic carriers now known

The chamber now has two arithmetically different divisor-lattice
realizations.  The `(16,12,1)` side has an even stable overclosure
`U+D8(-1)+A1(-1)^2`, the exact very-general Picard lattice of the explicit
BH `0.4.4.2` complex-K3 family.  An exact 87-reflection Weyl framing makes
all seventeen `(-2)` mirrors effective (one irreducible and sixteen
reducible) and gives a primitive square-44258 ample class inside all 29
walls.  On the `(16,6,1)` side, the doubled dual of
the unique odd index-eight closure is explicitly
`U+E8(-2)+<-2>^2`, the actual saturated involution-fixed divisor lattice on
the supersingular Fermat quartic in characteristic three.  All 29 oriented
walls and reflections transport in both cases and generate the Fermat target
with index one.  Neither arrangement is asserted to be the complete nef cone
or automorphism-fixed.  See the [complex-K3 carrier](alexeev_h11_29_ud8_divisor_overlattice_carrier.md),
its [nef/effectivity strengthening](a16_h11_29_ud8a1x2_nef_effectivity.md),
and the [Fermat scaled-dual carrier](a16_h11_29_fermat_char3_scaled_dual_carrier.md).

## The two lattice normalizations

The distinction between overclosure and global rescaling matters.

| stock | original fixed lattice | actual stable overclosure | normalization used for the common system |
|---|---|---|---|
| `(16,6,1)` | determinant `-256`, Smith `1^6 2^4 4^2` | unique index 8, odd determinant `-4`, Smith `1^10 2^2` | the stable overclosure itself |
| `(16,12,1)` | determinant `-16,384`, Smith `2^10 4^2` | unique index 8, odd determinant `-256`, Smith `1^4 2^8` | the original Gram divided globally by `2`, giving odd determinant `-4`, Smith `1^10 2^2` |

For each original lattice, the full discriminant module was exhausted.  In
both cases the reflection-stable integral subgroup orders are

`[1,2,4,8]`,

and the stable even subgroup orders are

`[1,4]`.

The common isometry satisfies

`A^T (G_(16,12,1)/2) A = M_(16,6,1)`

with `det(A)=1`.  Under the already certified facet permutation it maps all
29 oriented roots exactly after the forced positive root rescalings, and it
conjugates every corresponding reflection matrix.  Hence “same normalized
integral reflection system” is literal here.

## The misleading determinant collision

The original `(16,6,1)` fixed lattice and the separate
`A18_4_0_C03` `H11/28` fixed lattice both have rank 12, determinant `-256`,
and Smith invariants `1^6 2^4 4^2`.  They are nevertheless not integrally
isometric.  Their discriminant quadratic forms already differ: the latter
has classes of square `3/4` and `7/4` modulo `2`, while the former does not.
Their maximal reflection-stable overclosures also have different
determinants, `-16` and `-4` respectively.

Thus determinant plus Smith data were only a coarse header, not a lineage
identification.

## Algebro-geometric content and limit

Alexeev's work computes large reflective 2-elementary lattices and applies
them to K3 and hyperkahler geometry; see
[Alexeev's paper](https://arxiv.org/abs/2209.09110).  The `(16,6,1)` lattice
lies on the `r+a=22` K3 line and can occur as the Néron--Severi lattice of a
complex 2-elementary K3 surface.  This gives the construction a genuine
divisor-lattice carrier rather than an arbitrary Lorentzian vector space.

This original theorem is a lattice-and-diagram theorem: it proves that the
transposition acts integrally and that its fixed chamber is cofinite.  A
subsequent exact overclosure calculation now supplies the missing explicit
surface carrier.  The determinant-256 wall lattice has a unique nontrivial
reflection-stable even overgroup, of index four, and an explicit isometry
identifies it with `U+D8(-1)+A1(-1)^2`, the very-general Picard lattice of
the BH `0.4.4.2` Weierstrass family.  All 29 reflections extend integrally;
see the
[`H11/29` divisor-overlattice theorem](alexeev_h11_29_ud8_divisor_overlattice_carrier.md).

The automorphism boundary remains unchanged.  The row-40 Picard involution
preserves the wall lattice but not the selected 29-wall chamber, so the new
result is a mixed divisor-reflection subsystem rather than an
automorphism-fixed or quotient-nef chamber.  The `(16,12,1)` stock remains
the exact dual-lattice realization, not an additional K3 Picard-lattice
claim.

## Reproduction

The permanent certificates are:

- [construction and complete subgroup denominator](../data/enumerations/alexeev_118_dual_h11_29_construction.json), generated by
  [certify_alexeev_118_dual_h11_29.py](../scripts/certify_alexeev_118_dual_h11_29.py);
- [literal Vinberg finite-volume certificate](../data/enumerations/alexeev_118_dual_h11_29_vinberg.json), generated by
  [certify_alexeev_118_dual_h11_29_vinberg.py](../scripts/certify_alexeev_118_dual_h11_29_vinberg.py);
- [bounded priority and source-lineage certificate](../data/enumerations/alexeev_118_dual_h11_29_lineage.json), generated by
  [certify_alexeev_118_dual_h11_29_lineage.py](../scripts/certify_alexeev_118_dual_h11_29_lineage.py);
- [stable-overlattice and common-integral-system certificate](../data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json), generated by
  [certify_alexeev_118_dual_h11_29_stable_overlattices.py](../scripts/certify_alexeev_118_dual_h11_29_stable_overlattices.py).

Run from the repository root with `PYTHONPATH=src:.`:

```bash
python3 scripts/certify_alexeev_118_dual_h11_29.py --check-archive
python3 scripts/certify_alexeev_118_dual_h11_29_vinberg.py --check-archive
python3 scripts/certify_alexeev_118_dual_h11_29_lineage.py --check-archive
python3 scripts/certify_alexeev_118_dual_h11_29_stable_overlattices.py --check-archive
pytest -q tests/test_alexeev_118_dual_h11_29_stable_overlattices.py
```
