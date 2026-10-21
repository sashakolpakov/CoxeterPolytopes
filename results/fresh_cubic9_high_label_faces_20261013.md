# Large dihedral labels and a necessary lower-dimensional test

The primary sources checked do **not** rule out a compact H9 Coxeter polytope with an `I₂(18)^4 × A1` vertex. They do imply that such a polytope has nested compact Coxeter faces in dimensions 7, 5, 3 and 1. For the direct-sum cubic-nine hosts, these faces use mirrors whose removed binary-block coordinates vanish. This is a useful necessary test of the construction, not a conclusion from the currently available bounded searches.

This test has subsequently been completed for α=Tβ: the [full H3 chamber-stabilizer certificate](fresh_cubic9_global_nonreflectivity_20261013.md) now rules out finite-covolume integral reflection completions for that direct-sum H9 host. The literature review and other-coefficient scope below are unchanged.

## What the angle bounds actually say

Bredon's Proposition 3.1 bounds finite-volume dihedral labels by `m≤6` in dimensions `n≥32`. Proposition 3.2 gives the same bound in `n≥7` only under the additional assumption that all facets mutually intersect. Proposition 3.3 treats ideal polytopes in `n≥4`. None applies to a general compact H9 polytope with disjoint facets. Indeed, her Theorem 2.1 recalls the stronger Felikson–Tumarkin result that compact Coxeter polytopes with mutually intersecting facets have dimension at most four. The label restriction with that hypothesis therefore supplies no new compact-H9 obstruction. [Bredon, *On ADEG-polyhedra in hyperbolic spaces*, §2–3](https://arxiv.org/pdf/2507.05153).

Results with a small number of facets should also be kept separate. In particular, nonexistence of compact `d+4`-facet polytopes for `d≥8` is a facet-count restriction, not a restriction on every high-dimensional diagram containing an edge labelled 18. [Felikson–Tumarkin, *On compact hyperbolic Coxeter d-polytopes with d+4 facets*](https://arxiv.org/pdf/math/0510238).

## The face criterion and its consequence

Allcock's Theorem 2.1, attributed there to Borcherds, says that a face corresponding to a spherical subdiagram with no `A_n` or `D5` component is itself Coxeter. Theorem 2.2(1) says that angles are unchanged when the two remaining nodes do not attach to that subdiagram. An `I₂(18)` component satisfies the hypothesis. No connected spherical diagram of rank at least three contains a bond labelled 18, so every node that extends this pair spherically is orthogonal to both of its nodes. Thus its face diagram is the unchanged subdiagram on those orthogonal nodes. [Allcock, *Infinitely many hyperbolic Coxeter groups through dimension 19*, Theorems 2.1–2.2](https://arxiv.org/pdf/0903.0138).

Consequently a compact H9 polytope with the stated vertex has the following necessary sequence. Compactness is inherited by faces.

| Removed binary blocks | Face dimension | Vertex diagram in that face |
|---|---:|---|
| 1 | 7 | `I₂(18)^3 × A1` |
| 2 | 5 | `I₂(18)^2 × A1` |
| 3 | 3 | `I₂(18) × A1` |
| 4 | 1 | `A1` |

This consequence is our application of the cited face theorem. It is not a classification of any of these lower-dimensional faces.

## Exact specialization to the searched lattices

Let `Lα = <-α> ⊕ B^4 ⊕ <1>`, where

```
B = [[2, -T], [-T, 2]],
T³ - 3T - 1 = 0,  T = 2cos(π/9),
D = 2 + T - T²,  β = 4 - T² > 0.
```

In one block the two simple zero-height roots are `a=-e₁` and `b=D(e₁-e₂)`. For any candidate face wall `r`, the necessary equations `(r,a)=(r,b)=0` imply `(r,e₁)=(r,e₂)=0`, because `D≠0`. Since `det(B)=β≠0`, both coordinates of `r` in that block vanish. No projection or change of its norm is required. Repeating the argument removes any chosen collection of blocks.

The orthogonal complementary lattice is exactly

```
Lα,k = <-α> ⊕ B^k ⊕ <1>,  k = 3, 2, 1, 0,
```

with the same coefficient `α`. Every integral reflection in this complementary lattice extends by the identity on the removed blocks to an integral reflection of `Lα`. Conversely, a host reflection orthogonal to those blocks restricts to an integral reflection of the complementary lattice. These are direct-sum identities, independent of any numerical search.

If the full H9 integral reflection chamber were compact, the face just described would supply a compact reflection polytope inside each corresponding lower-dimensional integral reflection arrangement. Its reflection group is a subgroup of the full lower-dimensional integral reflection group. The latter is discrete; containing a compact-lattice subgroup forces finite index and a compact quotient. Hence **a proved nonreflectivity result for one of these full lower hosts would obstruct a compact H9 chamber in this particular host**. A bounded lower-host prefix with a spacelike escape does not prove that premise: unenumerated mirrors may still close it.

The H3 and H5 cases (`k=1,2`) are supported directly by the present producer and are cheaper necessary tests. The H7 case requires adding the three-block join, which the current producer does not implement. These observations do not claim that an arbitrary compact H9 polytope must arise in one of our arithmetic hosts.

## Review of the cubic-nine report

I read the [cubic-nine report](fresh_cubic9_dihedral_report_20261013.md), its independent source review, and the producer. I found no mathematical correction needed in the integral-basis/class-number argument, unit-sign normalization, dual-exponent norm pool, elimination of primitive norms 4 and 4β, complete zero-height chamber, or distinction between bounded infinite-volume prefixes and full-host nonreflectivity. The angle alphabet is consistent with the real conductor-nine field. The fourth coefficient `2Tβ` has the same safe norm pool as `2T` because its dual exponent is `2β`.

This is an additional source-level review and literature check. It does not repeat the full height-300 spatial enumeration or replace the archived Fraction and enumeration replays linked by that report.
