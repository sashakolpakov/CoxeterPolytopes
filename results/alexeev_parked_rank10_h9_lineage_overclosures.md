# Exact lineage and reflection-stable overclosures of five parked H9 seeds

## Verdict

The five finite-volume noncompact H9 sections have unique maximal integral
overlattices on which every displayed facet reflection remains integral.  The
three resulting even 2-elementary carrier types are classical:

\[
(10,4,0)=U\oplus D_4(-1)^2,\qquad
(10,6,0)=U(2)\oplus D_4(-1)^2,\qquad
(10,8,0)=U\oplus E_8(-2).
\]

After global sign reversal or scaled duality, every carrier is one of two
strongly squarefree rank-ten classes in the Scharlau--Kirschmer census.  Their
full reflection chambers are the classical Vinberg--Kaplinskaja `P9` chamber
or Tumarkin's `9-tum04_11_01` chamber.  All five rational quadratic spaces are
similar to the odd unimodular form (I_{9,1}).  Consequently these seeds do
not define a new arithmetic commensurability class or a new maximal lattice
reflection group.

There is one exact literal diagram identification in the pinned CoxIter
corpus: `A14_8_0_C01` is Im Hof's doubly truncated orthoscheme
`9-imhof85_03`.  The other four diagrams are unmatched in that bounded
corpus; this is not a worldwide priority theorem.

## Exact carrier table

| seed | facets | maximal index | determinant | maximal carrier | normalization | SK class and full chamber | literal pinned-corpus match | disjoint doubling walls |
|---|---:|---:|---:|---|---|---|---|---|
| `A14_8_0_C01` | 12 | 4 | -64 | (U(2)+D_4(-1)^2), ((10,6,0)) | (-2M^{-1}) | List9 entry 23, Tumarkin H9/11 | Im Hof `9-imhof85_03` | none |
| `A14_8_1_C02` | 16 | 8 | -16 | (U+D_4(-1)^2), ((10,4,0)) | (-M) | List9 entry 23, Tumarkin H9/11 | none | `(15,16)` |
| `A14_10_1_C01` | 13 | 4 | -256 | (U+E_8(-2)), ((10,8,0)) | (-2M^{-1}) | List9 entry 7, P9 H9/10 | none | none |
| `A14_10_1_C02` | 15 | 8 | -64 | (U(2)+D_4(-1)^2), ((10,6,0)) | (-2M^{-1}) | List9 entry 23, Tumarkin H9/11 | none | `(6,15)` |
| `A16_14_1_C01` | 15 | 32 | -16 | (U+D_4(-1)^2), ((10,4,0)) | (-M) | List9 entry 23, Tumarkin H9/11 | none | `(1,14)`, `(8,15)` |

The exact orders of the reflection-stable integral subgroups in the five
discriminant modules are, respectively,

```text
A14_8_0_C01   1, 4
A14_8_1_C02   1, 2, 4, 8
A14_10_1_C01  1, 2, 4
A14_10_1_C02  1, 2, 4, 8
A16_14_1_C01  1, 16, 32
```

In each row there is one maximal subgroup, so the maximal stable overclosure
is unique under inclusion.  “Maximal” here has a precise local scope: all
integral overlattices of the frozen fixed lattice are enumerated through its
discriminant bilinear module, and each candidate must be preserved by every
displayed facet reflection.  It does not enumerate unrelated rescalings of
individual wall-normal lines.

## Identification with full reflection groups

The certificate replays the two relevant Scharlau--Kirschmer entries rather
than relying on determinant alone.

- Entry 7 has genus symbol `II_(9,1)(2^2)`, determinant `-4`, discriminant
  quadratic-value histogram `{0:3, 1:1}`, and full chamber `P9` with 10
  facets, 8 finite vertices, and 2 cusps.

- Entry 23 has genus symbol `II_(9,1)(2^4)`, determinant `-16`, discriminant
  quadratic-value histogram `{0:10, 1:6}`, and full chamber
  `9-tum04_11_01` with 11 facets, 24 finite vertices, and 2 cusps.

For P9, the ten archived root normals span a determinant-(-4) lattice with
the entry-7 discriminant form.  For the Tumarkin chamber, the eleven root
normals have rank ten and primitive null relation

```text
(-1,-2,-3,-2,-1,0,1,2,3,2,1);
```

nodes 2 through 11 span the determinant-(-16) entry-23 lattice.  Signature,
Smith form, the exact discriminant quadratic form, and Nikulin uniqueness for
indefinite even 2-elementary lattices identify the integral classes.  An
explicit (mathrm{GL}_{10}(\mathbb Z)) congruence matrix is not claimed.

The displayed section reflection groups have finite index in these known full
lattice reflection groups because both chambers have finite hyperbolic
volume.  Their numerical indices are not recovered from Coxeter Euler
characteristics in odd hyperbolic dimension: the relevant Euler
characteristics vanish.

Exact determinant square classes and local Hasse invariants at every place
also show that all five forms are rationally similar to
(-x_0^2+x_1^2+\cdots+x_9^2).  This is the precise basis for the “no new
arithmetic commensurability class” conclusion.

## Literal diagram lineage

The bounded comparison scanned 1,214 ordinary and 2,959 index-two/doubling
`.coxiter` files at CoxIter commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe`.  It compares complete
edge-labelled graphs, including parallel and ultraparallel relations.

`A14_8_0_C01` has three relation-histogram collisions, Im Hof records `_02`,
`_03`, and `_04`, but only `_03` is graph-isomorphic.  The certified
candidate-to-reference wall permutation is

```text
(6,11,10,9,8,1,2,3,4,5,7,12).
```

The other bounded negative results are:

- `A14_8_1_C02` (H9/16): no same-header record;
- `A14_10_1_C01` (H9/13): eight same-header records, no histogram collision;
- each H9/15 seed: eight same-header records, no histogram collision.

This cleanly separates an exact bounded-corpus statement from a literature
priority claim.  The conservative status for those four literal diagrams is
“worldwide priority open.”

The named-stock exclusions are also exact but narrow:

- Allcock's Figure 4.1 polytope (Q\subset H^9) has 14 facets, so none of
  the five is literally (Q) by facet count.
- The current ADEG definition requires mutually intersecting facets and only
  labels 2, 3, and 6.  Every seed here contains label 4, so none is an ADEG
  polyhedron; none has the 14 facets of the exceptional H9 polyhedron (P^*).
- Tumarkin's complete (n+2)-facet classification concerns H9/11, so none of
  the five is literally one of those chambers.  His compact (n+3) result
  does not apply because all five seeds are noncompact.
- The Vinberg--Kaplinskaja P9 chamber has 10 facets, so no seed is literally
  P9, even though the scaled dual of `A14_10_1_C01` has its lattice class.

Thus “the arithmetic envelope is classical” and “the displayed chamber is a
known literal diagram” are different statements.

## Redoubling consequences

Under Allcock's definition, parallel walls are disjoint in hyperbolic space
even though their closures meet at infinity.  The exact parallel doubling
pairs in the table therefore make three seeds redoublable:

```text
A14_8_1_C02   (15,16)
A14_10_1_C02  (6,15)
A16_14_1_C01  (1,14), (8,15)
```

Allcock's Lemma 1.3 gives iterated finite-volume doubles, with volume
multiplied by (2^j) after (j) doublings.  Corollary 4.2 then gives
reflection-generated subgroups of every positive finite index for each of
these three chamber groups.

No seed has three certified pairwise disjoint doubling walls.  Therefore the
exponential subgroup-count conclusion of Allcock's Theorem 4.3 does not
follow from the displayed walls.  In particular, the three pairwise
orthogonal doubling walls of `A14_10_1_C01` meet rather than being disjoint,
so they do not make that seed redoublable in Allcock's sense.

## Conditional algebro-geometric interpretation

The carrier names have a classical K3 meaning, but it is conditional.  If one
uses a carrier as the full invariant lattice
(S=H^2(X,\mathbb Z)^+) of a nonsymplectic involution on a K3 surface, then
Nikulin's classification and the Alexeev--Engel account give:

| invariant type | lattice | moduli dimension | fixed locus of the involution |
|---|---|---:|---|
| ((10,4,0)) | (U+D_4(-1)^2) | 10 | one genus-4 curve and 3 disjoint rational curves |
| ((10,6,0)) | (U(2)+D_4(-1)^2) | 10 | one genus-3 curve and 2 disjoint rational curves |
| ((10,8,0)) | (U+E_8(-2)) | 10 | two disjoint smooth elliptic curves |

The dimension is (20-r).  Outside the exceptional ((10,8,0)) case, the
fixed locus is (C_g\sqcup k\mathbb P^1), where

\[
g=11-\frac{r+a}{2},\qquad k=\frac{r-a}{2}.
\]

For ((10,8,0)), the source gives a rational elliptic quotient whose
ramification divisor consists of two smooth fibers, hence a concrete double
cover construction.  These are known ten-dimensional K3-involution moduli
types; no new deformation type and no explicit polynomial normal form are
claimed here.

Most importantly, the arithmetic overclosure computation by itself does not
promote the original five fixed sections to invariant lattices of actual K3
involutions.  The table describes what the same abstract maximal lattices mean
if such a K3 promotion is separately supplied.

## Novelty boundary

Alexeev's published input is the ambient mixed-root reflection data.  The
project contribution is the exact fixed-section/divisor-reflection extraction
and its finite-volume certification.  The safe wording is:

> new exact Alexeev fixed-section/divisor-reflection extractions; one abstract
> chamber is Im Hof's and all five arithmetic classes have classical full
> reflection overgroups.

There is no claim of a new commensurability class, maximal reflection group,
doubling mechanism, K3 involution deformation type, or worldwide diagram
priority.

## Frozen artifacts and replay

- source five-seed archive:
  `data/enumerations/alexeev_parked_rank10_h9_examples.json`, SHA-256
  `e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430`;
- certificate:
  `data/enumerations/alexeev_parked_rank10_h9_lineage_overclosures.json`;
- builder:
  `scripts/certify_alexeev_parked_rank10_h9_lineage_overclosures.py`;
- replay test:
  `tests/test_alexeev_parked_rank10_h9_lineage_overclosures.py`.

```bash
PYTHONPATH=.:src python3 scripts/certify_alexeev_parked_rank10_h9_lineage_overclosures.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_alexeev_parked_rank10_h9_lineage_overclosures.py
```

The external corpora and primary-source snapshots are hash-pinned inside the
certificate, including CoxIter, Scharlau--Kirschmer `List9/Gram9`, Allcock's
PDF/source/Figure 4, and the Alexeev--Engel source tarball.

## Primary references

- V. Alexeev, [*Reflective hyperbolic 2-elementary lattices*](https://arxiv.org/abs/2209.09110).
- V. Alexeev and P. Engel, [*Compactifications of moduli spaces of K3 surfaces with a nonsymplectic involution*](https://arxiv.org/abs/2208.10383).
- V. Alexeev and V. Nikulin, [*Classification of log del Pezzo surfaces of index at most 2*](https://arxiv.org/abs/math/0406536).
- V. Nikulin, [*Integral symmetric bilinear forms and some of their geometric applications*](https://www.mathnet.ru/eng/im1677).
- D. Allcock, [*Infinitely many hyperbolic Coxeter groups through dimension 19*](https://doi.org/10.2140/gt.2006.10.737).
- P. Tumarkin, [*Hyperbolic Coxeter n-polytopes with n+2 facets*](https://arxiv.org/abs/math/0301133).
- N. Bredon, [*On ADEG-polyhedra in hyperbolic spaces*](https://arxiv.org/abs/2507.05153).
- R. Scharlau and M. Kirschmer, [reflective Lorentzian lattice census](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html).
- H.-C. Im Hof, *A class of hyperbolic Coxeter groups*, Expositiones
  Mathematicae 3 (1985), 179--186; and *Napier cycles and hyperbolic Coxeter
  groups*, Bull. Soc. Math. Belgique A 42 (1990), 523--545.
