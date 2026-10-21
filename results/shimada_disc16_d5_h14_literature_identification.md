# Literature identification of the Shimada `disc16` `D5` fold in `H14`

## Verdict

The 60-facet finite-volume noncompact chamber obtained from the
discriminant-16 singular-K3 divisor data is **not a new hyperbolic Coxeter
polytope**.  It is exactly Valery Alexeev's published reflective
two-elementary case

\[
(r,a,\delta)=(15,7,1).
\]

This is an exact identification, not a match by dimension, facet count, or
commensurability.  After a permutation of the 60 facets, all 3,600 entries of
the integral primitive-root Gram matrix agree with the matrix distributed in
the source of

- V. Alexeev,
  [*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler
  varieties*](https://arxiv.org/abs/2209.09110), arXiv:2209.09110v4;
- published as Chapter 9, pp. 130--150, in *Higher Dimensional Algebraic
  Geometry* (Cambridge University Press, print year 2025),
  [DOI 10.1017/9781009396233.010](https://doi.org/10.1017/9781009396233.010).

The result is still useful: the fishing mechanism independently reconstructs
a difficult 60-wall reflective chamber from concrete Néron--Severi chamber
data.  What it does not supply is a new polytope.

## Exact intrinsic lattice

Starting only with Shimada's matrix and chamber data, the saturated fixed
lattice is

\[
S=\operatorname{Fix}(M)
 \cong U(2)\oplus E_8(-1)\oplus A_1(-1)^5.
\]

The exact invariants are

| invariant | value |
|---|---:|
| rank | 15 |
| signature | `(1,14)` |
| determinant | 128 |
| Smith diagonal | `1^8, 2^7` |
| discriminant group | `(Z/2Z)^7` |
| two-elementary triple | `(15,7,1)` |

The integral-kernel Smith calculation proves saturation.  Nikulin's
uniqueness theorem for indefinite even two-elementary lattices then identifies
this with Alexeev's lattice, not merely with a rational quadratic space.  See
V. V. Nikulin,
[*Integral symmetric bilinear forms and some of their geometric
applications*](https://doi.org/10.1070/IM1980v014n01ABEH001060), especially
Theorem 3.6.2.

The 60 primitive mirrors split as

| primitive square and divisibility | count |
|---|---:|
| `(-2, div 1)` | 15 |
| `(-2, div 2)` | 35 |
| `(-4, div 2)` | 10 |

All 60 reflections are integral involutive isometries of `S`.  The 15
`(-2, div 1)` roots form an integral basis: their coordinate determinant is
`-1`, their Gram determinant is 128, and their graph is the edge-subdivision
`K5^(2)`.

This reproduces Alexeev's description exactly.  His case `(15,7,1)` has 15
main roots forming a basis with diagram `K5^(2)`, followed by 10 `(-4)` roots,
15 additional `(-2)` roots, and 20 further `(-2)` roots.  His Theorem 2 also
gives

\[
\operatorname{Aut}\Gamma_r=S_5,
\qquad O^+(S)=S_5\ltimes W_r.
\]

Consequently the reflectivity of this exact integral lattice was already a
theorem in Alexeev's work.

## Exact chamber match

Alexeev's arXiv source includes the full matrix in
`all_cases_in_2el_paper.txt`.  The pinned files used here have hashes

```text
arXiv source tar:
10e4002f6518f8bf6d669946775bbe42731f69277a52847c5602fa60ea33c2a6

all_cases_in_2el_paper.txt:
4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457
```

The comparison treats each matrix as a complete edge-coloured graph, with the
root square as the vertex colour and every inner product, including zero, as
an edge colour.  A colour-preserving isomorphism is found and then checked
entry by entry.  Before relabelling, both matrices already have the same exact
histograms:

| entries | value: count |
|---|---|
| diagonal | `-2: 50`, `-4: 10` |
| off diagonal | `0: 865`, `1: 20`, `2: 395`, `4: 120`, `6: 310`, `10: 60` |

The archived permutation is therefore a compact, replayable certificate that
the two weighted Coxeter diagrams—and hence the hyperbolic chambers—are
isometric.

## Compactness and volume status

The chamber is finite-volume and **noncompact**.

The independent rational double-description found 28,872 projective extreme
rays:

| ray type | count |
|---|---:|
| timelike (finite vertices) | 27,845 |
| lightlike (ideal vertices) | 1,027 |
| spacelike (hyperideal escapes) | 0 |

Thus there is no hyperideal escape and the chamber has finite volume; the
1,027 ideal rays prove noncompactness.  Independently, Alexeev reports exactly
1,027 maximal parabolic subdiagrams, in 20 orbits under `S5`.  The equality of
these counts is a strong cross-check, although the full Gram-matrix
isomorphism is the decisive identification.

## What was and was not already in the sources

Shimada's primary input is

- I. Shimada,
  [*An algorithm to compute automorphism groups of K3 surfaces and an
  application to singular K3 surfaces*](https://arxiv.org/abs/1304.7427);
- the official
  [`threediscs.zip`](https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AlgoAutK3/threediscs.zip)
  supplement, whose relevant `disc16.txt` member has SHA-256
  `a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6`.

Section 10.5 of the paper records the discriminant-16 surface, 4,539 chamber
classes, maximum level 17, 3,308 generators, and the curve-orbit bound.  The
supplement records the `D5` chamber, its 165 walls, and its stabilizer.  In the
bounded audit performed here, neither source states the saturated rank-15
fixed lattice, the irredundant 60-wall restriction, or its finite-volume
certificate.

Accordingly, the narrow Shimada-to-Alexeev reconstruction may be an unrecorded
realization.  That is not yet a priority theorem, and Alexeev already gives an
algebro-geometric K3 interpretation of the `(15,7,1)` lattice.  No surface
automorphism hypothesis is required for the chamber identification: the
construction uses the divisor lattice and its chamber symmetry as a factory.

## Position in the `H14` landscape

As of 6 August 2026, none of the following novelty claims survives:

| possible claim | status |
|---|---|
| new weighted Coxeter diagram | **No:** exact Alexeev match |
| new reflective integral lattice | **No:** Alexeev proves `(15,7,1)` reflective |
| new reflection group | **No:** same chamber, hence conjugate reflection group |
| new commensurability class | **No:** equality is stronger than commensurability |
| first finite-volume example in `H14` | **No** |
| compact example in `H14` | **No:** this chamber has 1,027 cusps |
| new Shimada-divisor derivation | possible, but not established |

There were already several finite-volume noncompact populations in `H14`:

1. Vinberg's arithmetic unit-form chambers.  The current CoxIter corpus has a
   17-wall chamber for `-x0^2+x1^2+...+x14^2` and a 20-wall chamber for
   `-2x0^2+x1^2+...+x14^2`, both noncocompact.  The primary lineage is E. B.
   Vinberg,
   [*On groups of unit elements of certain quadratic forms*](https://doi.org/10.1070/SM1972v016n01ABEH001346),
   *Math. USSR-Sb.* 16 (1972), 17--35.
2. Vinberg's nonarithmetic construction gives a cofinite noncocompact
   reflection group in `H14` as a `D4` Coxeter face of his `H18` mixture.  See
   E. B. Vinberg,
   [*Non-arithmetic hyperbolic reflection groups in higher
   dimensions*](https://doi.org/10.17323/1609-4514-2015-15-3-593-602),
   *Moscow Math. J.* 15 (2015), 593--602, Section 8.
3. Alexeev's exact 60-wall `(15,7,1)` chamber has been on arXiv since 2022
   and is now published by Cambridge University Press.

The general compact problem remains very far from this example.  The current
Felikson--Tumarkin
[hyperbolic Coxeter polytope status page](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html)
records compact examples only through dimension 8, whereas finite-volume
examples are known through dimension 19 and in dimension 21.  This is a
known-examples status statement, not a nonexistence theorem for compact
`H14` polytopes.

## Why the CoxIter-only audit missed it

The pinned official CoxIter snapshot at commit
[`b800be48`](https://github.com/rgugliel/CoxIter/commit/b800be48d8240fbabd2d1c9542e61054584ad7fe)
contains 1,214 `.coxiter` files.  Its three `H14` records have respectively
17, 17, and 20 facets, and the largest graph anywhere in that snapshot has 50
facets.  Therefore the exact 60-facet chamber is absent from that corpus.

That bounded absence was correct but was never evidence of global novelty:
Alexeev distributes large matrices with his paper rather than through the
CoxIter graph collection.  This is precisely why the fishing pipeline must
perform an invariant-guided literature comparison, not just a corpus-name or
facet-count lookup.

## Reproducible certificates

- [fixed-lattice certificate](../data/enumerations/shimada_disc16_d5_h14_fixed_lattice_certificate.json)
- [exact Alexeev match](../data/enumerations/shimada_disc16_d5_h14_alexeev_exact_match.json)
- [independent finite-volume audit](shimada_disc16_d5_h14_independent_audit.md)
- [bounded CoxIter corpus audit](../data/enumerations/shimada_disc16_d5_h14_coxiter_corpus_audit.json)
- [matrix-matching script](../scripts/match_shimada_disc16_d5_h14_alexeev_15_7_1.py)
- [fixed-lattice script](../scripts/certify_shimada_disc16_d5_h14_fixed_lattice.py)

With the pinned Alexeev source attachment unpacked at the archived path, the
new identification replays as

```sh
PYTHONPATH=src:. python3 scripts/match_shimada_disc16_d5_h14_alexeev_15_7_1.py --check-archive
python3 -m pytest -q tests/test_shimada_disc16_d5_h14_alexeev_exact_match.py
```
