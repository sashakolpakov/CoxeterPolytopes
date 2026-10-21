# Broader literature lineage audit for `D0004/C007` in `H10`

## Bottom line

The lineage boundary is sharper than the bounded CoxIter audit alone showed.

The following objects were already known:

1. the rank-11 even 2-elementary lattice
   `M = U + D6(-1) + A1(-1)^3`, with `(r,a,delta)=(11,5,1)`;
2. its finite-volume 14-wall chamber for the full integral root-reflection
   group;
3. its 22-wall `(-2)`-reflection, or K3 nef, chamber; and
4. concrete K3 surfaces having `NS(X) = M`, including an explicit double-plane
   model.

Moreover, finite-volume Coxeter groups in `H10` were plentiful long before
this computation.  The candidate is therefore **not** a new dimension, stable
lattice, maximal arithmetic reflection group, or first algebro-geometric
realization of that lattice.

What has not been located is the literal object certified here: the mixed-root
27-wall chamber cut out of Shimada chamber `D0004`, its reflection group of
index `2016` in `W(M)`, or that particular divisor-fixed extraction.  That is a
bounded priority question, not a worldwide novelty theorem.

## Exact prior identification of the 14-wall envelope

Alexeev--Nikulin's [*Classification of log del Pezzo surfaces of index at most
2*](https://arxiv.org/abs/math/0406536), later published as *Del Pezzo and K3
surfaces*, gives the decisive earlier identification.

- Table 1, entry `N=32`, is exactly `(r,a,delta)=(11,5,1)`.
- Theorem 3.1.1 says that the displayed diagram is the fundamental chamber of
  the full reflection group `W^(2,4)(M)=W(M)` and has finite volume.
- In the appendix, the entry is constructed as
  `M = U + D6(-1) + 3 A1(-1)`.  For `m=6,l=3`, its displayed root list has
  `7` alpha roots, `4` beta roots, the root `e`, and `2` gamma roots: exactly
  `14` walls.
- Proposition 3.4.1 proves that every elliptic triplet, hence `(11,5,1)`, is
  realized by a standard super-extremal K3 pair whose diagram is this full
  reflection chamber.

This is the same maximal chamber independently identified in the project with
Scharlau--Kirschmer `List10/Gram10` entry 20.  The maintained
[Scharlau--Kirschmer table](https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html)
records fundamental domains of **full** lattice reflection groups.  It does
not claim to enumerate their finite-index reflection-subgroup chambers.

Thus the exact 14-wall envelope was public no later than the 2004 arXiv
version and 2006 memoir.  Alexeev--Nikulin describe that work as an expansion
of earlier joint work, but the 2004/2006 source is the conservative checkpoint
used here.

## A concrete surface with this Neron--Severi lattice

Roulleau's [*An atlas of K3 surfaces with finite automorphism group*](https://arxiv.org/abs/2003.08985)
gives a particularly concrete geometric incarnation.  Its rank-11 section
treats

`NS(X) = U + D4(-1)^2 + A1(-1)`.

This lattice has signature `(1,10)`, 2-elementary invariants `(11,5,1)`, and
determinant `32`.  By Nikulin's uniqueness theorem for indefinite even
2-elementary lattices, it is isometric to

`U + D6(-1) + A1(-1)^3 = M`.

Roulleau gives all `22` smooth rational `(-2)`-curves, says that eleven of
them generate `NS(X)`, and exhibits both an ample divisor of square `48` and a
nef divisor of square `2`.  The latter realizes `X` as a double cover of
`P^2` branched over

`L + L' + Q4`,

where `L,L'` are lines and the residual quartic passes through their
intersection, producing a `D4` singularity of the sextic.  The paper also
describes the images of the 22 rational curves.

This is a direct answer to the algebro-geometric design question for the
**maximal lattice**: a concrete family of surfaces already realizes it.
However, the relevant K3 nef chamber has 22 `(-2)`-curve walls.  It is not the
27-wall chamber, whose primitive wall types in `M` are

`19*(-2) + 8*(-4)`.

Alexeev--Nikulin make the group-theoretic distinction explicit: the K3 nef
chamber is the chamber of `W^(2)(M)`, whereas their 14-wall chamber is the
chamber of `W^(2,4)(M)`.  The four black roots in entry `N=32` form
`A2+A2`, so the nef chamber is a union of `|W(A2+A2)|=36` maximal 14-wall
chambers.  Roulleau's 22 curves are consistent with this larger nef chamber.
The candidate's index `2016` and its eight primitive `(-4)` reflections show
that it is neither of these two standard chambers.

There is also an exact scoped absence here: Roulleau's complete rank-11 table
has `12,14,22,90` rational curves in its four finite-automorphism cases, not
27.  This excludes the candidate as a rank-11 **nef chamber in that atlas**;
it does not exclude a mixed-root finite-index subgroup chamber.

## The Shimada stock was known; the fixed section was derived here

Shimada's [*An algorithm to compute automorphism groups of K3 surfaces and an
application to singular K3 surfaces*](https://arxiv.org/abs/1304.7427)
records the ambient source of the construction.  Table 1 lists the singular
K3 surface with transcendental Gram data `[2,0,8]`, discriminant `16`, and
`4539` equivalence classes of induced chambers.  The accompanying
`disc16.txt` data is the source of chamber `D0004`.

Shimada's paper and supplement therefore establish the prior existence of the
ambient divisor lattice and the `D0004` chamber.  They do not, in the located
material, state the order-four `C007` fixed-section construction, its 27
strict facets, the stable overlattice `M`, or the index `2016`.  This is a
source inspection result, not a theorem that no equivalent construction was
published elsewhere.

## `H10` itself was already densely populated

Allcock's [*Infinitely many hyperbolic Coxeter groups through dimension
19*](https://arxiv.org/abs/0903.0138) proves that there are infinitely many
finite-covolume Coxeter groups in every dimension below 20.  In the
construction section he specifically notes redoublable finite-volume
polyhedra in `H10` coming from Vinberg's chambers for both `I_(10,1)` and the
form `-2 x0^2+x1^2+...+x10^2`.

Consequently:

- `H10` was filled by classical arithmetic examples;
- infinite families in `H10` were known by 2009; and
- a 27-facet count, even if new as a literal diagram, would not constitute a
  dimension record or a new existence theorem for `H10`.

The standard few-facet classifications do not decide this candidate.  For
example, Tumarkin's [finite-volume `n+2` classification](https://arxiv.org/abs/math/0301133)
and [compact `n+3` classification](https://arxiv.org/abs/math/0406226) have
scope `12` and `13` facets in dimension ten, whereas this noncompact chamber
has `27=n+17` facets.  The maintained
[Felikson survey](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html)
also presents a map of partial classifications, not a complete census of all
finite-volume `H10` subgroup chambers.

## What the searches do and do not prove

The current rigorous negative statements are only these:

- no 27-facet `H10` item occurs in the pinned 1,214-file ordinary CoxIter
  corpus;
- no such item occurs in the pinned 2,959-file reflection-doubling corpus;
- the candidate is not the 14-wall full chamber in Alexeev--Nikulin or
  Scharlau--Kirschmer;
- it is not the 22-wall K3 nef chamber in Roulleau's rank-11 atlas; and
- exact web searches for `H10`, `27 facets`, the lattice/genus, and index
  `2016` did not locate a primary source defining the candidate.

The last item is only a keyword non-location.  Diagrams can be published as
images, under different lattice bases, or as unexpanded doubling/gluing
constructions, so it is not an absence theorem.

The strongest safe priority statement is therefore:

> `D0004/C007` is a proved positive algebraic extraction and an unresolved
> priority candidate for the literal 27-wall, index-2016 reflection-subgroup
> chamber and its particular divisor-fixed construction.  Its dimension,
> stable lattice, full reflection group, maximal chamber, and a concrete K3
> realization of that lattice were known previously.

Worldwide lineage would close positively by finding a published diagram or
root list with an exact normalized-Gram permutation/isometry, or a published
index-2016 subgroup construction that reproduces the 27 walls.  A global
negative result would require a classification of finite-index reflection
subgroups of this `W(M)` at the relevant index; the surveyed maximal-lattice
and few-facet tables do not provide one.

## Frozen source checkpoints

The primary arXiv source archives were downloaded on `2026-08-06`.  Their
SHA-256 hashes were:

```text
math/0406536  8fa1248ecc6aeb7bf8317a12a151a662d88cde5d2d8fb978eb21bd6145161fbf
2003.08985    3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9
1304.7427     4a6406a092d245c81f2fee099cbad4220a926a86e5ed12907b23ba8457814c63
0903.0138     dbf49889ab9a2529aafeff69ee627c7ea4920c9f94d56f093e5db18a52055113
```

The exact candidate fingerprint and all computational proofs remain in the
[main lineage report](shimada_disc16_d4_c7_h10_27_lineage.md) and its linked
certificate.
