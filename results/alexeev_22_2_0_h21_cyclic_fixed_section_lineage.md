# Lineage audit for the `(22,2,0)` 210-wall chamber and its cyclic sections

Last updated: 2026-08-06.

## Bottom line

The ambient `H21` chamber is classical, not an Alexeev discovery and not a
candidate new example.  Borcherds gave the exact determinant-four lattice,
the `42+168=210` walls, and the finite-index reflection statement in his 1985
thesis and in journal form in 1987.  Dolgachev--Kondō gave the characteristic-
two K3 realization in 2003.  Alexeev's new contribution for this row is the
complete maximal-parabolic census: `1,095,990` diagrams in `18` symmetry
orbits.

The cyclic fixed sections are different constructions, but the first positive
one is **not a new abstract chamber**.  The unique order-three class `3b` with
fixed rank `12` cuts out a 22-facet finite-volume noncompact Coxeter polytope
in `H11`.  Its saturated fixed lattice is the class-number-one lattice
`U+A2(-1)^5`, and its 22 walls are the complete reflection chamber recorded
as Scharlau--Kirschmer `List11/Gram11` entry 27.  Thus its reflection-subgroup
index is one.  What remains potentially unrecorded is the particular
characteristic-two K3 fixed-section realization, not the lattice, chamber, or
abstract hyperbolic polytope.  The `H10` cast has one exact positive volume
audit, a now-identified classical lattice, and exact subgroup index
`1,566,720`, while `H13` and `H14` remain unscreened.

There is one especially important geometric qualification.  The full
weighted-diagram group has order `241,920`, whereas Dolgachev--Kondō's finite
subgroup realized by K3 automorphisms has order `120,960`.  Thus a cyclic
section obtained from the former is automatically a divisor-lattice section,
but it is a fixed construction on the surface itself only when its conjugacy
class meets the geometric subgroup.

## Exact identity of the ambient object

Let

`L = U + D20(-1)`.

Equivalently, `L` is the even sublattice of the odd unimodular Lorentzian
lattice `I_(21,1)`.  It has rank `22`, signature `(1,21)`, determinant `4`,
and 2-elementary invariants `(r,a,delta)=(22,2,0)`.

The full mixed-root reflection chamber has

- `42` short walls from primitive roots of square `-2`;
- `168` long walls from primitive reflecting roots of square `-4`; and
- weighted-diagram symmetry group
  `PSL(3,F4).D12`, of order `241,920`, where `D12` means the dihedral
  group of order `12`.

The archived Alexeev matrix has exact rank `22` and inertia
`(1,21,188)` in the convention `(positive,negative,zero)`.  This is a useful
independent replay of the Lorentzian realization, not the source of its
priority.

## Publication timeline

### Borcherds, 1985 thesis

Borcherds' [*The Leech lattice and other lattices*](https://math.berkeley.edu/~reb/papers/thesis/thesis.pdf),
Chapter 6, Example 5, takes `R=D4` in `II_(25,1)` and identifies its relevant
orthogonal lattice with the even sublattice of `I_(21,1)`.  It states that the
fundamental domain has `168` norm-four walls and `42` norm-two walls, that its
diagram automorphism group is `L3(4).D12` of order
`2^8*3^3*5*7=241,920`, and that the reflection group has finite index in the
lattice automorphism group.  The thesis introduction already summarizes this
as a `168+42`-sided fundamental domain.  The public PDF is a corrected 1999
copy of the Cambridge 1985 thesis; the result itself belongs to the original
thesis.

### Borcherds, 1987 journal paper

Example 5 of [*Automorphism groups of Lorentzian lattices*](https://math.berkeley.edu/~reb/papers/al/al.pdf),
J. Algebra **111** (1987), 133--153, pp. 148--149, repeats the exact
construction and explicitly calls the 210-wall domain a fundamental domain
for the reflection group.  It also records a useful negative contrast: `L`
lies in three copies of `I_(21,1)`, but the reflection groups of those odd
overlattices are not of finite index.  The cofinite object is the even
determinant-four lattice, not any of those three odd reflection groups.

Finite index in the arithmetic orthogonal group proves finite hyperbolic
covolume.  The chamber is noncompact: Alexeev's later enumeration exhibits
maximal parabolic subdiagrams, hence cusps.

### Esselmann, 1996

Esselmann's rank theorem makes this ambient example extremal on the theorem
side: up to the qualifications in his classification, it is the unique
reflective even hyperbolic lattice at rank `22`, and there is no even
reflective hyperbolic lattice of elliptic type in rank `21` or rank at least
`23`.  This concerns the ambient reflective lattice.  It does not prohibit
lower-rank fixed sublattices or finite-index reflection-subgroup chambers in
dimensions `10`--`14`.

### Dolgachev--Kondō, 2003

Dolgachev--Kondō's
[*A supersingular K3 surface in characteristic 2 and the Leech lattice*](https://arxiv.org/abs/math/0112283)
identifies `L` with the Picard lattice of the unique supersingular K3 surface
of Artin invariant one in characteristic two.  Among its equivalent concrete
models are the quasi-elliptic Weierstrass model, an inseparable double plane,
a quartic model, and a complete intersection in `P2 x P2`.

Their Sections 3 and 5 provide the exact geometric dictionary:

- the `42` short roots are classes of actual smooth rational `(-2)`-curves,
  in two families of `21`, with the point--line incidence graph of
  `P2(F4)`;
- the `168` long reflections are realized by `168` explicit involutive
  Cremona transformations; and
- the 210-wall domain is the intersection of a Conway chamber in
  `II_(1,25)` with the K3 positive cone.

Thus this is a genuinely algebro-geometric source of the mixed reflection
chamber.  It is nevertheless not correct to call all 210 walls divisor walls
of the K3 nef cone.  The 42 short mirrors come from effective rational curves;
the 168 long mirrors are automorphisms acting as `(-4)` reflections.  The
mixed chamber is a Picard-lattice reflection chamber inside the ample-cone
story, rather than the K3 nef chamber itself.

### Alexeev, 2022/2023

In the `(22,2,0)` subsection of
[*Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler varieties*](https://arxiv.org/abs/2209.09110),
Alexeev explicitly attributes the 42-plus-168 diagram to Borcherds and
Dolgachev--Kondō.  His new exact computation is the maximal-parabolic
census: `1,095,990` maximal parabolic subdiagrams, `18` modulo the diagram
group, all in the odd genus.  Every listed affine-root type has total finite
root rank `20`, the cusp rank appropriate to the rank-22 Lorentzian lattice.

Accordingly the safe attribution is:

| datum | earliest source pinned here |
|---|---|
| determinant-four rank-22 lattice and 210-wall cofinite chamber | Borcherds thesis, 1985 |
| journal publication of the chamber | Borcherds, 1987 |
| maximal-rank uniqueness theorem | Esselmann, 1996 |
| concrete characteristic-two K3 and geometric meaning of both root orbits | Dolgachev--Kondō, 2003 |
| complete maximal-parabolic census | Alexeev, 2022/2023 |

## Full diagram symmetry is not the geometric finite subgroup

Dolgachev--Kondō identify

`Aut(D(X)) = M21.D12 = PSL(3,F4).D12`,

of order `241,920`.  They then determine which part is realized on the
surface.  The dihedral factor is generated by a field-Frobenius involution, a
switch interchanging the point and line families, and a projective element of
order three.  They prove that the Frobenius isometry is **not** induced by an
automorphism of the K3 surface.  The geometric finite subgroup is

`PGL(3,F4).2`,

of order `120,960`, and every finite subgroup of `Aut(X)` is conjugate into
it.  The full automorphism group is the extension by the normal subgroup
generated by the 168 Cremona involutions.

This gives a mandatory label for every cyclic fishing result:

1. **diagram-lattice section:** the generator lies in the full 241,920-element
   symmetry group, so `L^g` and the literal chamber section are defined;
2. **geometric section:** the conjugacy class meets `PGL(3,F4).2`, so after
   conjugacy the fixed lattice comes from a finite K3 automorphism; or
3. **non-geometric diagram section:** the class lies outside every conjugate
   of that subgroup.  It remains an algebro-geometric *lattice/chamber*
   extraction from `NS(X)`, but not a fixed locus of a surface automorphism.

Merely observing that the generator permutes the 210 walls does not decide
between (2) and (3).

## Exact prime-order class dictionary

An independent exact reconstruction starts from the point--line Levi graph
of `P2(F4)` on its `21+21=42` vertices.  Its graph automorphism group has
order `241,920`.  The Gram matrix `A-2I` of the 42 short roots has rank `22`.
In the GAP character table `L3(4).D12`, the degree-42 permutation character
and its Neron--Severi quotient are

`chi_42 = X.1 + X.5 + X.7 + X.11`,

`chi_NS = X.1 + X.5 + X.11`.

The resulting prime-order class dictionary is exact.  Here “fixed short
walls” counts fixed vertices in the 42-vertex Levi action, and hyperbolic
dimension is fixed Neron--Severi rank minus one.

| GAP/ATLAS class | class size | fixed short walls | `tr(g|NS)` | fixed rank | section | meets `PGL(3,F4).2`? |
|---|---:|---:|---:|---:|---:|---|
| `2a` | 315 | 10 | 6 | 14 | `H13` | yes |
| `2b` | 280 | 0 | -2 | 10 | `H9` | no |
| `2c` | 360 | 14 | 8 | 15 | `H14` | no: field Frobenius |
| `2d` | 1,008 | 0 | 0 | 11 | `H10` | yes: point--line switch class |
| `3a` | 2,240 | 6 | 4 | 10 | `H9` | yes |
| `3b` | 672 | 12 | 7 | 12 | `H11` | yes |
| `3c` | 1,920 | 0 | 1 | 8 | `H7` | yes |

Thus the high-dimensional fishing queue has a precise geometric split:
`2d/H10`, `3b/H11`, and `2a/H13` are geometric, while `2c/H14` is the
non-geometric Frobenius diagram-lattice section.  In particular, `3b` is the
unique full-group order-three class with fixed rank `12`.  Its geometricity
does not depend on a delicate fusion choice: every element of odd order lies
in every index-two subgroup, hence in Dolgachev--Kondo's index-two geometric
subgroup.

## Positive section: `SS2-3B-H11-22` (exact lineage)

Take a generator in class `3b`.  Its 210 ambient inequalities form 78 wall
orbits on the rank-12 fixed space.  Exact restriction and global
irredundancy leave

- 22 negative genuine facet normals;
- 50 isotropic inactive restrictions; and
- 6 timelike inactive restrictions.

Every one of the `binom(22,2)=231` genuine-facet pairs is Coxeter.  Their
complete angle census is

`m_2=186,  m_3=35,  m_6=10`,

with no parallel or ultraparallel facet pair.

Two independent exact Vinberg enumerators agree on the entire edge and
vertex census:

| datum | exact count |
|---|---:|
| elliptic rank-10 edge diagrams | 14,580 |
| finite vertices | 2,304 |
| ideal vertices / full cusp neighborhoods | 72 |
| finite--finite edges | 10,965 |
| finite--ideal edges | 3,414 |
| ideal--ideal edges | 201 |

The three endpoint counts sum to all 14,580 edges, and every edge has exactly
two endpoints.  A further 180 exact affine-neighborhood audits find zero bad
neighborhoods.  These are the literal Vinberg checks that distinguish finite
volume from the weaker facts “Lorentzian Gram matrix” or “some parabolic
subdiagrams exist.”  They prove that the 22-facet section is finite-volume
and noncompact.

The saturated fixed lattice has rank 12, signature `(1,11)`, determinant
`-243=-3^5`, and Smith diagonal `1^7 3^5`.  The Scharlau--Kirschmer census
lists exactly one class in the corresponding genus,
`II_(11,1)(3^-5)`: one-based `List11/Gram11` entry 27.  Its displayed Gram is
equivalently `U(3)+2A2+E6` in the opposite sign convention; the unique
integral class can also be written `U+A2(-1)^5`.  The census entry records a
full chamber with 22 walls, 72 cusps, and elliptic-subdiagram vector

`[22,231,1480,6320,18704,38986,57182,57785,38140,14580,2304]`.

Direct enumeration of every elliptic subdiagram of ranks 1 through 11 for
the section gives that vector entry for entry, independently of the two
endpoint enumerators above.

The section walls are primitive integral roots of exactly two types:

- 12 singleton-orbit roots of type `(-2,div 1)`; and
- 10 triple-orbit roots of type `(-6,div 3)`.

All 22 reflections preserve the saturated lattice.  Their normalized facet
Gram has SHA-256
`2944c4ce6ebab9a5fb1b52c9a0f7c2a105656863036c3e2badbe52b15c9763cb`.

There is also a direct all-root Vinberg stopping proof, rather than merely a
match of numerical invariants.  If `p` is one third of the sum of the 22
primitive walls, then `p` is integral and primitive, `p^2=14`, and
`div(p)=1`.  The 12 short walls have `(p,r)=1`, while the 10 long walls have
`(p,r)=3`.  The exponent-three discriminant form permits only reflective-root
types `(-2,div 1)` and `(-6,div 3)`.  Exact positive-definite enumeration in
`p^perp` finds no reflective root at height zero; the two first affine shells
contain exactly the 12 and 10 displayed oriented roots and no additional
root.  Since the resulting chamber is already finite-volume, Vinberg's
stopping criterion proves completeness.  Hence the section group is the full
reflection group of the census lattice and the index is **one**.
The replayable [production certificate](../data/enumerations/borcherds_char2_3b_h11_22_certificate.json)
and independent [all-root lineage certificate](../data/enumerations/borcherds_char2_3b_h11_22_lineage.json)
freeze the two assertions separately.

Consequently the hyperbolic construction is positive but its abstract
lineage is closed: `SS2-3B-H11-22` is precisely the known Scharlau--Kirschmer
entry-27 chamber.  The bounded Borcherds--Dolgachev--Kondo--Alexeev audit did
not locate this *K3 fixed-section derivation*, so only that algebro-geometric
realization remains priority-open.

## Positive section: `SS2-2D-H10-47`

The geometric point--line switch class `2d` has fixed rank 11 and produces 47
genuine facets in `H10`.  One exact graph-based Vinberg enumeration currently
gives

- 19,910 elliptic rank-9 edge diagrams;
- 3,136 finite vertices and 216 ideal vertices;
- endpoint types 12,050 finite--finite, 7,260 finite--ideal, and 600
  ideal--ideal; and
- 568 passing affine-neighborhood audits, with zero bad edges or audits.

An independent exact Schur-complement replay returns the same edge, endpoint,
cusp, and affine-neighborhood counts.  The two implementations therefore
certify finite-volume noncompactness at theorem grade.

The lattice lineage is already much sharper.  The saturated fixed lattice
has rank 11, signature `(1,10)`, determinant `2^12`, and Smith diagonal
`2^10 4`.  Its form is divisible by two; after global rescaling by `1/2` it
is the odd determinant-two lattice in Scharlau--Kirschmer `List10/Gram10`
entry 2, equivalently the classical `f_2^10` class.  The full reflection
chamber there has 13 walls and 3 cusps.  Thus the 47-wall object is a
finite-index subgroup chamber in an old arithmetic lattice, not a new lattice
or first example in `H10`.

A complete Solomon sum over all finite special subgroups of the 47-wall
Coxeter system gives `chi(W_47)=-527/120960`.  The full 13-wall `f_2^10`
group has `chi(W_full)=-31/11147673600`; hence the exact index is

`chi(W_47)/chi(W_full)=1,566,720=2^11*3^2*5*17`.

The finite-special counts by rank are
`[1,47,540,3230,12215,31212,54582,64050,47760,19910,3136]`, independently
recovering the 19,910 elliptic edges and 3,136 finite vertices in the first
volume audit.  The complete face vector, including the 216 ideal vertices, is
`(3352,19910,47760,64050,54582,31212,12215,3230,540,47,1)`, and
Gauss--Bonnet gives `vol(P_47)=pi^5*527/3572100`.  The priority of the K3
fixed-section realization remains open.

The exact [production certificate](../data/enumerations/borcherds_char2_2d_h10_47_certificate.json)
freezes the 47-facet Gram with SHA-256
`454e34ef2358603578a4360145ac6bdd7abfe6e0b539344323ffa5e354038881`
and replays the lattice, irredundancy, integrality, and Vinberg tests from the
pinned 210-wall source.  A separate [lineage/index certificate](../data/enumerations/borcherds_char2_2d_h10_47_lineage.json)
gives an explicit determinant-minus-one integral basis change
`C^T (G_fixed/2) C = diag(2,-1^10)` and independently recomputes both
Solomon sums.  Thus the `f_2^10` identification is an exact isometry, not an
inference from determinant and Smith data alone.

The pinned CoxIter ordinary corpus has no `H10` entry above 14 facets, and its
separate reflection-doubling archive has no `H10` entry above 44 facets.
Consequently neither bounded corpus contains a literal 47-facet candidate.
This narrows the priority question for the subgroup chamber but is not a
worldwide non-location theorem.

## What is and is not known about the cyclic sections

The fixed-space mechanism itself is not new.  Borcherds' 1998 normalizer
paper gives a general Coxeter-theoretic construction of reflection groups on
fixed spaces of finite parabolic subgroups, and arithmeticity gives finite
covolume for stabilizers of suitable rational fixed subspaces.  Neither fact
identifies the literal polytope `D(X) intersect H^g`, its irredundant walls, or
whether its wall reflections generate the full stabilizer.

The following source inspection has so far found no explicit statement of the
*K3 fixed-section derivations* for the cyclic `H10`, `H11`, `H13`, or `H14`
casts.  The abstract `3b/H11` chamber itself has now been identified exactly
in the Scharlau--Kirschmer census:

- Borcherds' 1985 thesis and 1987 journal article;
- Borcherds' 1998 normalizer paper;
- Dolgachev--Kondō 2003;
- Alexeev's article and attached 210-by-210 matrix; and
- exact web searches combining `42`, `168`, `210`, `PSL(3,F4)`, `U+D20`,
  “fixed lattice”, “fixed subspace”, “cyclic”, and the target dimensions.

This establishes only a bounded non-location of the algebraic derivation.
It is no longer evidence for abstract novelty of the `H11` chamber.  The
`H10` literal subgroup chamber still requires exact comparison after its
second volume audit and index computation.

## Section-level admission gates

For each cyclic conjugacy class retained by the fishing computation, archive
the following before adding a positive row to the main ledger:

1. the exact generator or stable conjugacy-class identifier in
   `Aut(G_210)` and its order;
2. the saturated lattice `L_g = L intersect Fix_Q(g)`, its Gram matrix,
   signature, determinant, Smith invariants, parity, and discriminant form;
3. the complete restricted inequality stock and an exact irredundancy proof;
4. the normalized squared Gram matrix of the genuine facets, its canonical
   hash, and a strict Coxeter-angle check for every facet pair;
5. a literal Vinberg finite-volume certificate or an exact cone enumeration
   proving that every extremal ray is timelike or lightlike and none is
   hyperideal;
6. the finite/ideal endpoint and maximal-parabolic census, so “finite volume”
   is not inferred merely from the existence of some cusps;
7. the reflection-stable lattice when it differs from the divisor fixed
   lattice, together with any containing maximal reflection chamber and
   subgroup index;
8. exact normalized-Gram comparison with the pinned CoxIter ordinary and
   doubling corpora; and
9. membership, up to conjugacy, in the geometric subgroup
   `PGL(3,F4).2`.

Only gates 1--6 establish a new literal hyperbolic Coxeter polytope.  Gates
7--9 establish its arithmetic and algebro-geometric lineage.

## Frozen source checkpoints

The primary source files inspected on 2026-08-06 had these SHA-256 hashes:

```text
Borcherds thesis PDF          d0b3ee3d8b1e2c17e1736b7bfae7809cbea196202c8ab845c76aab92d291572e
Borcherds 1987 PDF            78b89beb3cffc077d38766644a59dd0d52ae54781cb2e9eb9d9d4b76ca0a851f
Dolgachev--Kondo arXiv source 1a5d39a3332308c0325f201e8f615ee84c409eb8d57b8cf7639bade72931ffe6
Alexeev attached matrix       4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457
Scharlau--Kirschmer Lorentz.zip 74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5
Scharlau--Kirschmer List11    2ca0725232cb101e60380f14ad0ab7777478c60eff06991ba591b6d3c71e3184
Scharlau--Kirschmer Gram11    af8d52984fc9caabc00ffe670ae646d2151e68d141082544e18f4827ba142540
Scharlau--Kirschmer List10    9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21
Scharlau--Kirschmer Gram10    b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18
```

The files inside the Scharlau--Kirschmer archive carry timestamps of
2022-07-29, and the public archive has an HTTP last-modified date of
2023-11-20.  Thus both census matches indisputably predate this project.
For `H10` the lineage is much older still: Vinberg's 1972 groups-of-units
paper explicitly treats the forms `-2*x0^2+x1^2+...+xn^2` for `n<14` and
computes their reflection polyhedra.

The absence statement above is frozen to these source inspections and the
named searches.  It is not an assertion that no equivalent section exists
anywhere in the literature.
