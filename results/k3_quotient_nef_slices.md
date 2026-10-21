# Automorphism quotients of K3 surfaces give exact Coxeter divisor cones

## Result

The algebro-geometric construction does not have to stop at the full nef cone
of a K3 surface.  If a finite group (G\subset\operatorname{Aut}(X)) acts on
a projective K3 surface, then the fixed nef section is itself the divisor cone
of the normal projective quotient surface

\[
Y=X/G.
\]

Applied to two of Kondō's high-rank finite-automorphism K3 surfaces, this
gives exact finite-volume Coxeter divisor cones in dimensions 9 and 10:

| (NS(X)) | group element | (\operatorname{rk}NS(X)^G) | quotient-nef polytope | vertices | cusps |
|---|---|---:|---|---:|---:|
| (U\oplus E_8(-1)^2) | symplectic involution | 10 | Coxeter simplex in (\mathbb H^9), 10 facets | 10 | 2 |
| (U\oplus E_8(-1)^2\oplus A_1(-1)) | a transposition in (\mathfrak S_3) | 11 | Coxeter polytope in (\mathbb H^{10}), 13 facets | 38 | 4 |

All remaining vertices are finite: respectively 8 and 34.  Thus both
polytopes have finite hyperbolic volume and are noncompact.  This conclusion
uses the literal all-subset audit of Vinberg's criterion given below; it is
not inferred from the displayed vertex and cusp counts alone.  The polytopes
themselves are **not claimed to be new**; the result here is their exact
realization as divisor cones of algebraic quotient surfaces and the
reproducible folding certificate.

Kondō proves that the graph actions used here are induced by actual
symplectic automorphisms: the unique nontrivial symmetry in his case
(U\oplus E_8\oplus E_8), and the full (\mathfrak S_3) symmetry in his
case (U\oplus E_8\oplus E_8\oplus A_1).  See
[Kondō, Sections 4.1--4.2](https://doi.org/10.1017/S0027763000001653).
The complete curve configurations are also printed in
[Roulleau, Sections 15.4--15.5](https://epiga.episciences.org/10377/pdf).

## Why this is genuinely a divisor-cone construction

Let (X) be a normal projective variety, let (G) be finite, and put
(\pi:X\to Y=X/G).  Pullback gives

\[
\pi^*:N^1(Y)_{\mathbb R}\xrightarrow{\sim}N^1(X)_{\mathbb R}^{G}.
\]

Indeed, the norm identity

\[
\pi^*\pi_*D=\sum_{g\in G}gD
\]

gives descent over (\mathbb R).  The projection formula gives

\[
D\text{ nef on }Y
\iff \pi^*D\text{ nef on }X.
\]

Consequently

\[
\boxed{
\pi^*\operatorname{Nef}(Y)
=\operatorname{Nef}(X)\cap N^1(X)_{\mathbb R}^{G}.}
\]

Intersection forms differ only by the degree:

\[
(\pi^*D)^2=|G|D^2.
\]

After normalization, this constant rescaling does not alter the hyperbolic
metric.  The fixed section is therefore not an artificial linear slice: it
is the normalized nef cone of the algebraic surface (Y).  Equivalently it
is the polarization cone of the quotient stack ([X/G]).

The two covering involutions are symplectic, so (Y) has eight (A_1)
quotient singularities and trivial canonical class in the orbifold sense.
Moreover (X) is a Mori-dream K3 surface, and a surjective image of a Mori
dream space is again a Mori dream space by
[Okawa](https://arxiv.org/abs/1104.1326).  Thus finite generation of the
quotient Cox ring and polyhedrality of its divisor cones are built into the
construction.  Polyhedrality of this **complete quotient-nef cone** already
implies finite hyperbolic
volume in relative dimension at least two by the convex--Hodge argument.  One
may apply Hodge index on the normal quotient, or equivalently pull back to the
complete invariant-nef section of `X`.  The exact vertex audit independently
certifies the extracted wall cone and its finite/ideal census; it is not being
used to compensate for an arbitrary incomplete divisor subsystem.

## Folding the curve walls

Let (C) be a ((-2))-curve wall.  On the invariant subspace, all
inequalities in one (G)-orbit coincide.  The corresponding normal is the
orbit sum

\[
R_{\mathcal O}=\sum_{C\in\mathcal O} C.
\]

For the two involutions used here every contributing orbit is either:

1. a fixed curve, with (R_{\mathcal O}^2=-2); or
2. two disjoint curves, with (R_{\mathcal O}^2=-4).

In the second case the ambient root reflections commute and

\[
(s_Cs_{gC})|_{NS(X)^G}=s_{C+gC}.
\]

Hence every folded wall reflection preserves the invariant divisor lattice.
The mixture of norm (-2) and norm (-4) roots naturally creates Coxeter
edges of label 4; this is the usual crystallographic folding mechanism, now
performed on a complete algebraic curve configuration.

If an orbit sum has positive square and pairs positively with an ample
class, it lies in the chosen timelike cone.  Its orthogonal hyperplane misses
hyperbolic space and its nef inequality is automatic.  This occurs exactly
once in the rank-19 case: the orbit ({A_{23},A_{24}}) has

\[
(A_{23}+A_{24})^2=-4+2\cdot6=8
\]

and is therefore redundant.

## The dimension-9 quotient-nef simplex

For (NS(X)=U\oplus E_8(-1)^2), Kondō constructs the involution by swapping
the two (E_8) factors and fixing (U).  Thus

\[
NS(X)^G\simeq U\oplus E_8(-2),
\]

of rank 10 and signature ((1,9)).  In Roulleau's wall labels, the ten
orbits are

\[
\begin{gathered}
18|1, 17|2, 15|3, 16|4, 14|5,
13|6, 12|7, 11|8, 9, 19|10.
\end{gathered}
\]

Writing these orbit normals as (R_1,\ldots,R_{10}), their Gram matrix is

\[
\begin{pmatrix}
-4&2&0&0&0&0&0&0&0&0\\
2&-4&2&0&0&0&0&0&0&0\\
0&2&-4&2&2&0&0&0&0&0\\
0&0&2&-4&0&0&0&0&0&0\\
0&0&2&0&-4&2&0&0&0&0\\
0&0&0&0&2&-4&2&0&0&0\\
0&0&0&0&0&2&-4&2&0&0\\
0&0&0&0&0&0&2&-4&0&2\\
0&0&0&0&0&0&0&0&-2&2\\
0&0&0&0&0&0&0&2&2&-4
\end{pmatrix}.
\]

It has determinant (-1024) and inertia ((1,9)).  The label-3 edges are

\[
12,23,34,35,56,67,78,8,10,
\]

and (9,10) has label 4; all other pairs are orthogonal.  This is a
finite-volume hyperbolic Coxeter 9-simplex.  Exact inversion gives ten
projective vertices: eight timelike and two isotropic.  At each isotropic
vertex the active Gram matrix is negative semidefinite of rank 8 and
corank 1, giving the two affine cusp links.

## The dimension-10 quotient-nef polytope

For (NS(X)=U\oplus E_8(-1)^2\oplus A_1(-1)), choose a transposition in
Kondō's symplectic (\mathfrak S_3).  Its fixed divisor lattice has rank 11.
After deleting the one redundant timelike orbit, thirteen wall orbits remain.
Their Gram matrix has rank 11 and inertia ((1,10,2)).  Number the facets in
the orbit order stored in the certificate.  The nonorthogonal Coxeter edges
are:

- label 3:
  (1,11,23,34,45,56,67,78,89,7,12);
- label 4: (12,9,10);
- parallel: (10,13,11,13).

There are no ultraparallel facet pairs.  Exhausting all rank-10 active wall
sets gives 38 projective extreme rays: 34 timelike and four isotropic.  The
four ideal links are negative semidefinite of rank 9.  They agree with the
four cusps in the finite-covolume certificate below.

## Literal Vinberg--CoxIter finite-covolume audit

The applicable statements are Guglielmetti's cocompactness and
finite-covolume criteria,
[CoxIter, p. 761, Propositions 1.9--1.10](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/A117FD150E108CACCACF9EF00C61A23C/S1461157015000273a.pdf/div-class-title-coxiter-computing-invariants-of-hyperbolic-coxeter-groups-div.pdf#page=8).
For the Coxeter diagram of a finite-sided chamber in \(\mathbb H^n\), it
requires:

1. at least one spherical subdiagram of rank \(n\), or one parabolic
   subdiagram of rank \(n-1\); and
2. every spherical subdiagram of rank \(n-1\) to have exactly two extensions,
   each either spherical of rank \(n\) or parabolic of rank \(n-1\).

The input hypotheses are checked before these two bullets are applied.  The
published curve lists are complete; the actual surface involutions preserve
their Gram matrices; the full orbit fold is used; every positive-square orbit
normal is removed only after its mirror is shown to miss hyperbolic space;
and every remaining inequality is proved to be an irredundant facet.  After
changing from the surface sign convention, the folded matrices are
indecomposable acute Coxeter Gram matrices with inertia
\((n,1,m-n-1)\), and every normalized pair is an allowed Coxeter angle or a
parallel pair.  These facts, rather than the subdiagram counts in isolation,
put the two matrices under Proposition 1.10.

The verifier then performs a second audit independent of its optimized
zero-neighborhood routine: it scans every spherical rank-\((n-1)\) and
rank-\(n\) principal subdiagram and **all** componentwise-parabolic
rank-\((n-1)\) principal subdiagrams, not merely a supplied list of maximal
parabolics, and counts containment extensions literally.  The result is:

| chamber | spherical rank \(n-1\) | spherical rank \(n\) | parabolic rank \(n-1\) | endpoint-type histogram |
|---|---:|---:|---:|---|
| \(\mathbb H^9\) | 45 | 8 | 2 | \((2,0):28,\ (1,1):16,\ (0,2):1\) |
| \(\mathbb H^{10}\) | 208 | 34 | 4 | \((2,0):137,\ (1,1):66,\ (0,2):5\) |

Here \((f,i)\) means that an edge has \(f\) finite and \(i\) ideal
extensions.  Thus all 45 and all 208 edges have exactly two endpoints.  In
the \(\mathbb H^{10}\) case, seven connected affine subdiagrams collapse to
four unique full parabolic neighborhoods; three of the four are reducible
with two affine components.  Every connected affine diagram has one of those
full rank-nine neighborhoods, and no finite component or other orthogonal
wall contaminates a putative cusp.

This edge-extension criterion is itself an if-and-only-if theorem, so it does
not require a separate Lannér test.  Lannér cones occur as the independent
second condition in the alternative full-neighborhood formulation of
[Heckman's Theorem 5.12](https://www.math.ru.nl/~heckman/CoxeterGroups.pdf#page=106);
checking every edge extension is an equivalent route, not a deletion of that
condition.  The direct extreme-ray enumeration and the full-neighborhood
audit agree exactly with the literal Proposition 1.10 endpoint sets and are
retained as independent geometric cross-checks.

## A sharp rank obstruction for nonrational eigenspaces

The quotient construction uses the rational (+1)-eigenspace of an
involution.  It cannot replace the real-quadratic construction needed for a
compact chamber.  The distinction has an exact rank formulation.

Let (V_{\mathbb Q}) have signature ((1,\rho-1)), and let a totally real
field (K) of degree (e>1) act self-adjointly on a nondegenerate
(K)-submodule of (V).  If its distinguished real eigenspace has signature
((1,d)), then every conjugate eigenspace has the same real dimension
(d+1).  Distinct eigenspaces are orthogonal.  Because (V) has only one
positive direction, all conjugate eigenspaces except the distinguished one
must be negative definite.  Therefore

\[
\boxed{\rho\ge e(d+1).}
\]

Equivalently, for a self-adjoint semisimple correspondence (\tau), an
irreducible degree-(e) factor of its minimal polynomial with a Lorentzian
eigenspace of dimension (d+1) consumes at least (e(d+1)) divisor ranks.

For a complex K3 surface, (\rho\le20).  Hence a nonrational eigenspace in
hyperbolic dimension at least 9 forces

\[
e=2,qquad d=9,qquad \rho=20,
\]

with no negative-definite complement.  In particular:

- none of the rank-15 through rank-19 finite-automorphism K3 lattices can
  carry a real-quadratic Lorentzian eigenspace of dimension at least 9;
- a K3 real-multiplication construction cannot reach (\mathbb H^{10});
- the rank-20 golden construction is not merely convenient but forced.

Rational involution-fixed spaces evade the doubling bound and produce the
two quotient examples above.  However they again live over an integral
rational Lorentzian lattice.  If their reflection chamber has finite volume,
Meyer's theorem makes the arithmetic group nonuniform in dimension at least
4, so this rational quotient route cannot produce a compact high-dimensional
chamber.  The two and four explicit cusps are the geometric manifestation.

Thus the two algebro-geometric routes now have cleanly separated jobs:

1. **finite quotients of Mori-dream K3 surfaces** automatically give finite
   algebraic wall data and can produce high-dimensional finite-volume
   Coxeter polytopes;
2. **totally-real correspondence eigenspaces** are necessary for the compact
   dimension-9 target, and on a K3 they force Picard rank 20.

## Reproduction

The source action data are in
[`k3_quotient_nef_slices.json`](../data/known_examples/k3_quotient_nef_slices.json).
The verifier consumes that file and the independently transcribed complete
K3 wall matrices:

```sh
PYTHONPATH=src:. python3 scripts/verify_k3_quotient_nef_slices.py
PYTHONPATH=src:. pytest -q tests/test_k3_quotient_nef_slices.py
```

It checks the wall-graph automorphisms, invariant ranks and signatures,
integrality of every folded reflection, all Coxeter angles, facet
nonredundancy, every projective extreme ray, and every parabolic cusp link
using exact rational arithmetic.  It also performs and freezes the literal
all-subset Proposition 1.10 audit and checks that its endpoint sets agree
with the independently optimized full-neighborhood calculation.
