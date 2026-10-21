# The toric `H6/9` endpoint is a 216-fold congruence refinement of `P^6`

**Status (2026-08-09): exact independent audit, no cdd.** This adds a
group-theoretic and smoothing interpretation to the smooth-toric endpoint; it
does not produce a new abstract Coxeter polytope.

Let \(L=I_{1,6}\), let \(W_{\rm tor}\) be the reflection group of the toric
nine-cycle `(4,4,3)^3`, and let \(\Gamma_2^6\) be the principal level-two
reflection group whose chamber is the Potyagailo--Vinberg polytope \(P^6\).
Then

\[
 W_{\rm tor}=\Gamma_2^6\rtimes W(A_2)^3,
 \qquad [W_{\rm tor}:\Gamma_2^6]=6^3=216.                 \tag{1}
\]

Equivalently, \(W_{\rm tor}\) is the full inverse image of its order-216
mod-two image, and one \(P^6\) chamber is tiled by exactly 216 toric endpoint
chambers.

## The missing kernel argument

The six square-minus-two boundary roots split into three orthogonal \(A_2\)
pairs. Their integral reflections generate a finite subgroup

\[
 F=W(A_2)^3,
\]

of order 216. Direct enumeration of its integer \(7\times7\) matrices finds
216 distinct reductions modulo two, so \(F\cap\ker(\rho_2)=1\).

The three remaining boundary roots are \(E_4,E_5,E_6\), of square minus one.
Their \(F\)-orbits are three disjoint sets of nine roots. Their union is
coefficientwise exactly the classical stock

\[
 E_i,\qquad H-E_i-E_j,\qquad
 2H-\sum_{j\ne i}E_j,                                      \tag{2}
\]

of all 27 exceptional classes on a smooth cubic surface. Reflections in the
classes (2) generate \(\Gamma_2^6\): this is the congruence-two theorem used in
the exact \(P^6\) identification of Everitt--Ratcliffe--Tschantz. Therefore
all 27 kernel generators already lie in \(W_{\rm tor}\). Conversely, the
three square-minus-one toric generators lie in \(\Gamma_2^6\), while the other
six lie in \(F\). This proves (1), rather than inferring it from a numerical
volume ratio.

There is also an exact chamber-containment check. Each class in (2) is a
`0/1` sum of the nine toric boundary classes. Thus every toric-nef class is
nonnegative on all 27 cubic-line classes, so the chosen toric chamber is
inside the chosen \(P^6\) chamber. The subgroup index then gives the claimed
216-chamber tiling with no orientation ambiguity.

## Independent Euler and volume check

Solomon's formula, summed over all 322 spherical standard parabolics of the
nine-cycle, gives

\[
 \chi(W_{\rm tor})=-\frac1{1728}.
\]

For the right-angled \(P^6\) group, the orthogonal-clique vector is

```text
1, 27, 216, 720, 1080, 648, 72,
```

and hence \(\chi(\Gamma_2^6)=-1/8\). Their ratio is 216. Six-dimensional
Gauss--Bonnet, \(\operatorname{vol}/\pi^3=-(8/15)\chi\), gives

\[
 \operatorname{vol}(P^6)=\frac{\pi^3}{15},\qquad
 \operatorname{vol}(P_{\rm tor})=\frac{\pi^3}{3240}.       \tag{3}
\]

## What the `3A2` smoothing statement means

The anticanonical contraction of the toric surface is

\[
 S_0:\ XYZ=W^3\subset\mathbf P^3.
\]

In characteristic different from three its singular locus consists of the
three coordinate points. In the chart at each point the equation is
\(uv=w^3\), so the singularities are \(3A_2\). The six toric
square-minus-two curves are exactly the three \(A_2\) exceptional chains.
The only lines on \(S_0\) are

\[
 W=X=0,\qquad W=Y=0,\qquad W=Z=0.
\]

Indeed, on a line not contained in \(W=0\), restriction would factor \(s^3\)
as the product of the three restricted linear forms \(X,Y,Z\); unique
factorization forces all four coordinates to be proportional, a
contradiction. Inside \(W=0\), the equation \(XYZ=0\) leaves exactly the
three displayed lines.

An explicit characteristic-zero smoothing is

\[
 XYZ-W^3-t(X^3+Y^3+Z^3)=0.                                 \tag{4}
\]

The total space is smooth at the three singular points because the
\(t\)-derivative is nonzero there. For \(t\ne0\), the gradient equations can
have a projective solution only if \((3t)^3=1\); hence a generic fibre of (4)
is a smooth cubic. On the marked simultaneous-resolution cover, the local
Weyl deck group is \(W(A_2)^3\), and the 27 nearby exceptional classes split
into the three size-nine orbits computed above. This is the precise sense in
which the three central toric lines unfold to the 27-line/\(P^6\) wall stock.

One caveat is essential: the full product \(W(A_2)^3\) belongs to the product
local simultaneous-resolution marking cover. A single one-parameter pencil
such as (4) need not have the full product as its monodromy image, and
individual nearby lines need not extend without base change. The assertion
here is a marked Picard-lattice/chamber statement, not an unqualified claim
about the monodromy of one pencil. The simultaneous-resolution input in this
form is the theorem of J. Wahl,
[*Simultaneous resolution of rational singularities*](https://www.numdam.org/article/CM_1979__38_1_43_0.pdf),
whose Galois group is the direct product of the Weyl groups attached to the
square-minus-two configurations.

## Novelty boundary and replay

The labelled endpoint is Im Hof's known dimension-six Napier-cycle polytope.
The abstract polytope is therefore not new. Project-new here are the smooth
toric full-nef realization, its exact congruence refinement of \(P^6\), and
the marked `3A2 -> 27 lines` chamber interpretation. No worldwide priority
is claimed for those interpretations.

Primary lineage inputs are H.-C. Im Hof,
[*Napier cycles and hyperbolic Coxeter groups*](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/poly/napier.pdf)
(1990), and Everitt--Ratcliffe--Tschantz,
[*Right-angled Coxeter polytopes, hyperbolic 6-manifolds, and a problem of
Siegel*](https://arxiv.org/abs/1009.3231), especially their congruence-two
description of \(P^6\).

The independent certificate uses only tiny integer matrices and at most
`2^9` Coxeter subsets:

```bash
PYTHONPATH=.:src python3 scripts/certify_toric_h6_congruence_smoothing.py
PYTHONPATH=.:src pytest -q tests/test_toric_h6_congruence_smoothing.py
```
