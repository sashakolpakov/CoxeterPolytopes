# TOOT: an equation-level K3 carrier for `SS2-2A H13/50`

## Positive theorem

There is an explicit six-parameter Jacobian elliptic K3 family whose very
general member has

\[
 \operatorname{NS}(X)=U\oplus D_8(-1)\oplus D_4(-1)
\]

and carries all fifty divisor mirrors of the certified finite-volume
noncompact `SS2-2A H13/50` Coxeter chamber.

For parameters \(a,b,c,d,e,f\), set

\[
\begin{aligned}
P(t)&=-3t^4+6t^3+at^2+bt+c,\\
Q(t)&=2t^6-6t^5+(3-a)t^4+(1+a-b)t^3+dt^2+et+f.
\end{aligned}
\]

Take the minimal resolution of

\[
 \boxed{\quad y^2=x^3+t^2P(t)x+t^3Q(t).\quad}
\]

On a nonempty Zariski-open parameter set its singular fibers are

\[
 I_4^*+I_0^*+8I_1.
\]

The \(I_4^*\) and \(I_0^*\) components supply \(D_8(-1)\) and \(D_4(-1)\),
while the fiber and zero section supply \(U\).  The certificate proves that
this trivial lattice is saturated and is the entire very-general Picard
lattice.  A determinant-minus-one 14-by-14 marking then attaches the full
50-mirror chamber, reflection by reflection.

This is an equation-level realization of an existing chamber and a published
rank-14 K3 family.  It is not a new abstract Coxeter polytope and no worldwide
priority is claimed.

## Fiber calculation

At \(t=0\),

\[
 \operatorname{ord}(a_4,a_6,\Delta)=(2,3,6),
 \qquad
 \Delta/t^6\big|_{t=0}=-16(4c^3+27f^2).
\]

Away from \(cf(4c^3+27f^2)=0\), the
[Kodaira fiber criterion](https://doi.org/10.2307/1970131) gives an exact
\(I_0^*\) fiber.

At infinity put \(z=1/t\), \(X=z^4x\), and \(Y=z^6y\).  The equation becomes

\[
Y^2=X^3+z^2\bar P(z)X+z^3\bar Q(z),
\]

where

\[
\begin{aligned}
\bar P&=-3+6z+az^2+bz^3+cz^4,\\
\bar Q&=2-6z+(3-a)z^2+(1+a-b)z^3+dz^4+ez^5+fz^6.
\end{aligned}
\]

Direct expansion gives

\[
4\bar P^3+27\bar Q^2
=-9\bigl(a^2+6a+12b-12c-12d+9\bigr)z^4+O(z^5).
\]

Thus the orders at infinity are \((2,3,10)\), giving exactly \(I_4^*\),
off the displayed coefficient's zero locus.

Globally, the remaining discriminant factor

\[
D(t)=4P(t)^3+27Q(t)^2
\]

has degree eight.  If it is squarefree and coprime to \(P\), its eight roots
are \(I_1\) fibers.  The Euler sum is

\[
6+10+8=24.
\]

The open set is genuinely nonempty.  At
\((a,b,c,d,e,f)=(-2,-2,-2,-2,-2,-2)\),

\[
\begin{aligned}
D(t)={}&-225t^8+702t^7-833t^6+480t^5-300t^4\\
       &+172t^3+132t^2+120t+76.
\end{aligned}
\]

Its polynomial discriminant is

```text
-308186920225388788221712962876472295424
```

and \(\operatorname{Res}(D,P)=11489984737344\).  Both are nonzero; also
\(D(0)=76\) and the leading coefficient is \(-225\).  This single rational
witness verifies every required open condition.

The degrees of \(a_4=t^2P\) and \(a_6=t^3Q\) are 6 and 9, so after
homogenization they are sections of \(\mathcal O_{\mathbf P^1}(8)\) and
\(\mathcal O_{\mathbf P^1}(12)\).  The minimal elliptic surface therefore has
fundamental line bundle \(\mathcal O(2)\); its smooth minimal resolution is a
K3 surface.

## Why the very-general Picard lattice is exact

The six parameters are recovered directly from the normalized coefficients:

\[
(a,b,c,d,e,f)=([t^2]P,[t]P,P(0),[t^2]Q,[t]Q,Q(0)).
\]

The two uniquely typed fibers fix \(0\) and \(\infty\) on the base.  Under
\(t\mapsto\lambda t\) and a short-Weierstrass homothety
\((x,y)\mapsto(\mu^2x,\mu^3y)\), preserving the leading pair \((-3,2)\)
forces \(\mu^2=\lambda^3\).  Preserving the normalized \(6t^5\) coefficient
then forces \(\lambda=1\), leaving only the finite ambiguity \(\mu^2=1\).
The marked-to-unmarked map therefore has zero-dimensional fibers, so this is
a genuine six-dimensional K3 family.

The fiber divisors already give Picard rank at least 14.  The period space of
rank-14 lattice-polarized K3 surfaces has dimension six, while the locus of
Picard rank at least 15 has dimension at most five; this is also the standard
period-domain count in
[Dolgachev's lattice-polarized K3 construction](https://doi.org/10.1007/BF02362332).
Hence the very general member has \(\rho=14\).  Shioda--Tate now gives

\[
14=2+8+4+\operatorname{rank}MW,
\]

so the Mordell--Weil rank is zero.

Rank does not by itself exclude a finite-index saturation.  For a nonzero
section, [Shioda's height formula](https://doi.org/10.14992/00009994) is

\[
h(R)=4+2(R\cdot O)-\sum_v\operatorname{contr}_v(R).
\]

At \(I_4^*/D_8\), the multiplicity-one nonidentity components have correction
values \(1,2,2\); at \(I_0^*/D_4\), they have values \(1,1,1\).  Therefore

\[
h(R)\ge 4-(2+1)=1.
\]

A nonzero torsion section would have height zero, which is impossible.
Consequently the Mordell--Weil torsion is trivial, the trivial lattice is
saturated, and

\[
\operatorname{NS}(X)=U+D_8(-1)+D_4(-1)
\]

for a very general member.

This conclusion has an independent primary-source denominator.
[Clingher--Malmendier, *On K3 surfaces of Picard rank 14*](https://doi.org/10.1002/mana.202200197),
Lemma 2.2 and Theorem 3.7, classify the standard Jacobian fibration on a
general
\(P'_{14}=U+D_8(-1)+D_4(-1)\)-polarized K3 surface as precisely
\(I_4^*+I_0^*+8I_1\) with trivial Mordell--Weil group.  Their Theorem 1.5
gives the same six-dimensional coarse moduli space.  The present equation is
a normalized short-Weierstrass chart of that published standard fibration,
not a newly claimed K3 family.  The exact audited version is
[arXiv:2009.09635v6](https://arxiv.org/abs/2009.09635).

## Exact 14-by-14 bridge to the chamber

Let \(G_R\) be the Gram matrix in Roulleau's archived curve basis
\(A_1,\ldots,A_{14}\), and let \(S\) be the standard block Gram of
\(U+D_8(-1)+D_4(-1)\).  The certificate records an integral matrix \(T\) with

\[
\det T=-1,
\qquad
T G_R T^{\mathsf T}=S.
\]

The first standard basis vector is the affine \(D_8\) fiber

\[
F=2A_1+2A_2+2A_3+2A_4+A_5+2A_6+A_7+A_{16}+A_{17},
\]

where \(A_{16},A_{17}\) are expressed in the archived 14-element basis by
Roulleau's exact wall-coordinate table.  The second standard vector is
\(F+A_8\), with \(A_8\) a section, and these two vectors have Gram
\(\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)\).  The remaining
twelve rows split as the displayed standard \(D_8(-1)+D_4(-1)\) root basis.
All fourteen rows of \(T\) are frozen in the script and JSON certificate.

Composing this marking with the prior scaled-dual map gives a second explicit
unimodular matrix \(B\), of determinant 1, satisfying

\[
B^{\mathsf T}S B=\operatorname{Gram}(M^\vee(2)),
\]

where \(M\) is the original characteristic-two `SS2-2A` fixed lattice.
This is an actual coordinate bridge, not merely an isometry-class comparison.

Transporting the fifty normals through \(B\) gives

```text
40 primitive (-2), divisibility 1
10 primitive (-4), divisibility 2.
```

Their 14-by-50 coordinate matrix has Hermite normal form \(I_{14}\), so the
mirrors generate the entire Picard lattice.  Every one of the fifty
reflections is integral, unimodular, involutive, and Gram preserving in the
literal fiber-divisor basis.  The transported strict interior has square
2012 and minimum wall pairing 6.

## Literal Vinberg conclusion

The volume statement is inherited only after the exact lattice and halfspace
transport.  The prior literal Vinberg and independent exact-cone audits give:

| invariant | value |
|---|---:|
| hyperbolic dimension | `13` |
| facets | `50` |
| finite vertices | `45,954` |
| cusps | `2,598` |
| elliptic rank-12 edges | `391,614` |
| connected affine components | `6,210` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |
| hyperideal rays in the independent exact cone | `0` |

Thus the chamber is finite-volume and noncompact.  This conclusion does not
follow merely from the K3 equation or its signature.

## Exact boundary

This chamber is not the full nef cone.  Its ten primitive \((-4)\) mirrors
cannot be walls cut out by irreducible K3 curves, and additional \((-2)\)
roots may refine it.  After Weyl-transporting the marking to an ample chamber,
K3 Riemann--Roch orients the forty \((-2)\) classes as effective divisors;
they may be reducible and are not asserted to be all curve walls.

No surface automorphism is asserted to preserve the transported 50-wall
chamber.  The construction needs only algebraic divisor classes and their
integral reflections.

The underlying abstract `H13/50` polytope and the
\(U+D_8(-1)+D_4(-1)\) K3 family were already known.  The positive gain is the
explicit equation-to-lattice-to-fifty-mirror attachment.

## Reproduction

- [certificate generator](../scripts/certify_ss2_2a_h13_50_ud8d4_weierstrass_carrier.py)
- [exact JSON certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_weierstrass_carrier.json)
- [focused tests](../tests/test_ss2_2a_h13_50_ud8d4_weierstrass_carrier.py)
- [prior scaled-dual carrier](ss2_2a_h13_50_ud8d4_divisor_carrier.md)

Run from the repository root:

```bash
PYTHONPATH=src:. python3 -m scripts.certify_ss2_2a_h13_50_ud8d4_weierstrass_carrier --check-archive
PYTHONPATH=src:. pytest -q tests/test_ss2_2a_h13_50_ud8d4_weierstrass_carrier.py
```

For an optional byte-level replay of the audited primary preprint:

```bash
PYTHONPATH=src:. python3 -m scripts.certify_ss2_2a_h13_50_ud8d4_weierstrass_carrier \
  --primary-pdf /path/to/2009.09635v6.pdf
```
