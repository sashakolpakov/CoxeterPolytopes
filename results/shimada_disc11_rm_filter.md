# Shimada's Discriminant-11 K3 Data under the Real-Multiplication Filter

Status: the fixed real-multiplication eigenspace has empty intersection with
the K3 nef cone, proved by two effective divisor walls.  A separate exact
calculation has since found a nonempty conjugate marking; see
`results/shimada_disc11_conjugate_marking.md`.

## Exact surface match

Shimada's No. 5 singular K3 surface has transcendental Gram matrix

\[
T_X=\begin{pmatrix}2&1\\1&6\end{pmatrix},
\]

so it is exactly the singular K3 surface underlying the
`coupled_discriminant11` candidate, not merely another surface of the same
Picard rank.  Shimada computes 1,098 induced Borcherds chambers (maximum level
13), 767 generators of the image of \(\operatorname{Aut}(X)\), and an upper
bound of 347 orbits of smooth rational curves.  The computational supplement
contains 347 explicit \((-2)\)-vectors in its set `B` [Shimada2015,
ShimadaDisc11Data].

## Certified basis translation

The published Néron--Severi Gram matrix is

\[
U_\phi\oplus -T_X\oplus E_8(-1)^2,
\qquad
U_\phi=\begin{pmatrix}0&1\\1&-2\end{pmatrix}.
\]

The importer gives an integral determinant-one basis translation from this
form to the interleaved \((e_i,\phi e_i)\) trace basis of \(H_{11}\).  It
verifies both the published Gram matrix and all 347 trace squares exactly.
The source file used in the calculation has SHA-256

```text
fa39e5b59f4a4b2b1b7fe7632e5737d1734ee3eb7e7023ae26bf2b83387afb68
```

## Filter result

Of the 347 published curve representatives:

- 17 have field norm positive at both real places;
- all 17 satisfy the \(\mathcal O_F\)-reflection divisibility condition;
- eight have norm \(2\), and nine have norm \(4-2\phi\);
- the other three possible K3 wall shells do not occur in this representative
  set; and
- the 17 vectors have field span rank 9, not 10.

The seed has a simple exact description: it consists of an \(E_8\) simple
basis, its scalar multiple by \(\phi-1\), and one root in the coupled binary
block.  Among its 136 unordered pairs, 28 field pairings are negative, 100
are zero, and 8 are positive.  Thus these vectors are curve representatives
from different parts of the nef tessellation, not an outward-oriented simple
root system for one Coxeter chamber.  Rank 9 alone already prevents this set
from bounding a finite-volume polyhedron in \(\mathbb H^9\).

Applying each of the 767 published automorphism generators once to each of the
17 roots yields 2,567 distinct effective \((-2)\)-classes.  Exact refiltering
produces no additional field-reflection root: the surviving set is still the
same 17 vectors of field span rank 9.  None of those 767 displayed generators
individually commutes with the integral multiplication-by-\(\phi\) operator.

## Two-divisor empty-section certificate

Put \(\alpha=-\phi/\sqrt5\), let \(h\) be the field form, and let
\(q=\sigma(\alpha)h_\sigma\) be the form of signature \((1,9)\) on the
identity eigenspace.  With \(e=(1,0,\ldots,0)\), exact calculation gives

\[
q(e,e)=\frac25+\frac65\phi>0.
\]

Two oriented effective roots in Shimada's published set `B` project to

\[
d_+=(1+\phi,-1+2\phi,0^8)
\]

and

\[
d_-=(-3-\phi,1-3\phi,-2,-4,-6,-5,-4,-3,-2,-3).
\]

Their exact field norms and control-vector pairings are

\[
\begin{aligned}
h(d_+,d_+)&=6-4\phi,
&q(e,d_+)&=\frac85+\frac95\phi>0,\\
h(d_-,d_-)&=8-6\phi,
&q(e,d_-)&=-\frac{11}{5}-\frac{18}{5}\phi<0.
\end{aligned}
\]

Moreover

\[
q(d_+,d_+)=-\frac85+\frac65\phi>0,
\qquad
q(d_-,d_-)=-2+2\phi>0,
\]

and both trace squares are \(-2\).  Thus \(d_+\) and \(d_-\) are timelike and
have opposite time orientation.  Every positive \(x\) in the component of
\(e\) has \(q(x,d_-)<0\), while every positive \(x\) in the other component
has \(q(x,d_+)<0\), by reverse Lorentzian Cauchy--Schwarz.  Shimada's
orientation makes both classes effective, so nefness requires both pairings
to be nonnegative.  Therefore

\[
E_\sigma\cap\mathcal P\cap\operatorname{Nef}(X)=\varnothing.
\]

This is a complete obstruction for the fixed self-adjoint
\(\mathbb Z[\phi]\)-marking; computing its automorphism centralizer cannot
repair the empty intersection.  It does **not** prove the field lattice
nonreflective.  Conjugating the marking by a K3 Weyl isometry changes the
eigenspace relative to the nef cone, so such markings require a separate,
certified chamber-walk analysis.  The first successful chamber-interior
certificate is now archived for chamber 68 in
`results/shimada_disc11_conjugate_marking.md`.

## Reproduction

Download Shimada's `threediscs.zip`, extract `disc11.txt`, and run

```sh
PYTHONPATH=src:. python3 scripts/analyze_shimada_disc11.py --summary-only /path/to/disc11.txt
python3 -m pytest tests/test_shimada_disc11_import.py
```

The external 47 MB text file is not vendored; its URL and checksum are recorded
in `data/candidates/shimada_disc11_rm_filter.json`.
