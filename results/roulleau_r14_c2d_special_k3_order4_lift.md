# TOOT: a genuine order-four K3 realization of the R14 `H11/19` section

This is a positive result.  There is a complex projective K3 surface (X)
with

\[
 \operatorname{NS}(X)=N=U\oplus E_8(-1)\oplus A_1(-1)^{\oplus4}
\]

and a purely non-symplectic automorphism (sigma) of order four whose action
on (N) is the double transposition used in
`R14-UE8A1x4-C2D-H11-19`.  Consequently the archived finite-volume
noncompact (H^{11}) Coxeter polytope with 19 facets is literally the fixed
section of a K3 nef chamber under a genuine surface automorphism.

This is a new project-level algebraic realization of the section, not a new
abstract Coxeter polytope, dimension, or arithmetic commensurability class.

## Exact Picard action

An explicit unimodular change of basis puts Roulleau's rank-14 source pair in
the form

\[
 N=U\oplus E_8(-1)\oplus A_1(-1)^4,
 \qquad
 \sigma^*|_N=1_{U\oplus E_8(-1)}\oplus(12)(34).
\]

Thus

\[
 N^{\sigma}=U\oplus E_8(-1)\oplus\langle-4\rangle^{\oplus2},
 \qquad
 N_{-}=\langle-4\rangle^{\oplus2}.
\]

The fixed lattice has rank 12, signature ((1,11)), determinant (-16), and
Smith diagonal (1^{10},4,4), exactly as in the archived C2D section.  The
anti-invariant lattice has no vector of square (-2); this is the decisive
Kähler gate that failed in the attempted normal-(V_4) lift.

## K3-lattice glue

Take

\[
 T=U^{\oplus2}\oplus A_1(-1)^{\oplus4},
 \qquad \operatorname{sig}(T)=(2,6),\quad \det T=16.
\]

On each rank-four block use

\[
 J=
 \begin{pmatrix}
 0&0&0&-1\\
 0&0&-1&0\\
 0&1&0&0\\
 1&0&0&0
 \end{pmatrix},
 \qquad J^2=-1,
\]

and set (g_T=J\oplus J).  Modulo two, the second block is the required
double transposition on the four discriminant coordinates.  The certificate
checks all 16 values of the discriminant quadratic form and the equivariance
of the anti-isometry

\[
 P=
 \begin{pmatrix}
 1&1&1&1\\
 1&0&1&1\\
 1&1&0&1\\
 1&1&1&0
 \end{pmatrix}.
\]

Graph gluing (A_N) to (A_T) through (P) gives an integral even
unimodular lattice of signature ((3,19)).  Hence it is the K3 lattice.  The
glued isometry has

\[
 \chi_g(x)=(x-1)^{12}(x+1)^2(x^2+1)^4,
 \qquad g^4=1,
\]

and (g^2) fixes exactly the primitive copy of (N).

## An exact period with no Picard jump

The (i)-eigenspace of (g_T) has a basis whose Hermitian Gram matrix is

\[
 \operatorname{diag}(4,-4,-4,-4).
\]

Let (alpha>0) satisfy (alpha^4=2), and take the period coefficients

\[
 1,\qquad \frac{\alpha}{10},\qquad
 \frac{\alpha^2}{100},\qquad \frac{\alpha^3}{1000}.
\]

The period is isotropic, is an (i)-eigenvector, and has positive norm

\[
 4\left(1-\frac{\alpha^2}{100}-\frac{\alpha^4}{10000}
          -\frac{\alpha^6}{1000000}\right)
 > \frac{244949}{62500}>0.
\]

The polynomial (z^4-2) is irreducible over (mathbf Q(i)), so
(1,alpha,alpha^2,alpha^3) are linearly independent over
(mathbf Q(i)).  Pairing a rational vector of (T) with the period therefore
forces its pairing with every (i)-eigenbasis vector to vanish.  The
(i)- and (-i)-eigenbases together span (T\otimes\mathbf C), hence that
rational vector is zero.  It follows that the period has
(operatorname{NS}(X)=N) exactly, with Picard rank 14.

## Kähler cone and Strong Torelli

The archived graph action permutes all 27 walls of Roulleau's complete K3
nef chamber.  Averaging the archived ample class (D_{132}) gives the fixed
integral vector

\[
 (4,10,18,8,20,24,30,38,18,30,24,10,10,0),
\]

of square (528).  Its minimum pairing with the 27 walls is (2), so it is
strictly inside the chamber.  Thus the glued Hodge isometry preserves a
Kähler class.  Strong Torelli realizes it as a genuine order-four K3
automorphism.  This closes the step that a lattice or wall permutation alone
does not close.

## Vinberg finite-volume conclusion

The fixed chamber has the already certified literal Vinberg census:

| invariant | value |
|---|---:|
| hyperbolic dimension | 11 |
| facets | 19 |
| elliptic edges | 795 |
| finite vertices | 114 |
| cusps | 13 |
| bad edges | 0 |
| bad parabolic neighborhoods | 0 |

Every elliptic edge has exactly two finite/ideal endpoints, and every cusp
neighborhood passes the maximal-parabolic test.  Therefore the polytope is
finite-volume and noncompact, not merely a locally consistent Coxeter
diagram.

## Independent classification cross-check

The exhaustive Brandhorst--Hofmann order-four table has exactly two rows with
cohomological character ((r,l,m)=(12,2,4)).  Row `0.4.3.4` has fixed
determinant (-16) and Smith diagonal (1^{10},4,4); row `0.4.3.12` has
fixed determinant (-64).  The optional full replay constructs an explicit
integral conjugacy between row `0.4.3.4` and the R14 pair, including the split
(U\oplus E_8(-1)\oplus A_1(-1)^4) basis and the two A1 swaps.  See
[Brandhorst--Hofmann's archived dataset](https://doi.org/10.5281/zenodo.10410278).

This classification check is independent corroboration; the existence proof
above is already self-contained at the exact lattice-and-period level.

## Algebro-geometric normal form

The same order-four symmetry has the useful Weierstrass search chart

\[
 y^2=x^3+(a_1t+a_3t^3)x^2
       +(t^2-p^2)(t^2-q^2)x,
 \qquad
 \sigma(t,x,y)=(-t,-x,iy).
\]

Substitution gives (sigma^*F=-F), and
(sigma^2(t,x,y)=(t,x,-y)), elliptic inversion.  The visible section
((0,0)) is 2-torsion.  For example,

\[
 y^2=x^3+(t^3+t)x^2+(t^2-1)(t^2-4)x
\]

has an (I_4^*) fiber at infinity, four (I_2) fibers, and six (I_1)
fibers.  The trivial lattice (U\oplus D_8(-1)\oplus A_1(-1)^4) has
determinant (-64); adjoining the visible torsion section gives the expected
index-two even overlattice of determinant (-16).  This equation is recorded
as an algebraic search chart.  The exact Picard and Torelli assertions are
certified by the period construction, rather than by an unperformed
point-counting proof for that one numerical member.

## Reproduction

Run

```bash
PYTHONPATH=src:. python3 scripts/certify_roulleau_r14_c2d_special_k3_order4_lift.py --check-archive
```

For the independent 79-row classification replay, add

```bash
--brandhorst-hofmann-order4 /path/to/K3Groups/results/purely_ns/order04.txt
```

The exact certificate is
[`roulleau_r14_c2d_special_k3_order4_lift.json`](../data/enumerations/roulleau_r14_c2d_special_k3_order4_lift.json).
