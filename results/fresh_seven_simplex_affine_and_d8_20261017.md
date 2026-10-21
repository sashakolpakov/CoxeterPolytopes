# Four remaining seven-simplex stars with arbitrary supporting whites

Each of the four original-wall-star cases below has at most three ideal
vertices. This concerns actual seven-dimensional simplex faces of finite-sided
finite-volume hyperbolic Coxeter polytopes. Induced face angles need not be
Coxeter angles, and the supporting white block can have arbitrarily many roots.

This package completes the four cases left separate in the
[175-case analysis](fresh_finite_edge_seven_simplex_20261017.md). Alone it is
a four-case lemma, not an ambient dimension bound. The
[exact checker](../scripts/fresh_seven_simplex_affine_and_d8_20261017.py) and
[saved certificate](fresh_seven_simplex_affine_and_d8_20261017.json) use the
frozen standalone Fraction arithmetic over Q(sqrt2). Complete attachment
profiles are inputs from the
[exact classical census](fresh_high_rank_quasilanner_census_20261017.md).

Labels L20,L21,L22,L23 refer to that census. The eight black roots B are
original facet normals of the simplex face; W is its positive-definite
white supporting block. Every black deletion in B union W is PSD, and
is positive definite exactly when that actual vertex is finite. The
[original-star argument](fresh_multi_ideal_triangle_scope_20261017.md)
justifies this without assuming simple ideal vertices.

| Minimal diagram | Internal whites | Black block | Ideal-vertex bound |
|---|---|---|---:|
| L20 | 8 | affine B7 | 1 |
| L21 | 8 | affine D7 | 2 |
| L22 | 8 | affine E7 | 3 |
| L23 | 8,9 | finite D8 | 3 |

## General reduction for an affine black block

Normalize every root to squared norm2. If a white root u touches an affine
black block B, its pairing with B's positive null combination is strictly
negative, so B+u is indefinite. Deleting u leaves affine B; deleting any
black root leaves a principal block of an actual vertex Gram and is PSD.
Thus B+u is a minimal indefinite diagram of order9. The complete census
supplies every attachment, including all Coxeter labels:

| B | Supports of nonzero white-to-black contacts |
|---|---|
| affine B7 | {0}, {1} |
| affine D7 | {0}, {1}, {6}, {7} |
| affine E7 | {4}, {7}, {4,7} |

Every listed contact has magnitude1, hence label3. The certificate checks
the actual census provenance maps against each literal black Gram.

Let U be all whites touching B, and Z all other whites. The set U is
nonempty. At black deletion b, let C_b be the inverse Gram of B minus b,
and p_i the restricted nonnegative contact vector of direct white i. Put

\[
 (K_b)_{ij}=p_i^tC_bp_j,\qquad H_b=G_U-K_b.
\]

The finite Gram inverse C_b is entrywise nonnegative. Eliminating the
arbitrary remaining white block Z subtracts

\[
 D=R^tG_Z^{-1}R
\]

from every H_b, where R is the Z-to-U contact matrix. The matrix D is PSD
and entrywise nonnegative: R is nonpositive and the inverse of the
positive-definite white Z-matrix is nonnegative. Every actual vertex
condition is exactly H_b-D PSD; that vertex is finite exactly when this
residual is positive definite. Empty Z means D=0. This includes arbitrary
white trees and indirect contacts, with no rank or attachment-length cutoff.

## Complete direct-white configurations

For a single profile i, the largest projection cost (K_b)_ii is7/4 in
B7 and D7,3/2 for either single-leaf E7 profile, and2 for the double-leaf
E7 profile. All exceed1. Two copies of one profile already have a negative
two-by-two residual determinant at a maximizing vertex even at native
right angle. A nonright contact only worsens it. Thus profiles cannot repeat.

For different profiles the pair determinant is

\[
 (2-K_{ii})(2-K_{jj})-(a+K_{ij})^2,
\]

where a is the nonnegative magnitude of their native contact. A nonzero
Coxeter contact has a>=1. For each pair the certificate supplies a vertex
where this determinant is negative at a=0, or zero at a=0 or a=1.
The indicated K_ij is positive, so larger a is excluded strictly.
This covers every Coxeter label. The complete list is:

* B7: a single profile.
* D7: one profile, or one from {0,1} and one from {6,7}; a pair is right.
* E7: one profile, or the two single-leaf profiles {4},{7}, with native
  contact magnitude0 or1. The double-leaf profile {4,7} occurs alone.

Every three distinct profiles contain an excluded pair, excluding larger U.

For one direct white, D is a nonnegative scalar. Its maximum is the
minimum of H_b=2-(K_b)_11. Only minimizing deletions can become ideal.
There is one minimum for B7, two for D7, and one for each E7 single-profile
configuration. These give the stated single-profile bounds.

For a D7 pair, two black deletions are already ideal. At one of them H_b
has a strictly positive null vector h. PSD of H_b-D implies h^tDh=0;
entrywise nonnegativity forces D=0. The ideal count stays two. For the
E7 pair of native magnitude1, the same argument forces D=0; its three
ideal black deletions are1,4,7.

## The E7 pair with native right angle

Write the remaining correction as

\[
 D=\begin{pmatrix}x&z\\z&y\end{pmatrix}.
\]

At deletions4 and7, the uncorrected residuals are diag(2,1/2) and
diag(1/2,2). Thus0<=x,y<=1/2. PSD and entrywise nonnegativity of D imply
0<=z<=sqrt(xy)<=1/2.

At deletions0,2,3,5,6 the residual is diagonal before correction.
Its determinant lower bounds after correction are respectively

\[
 5/16,\quad11/36,\quad1/4,\quad11/36,\quad1/4.
\]

These and the diagonal lower bounds are strictly positive, so these
five vertices stay finite. At deletion1 the uncorrected residual is

\[
 H_1=\begin{pmatrix}9/8&-1/8\\-1/8&9/8\end{pmatrix}.
\]

After subtracting D its diagonals are at least5/8 and its offdiagonal
magnitude at most5/8. Its determinant could vanish only if x=y=z=1/2.
That choice gives determinant -1/4 at deletion4, violating its vertex
condition. Thus deletion1 is also finite. Only4 and7 can be ideal.
This completes all three affine-black cases.

## Finite D8 and the smaller-diagram alternative

For L23 with internal whites8,9 the black roots0,...,7 form finite D8.
Suppose some white u makes B+u indefinite. Every black deletion is PSD,
so a minimal indefinite subdiagram must contain all eight black roots.
It is the full order-nine diagram. Its affine deletion cannot be u,
since B is positive definite.

The checker separately reconstructs all31 order-nine cases whose sole
white deletion is positive definite. With two black affine deletions,
all further whites are orthogonal to the minimal diagram. With one,
eliminating further whites decreases only that black diagonal. Exact
minimum thresholds and ties give at most three ideals. This is the
positive-kernel/rank-one argument of the
[eight-simplex proof](fresh_finite_edge_cusp_tradeoff12_20261017.md), with
the smaller black set. It needs no unproved hard case.

## Finite D8 when every B+u is PSD

Assume instead B+u PSD for every white u. Adjoin existing internal white9
to B, forming affine E8, denoted F. The standard E8 numbering maps to L23 as

\[
 (0,1,2,3,4,5,6,7,8)\longmapsto(2,0,1,9,3,4,5,6,7).
\]

The internal white is standard node3. The other existing white, old8,
has profile B: magnitude sqrt2 at standard node8.

Any additional white u touching F makes F+u indefinite. Black deletions
are PSD by incidence; deleting internal white3 gives B+u, PSD by assumption;
deleting u leaves affine F. Thus F+u is minimal of order10. Its complete
E8 attachment list consists of the three profiles A,B,C in the eight-simplex
proof: magnitude1 at8, magnitude sqrt2 at8, or magnitude1 at7.

At black deletion1 the existing B-profile cannot coexist with any second
direct white. At native right angle the residual pair determinants for a
second A,B,C are respectively

\[
 -4/3,\qquad-28/9,\qquad-10/9.
\]

Nonright contacts worsen them. Thus the existing B-root is the only white
outside F touching F. Every other white Z is orthogonal to F. Eliminating
Z leaves a scalar K=2-D>0 for the sole direct white. At black deletion j
its residual is K-c_j, with projection costs

| Standard black j | 0 | 1 | 2 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| c_j | 5/3 | 16/9 | 7/4 | 8/5 | 3/2 | 4/3 | 1 | 0 |

Only1 attains the maximum16/9. PSD there forces strict positivity at all
other black deletions, giving at most one ideal vertex. Together with
the smaller-diagram alternative this case has at most three.

## Verification and scope

**FINAL/FROZEN.** Generation and a fresh saved replay both passed with
exit0. Source hashes, syntax and local links were checked. The separate
frontier_literature agent reviewed the arbitrary-white reductions and
independently reconstructed the critical E7 bounds and D8/E8 map,
costs and exclusion determinants with SymPy; that internal review passed.
The checker binds its literal census inputs and frozen arithmetic
implementation, checks all direct profile pairs, all31 smaller-diagram
alternatives, the E7 right-pair inequalities and the D8/E8 identification.
Arbitrary-white reductions are written proofs; this is not proof-assistant
verification. No external review, new ambient dimension ceiling or novelty
claim is made. Predecessor files were preserved.

