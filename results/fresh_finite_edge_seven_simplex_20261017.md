# Seven-simplex cases:162 affine-deletion cases and nine finite-black reductions

For an actual seven-dimensional simplex face S of a finite-volume acute
Coxeter polytope, assume S has at least two finite vertices. This package
proves that S has at most three ideal vertices in171 of the175 minimal
diagram cases. The four other cases are explicitly assigned to a separate
proof. Thus this package alone does not claim the complete seven-simplex
theorem, a new uniform triangle coefficient, or an ambient dimension bound.

The main additional result is an all-size treatment of nine cases with a
finite eight-root black block. It permits arbitrarily many supporting
white roots and arbitrary indirect white tails. No rank or label cutoff
is used in that treatment.

## Exhaustive partition

The [simplex-star theorem](fresh_multi_ideal_triangle_scope_20261017.md)
gives eight original black facet normals B. An inclusion-minimal
indefinite diagram L contains all eight. Its order is8,9, or10. The
[frozen high-rank census](fresh_high_rank_quasilanner_census_20261017.json)
has respectively4,4,3 such types. Choosing the internal white roots
therefore gives

\[
 4+4\binom91+3\binom{10}2=175
\]

cases. The [exact output](fresh_finite_edge_seven_simplex_20261017.json)
records every case, black/white choice, affine black deletion, applicable
threshold, and equality tie. The partition is:

| Case mechanism | Count |
|---|---:|
| At least two affine black deletions | 77 |
| Exactly one affine black deletion | 85 |
| No affine black deletion | 13 |

The13 cases in the last row, using literal census indices, are
L20/L21/L22 with internal white8; L23 with internal whites8,9; and L25
with internal whites a,9 for a=0,…,8. The first four are outside this
package's conclusion. All nine L25 cases are proved below.

## The162 cases with an affine black deletion

The [preceding proof](fresh_finite_edge_cusp_tradeoff12_20261017.md)
applies without change to one or two internal whites. For an affine
black deletion b, the positive kernel of L\{b}, PSD of its actual vertex
star, and nonpositive native contacts force every external white root to
be orthogonal to all of L\{b}.

Two such deletions cover L, forcing complete orthogonality between L and
the external whites. Their positive-definite block cannot change vertex
status. This handles77 cases.

For a unique affine black deletion b, eliminating the external white
block changes only the diagonal at b by −τ, τ≥0. For each other BLACK
deletion j, put H_j=G_{L\setminus\{j\}}. The necessary and sufficient
vertex PSD conditions are

\[
 \tau\le\frac1{(H_j^{-1})_{bb}}.                 \tag{1}
\]

The internal-white deletions are not included as vertex constraints.
At τ strictly below the minimum there is only the original ideal vertex.
At equality, precisely the tied minimizing black deletions become ideal.
The certificate recomputes every threshold and the matrix inertia at
that boundary.

Of the162 cases,130 have upper bound2 on ideal vertices,31 have upper
bound3, and one requires an explicit qualification. The exceptional
case L19 with internal white8 has seven simultaneous threshold ties.
Its allowed ideal counts are1 or8. The all-eight-ideal boundary is
outside the hypothesis that S has at least two finite vertices.
Consequently all162 cases under the stated hypothesis have at most
three ideal vertices.

This exception is a reason not to state an unqualified theorem that
every seven-simplex face has at most three ideal vertices.

## The nine finite-black cases

In L25, nodes0,…,8 form the following affine-E8 diagram F, with
norm-two roots and all displayed edges labeled3:

\[
 0-1,\qquad0-2-3,\qquad0-4-5-6-7-8.
\]

The chosen internal whites are a∈{0,…,8} and node9; the eight black roots
are B=F\{a}. They form a positive-definite Gram matrix. The original
white9 has only one F contact: label3 at node8. Call this profile A.

First suppose some supporting white root u has B∪{u} indefinite.
Deleting u leaves positive-definite B; deleting any black root leaves
a PSD principal submatrix of an actual vertex star. Therefore B∪{u}
is itself an order-nine minimal indefinite diagram containing every
black root. It must be one of the exact census types. An affine deletion
cannot be u because B is positive definite, so this new minimal diagram
has an affine BLACK deletion and belongs to the already proved162
cases. In particular it has at most three ideal vertices under f≥2.

It remains to assume

\[
 G_{B\cup\{u\}}\succeq0
 \quad\hbox{for every supporting white root }u.   \tag{2}
\]

Let U be all whites outside F having a nonzero contact with F, and let Z
be all remaining whites outside F. In particular U contains original
white9, and Z is orthogonal to the entire F. For every u∈U, the positive
affine kernel of F and its nonzero nonpositive contact with u make F∪{u}
indefinite. Deleting u gives affine F; deleting a gives the PSD matrix
in(2); deleting any black root gives a PSD actual-vertex submatrix.
Thus F∪{u} has all proper principal matrices PSD. The complete attachment
analysis can also be read directly from its deletion inverses below.

Deleting F node8 gives inverse diagonal

\[
 (30,8,14,4,20,12,6,2)
\]

on nodes0,…,7. All entries of these finite Gram inverses are nonnegative.
Every nonzero norm-two Coxeter contact has magnitude at least1. The
vertex/deletion budget d^tG^{-1}d≤2 therefore forces contacts at all nodes
except7,8 to vanish, and a node7 contact has magnitude exactly1.
Deleting F node3 gives inverse block

\[
 \begin{pmatrix}2&1\\1&1\end{pmatrix}
\]

on7,8. Hence both cannot be contacted, and a sole node8 contact has
magnitude at most sqrt2. The complete possibilities for U are

| Profile | Contact |
|---|---|
| A | node8, magnitude1 |
| B | node8, magnitude sqrt2 |
| C | node7, magnitude1 |

If the deleted node in this argument is the internal white a, the PSD
input is exactly assumption(2). Otherwise it is the actual black vertex
condition. Thus the list is complete for every a, without an arbitrary
contact alphabet or rank cutoff.

## The existing A root excludes profiles B and C

If a≠3, actual black deletion3 is available. Its projection-cost Gram
for profiles A,B,C is

\[
 \begin{pmatrix}
 1&\sqrt2&1\\\sqrt2&2&\sqrt2\\1&\sqrt2&2
 \end{pmatrix}.
\]

A B or C root has zero residual diagonal while its residual cross product
with the already existing A root is strictly negative. This contradicts
PSD. At native right contact, the A+B and A+C determinants are−2 and−1;
a negative native contact only makes them smaller.

If a=3, use actual black deletion1. Its finite inverse block on nodes7,8
is

\[
 \begin{pmatrix}14/9&7/9\\7/9&8/9\end{pmatrix}.
\]

The A+B and A+C residual determinants at native right contact are
respectively−4/3 and−1/9. Again every permitted nonpositive native product
only worsens these determinants. Thus B and C are impossible here too.

In all nine cases, every member of U has the same profile A. This uses
the fixed existing A root; it would be incorrect to omit that hypothesis
in a different star.

## Arbitrary white tails reduce to one positive-definite matrix

Let G_U,G_Z be the native white Gram blocks and R their cross products.
Eliminating Z gives

\[
 K=G_U-R^tG_Z^{-1}R\succ0.                        \tag{3}
\]

Positive definiteness holds because U∪Z is a subset of the supporting
white set W. No assumption on the number, graph, or label sizes of Z or
U is needed. For empty Z use K=G_U.

For any actual black deletion j≠a, eliminate the finite block F\{j}.
All rows of its U-contact matrix agree, since every U root has profile A.
Since Z is orthogonal to F, eliminating it commutes with this projection.
The residual vertex matrix is

\[
 K-c_j\mathbf1\mathbf1^t,
\qquad
 (c_0,\ldots,c_8)=
 \left(\frac56,\frac89,\frac78,1,
       \frac45,\frac34,\frac23,\frac12,0\right). \tag{4}
\]

The zero at j8 occurs because the sole contacted node is deleted. Each
coefficient is verified from a complete exact inverse. Among the actual
black vertices j≠a, the maximum is unique: j3 if a≠3, and j1 if a=3.

Let c_* be this maximum. The actual vertex condition gives
K−c_*11^t PSD. For c<c_*,

\[
 K-c11^t=(K-c_*11^t)+(c_*-c)11^t
\]

is positive definite. Indeed, a zero quadratic vector would be in the
kernel of the first PSD summand and orthogonal to1. It would then have
zero quadratic form under K, contradicting K>0. Therefore every actual
black vertex other than the unique maximizing one is finite. This
entire branch has at most one ideal vertex.

Combining this with the initial order-nine reduction proves at most
three ideal vertices, under f≥2, in all nine L25 cases. The proof covers
arbitrarily long white tails: they affect K but not the strictly ordered
scalar costs. It never infers a bound on supporting whites from the
number of nodes of a minimal indefinite subdiagram.

## Scope and certificate

The remaining cases L20/L21/L22 with white8 and L23 with whites8,9 are
not concluded here. A complete coefficient improvement below12 requires
those additional proofs. The current package's output explicitly lists
them as outstanding, even if another package later resolves them.

The [checker](../scripts/fresh_finite_edge_seven_simplex_20261017.py)
uses the frozen standalone Fraction Q(sqrt2) matrix engine from the
tradeoff12 certificate. It imports no arithmetic from the producer of
the high-rank census. It reconstructs every required matrix and all175
choices, saves all162 easy-case thresholds/statuses, checks the nine
finite-black inverses and both pairs of strict determinant obstructions,
and verifies the unique projection-cost maxima. It binds the frozen
census, finite-ended-port identity, simplex-star argument and prior
coefficient12 proof. The geometric reductions above are written proofs,
not finite tests of white tails.

Generation and fresh saved replay PASS. Syntax, bound source hashes and
local links PASS. No predecessor file is modified. This is a scoped
internally checked lemma, with no ambient dimension or full seven-simplex
classification claim and no assertion of external peer review or formal
proof-assistant verification.
