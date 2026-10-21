# Two large-affine counterexamples to the contact-only-black reduction

The following proposed graph implication is false:

> Let A be a connected affine classical diagram of order N≥18, with two
> marked roots v1,v2. Suppose both (A−vi)+z are spherical and A+z contains
> no connected affine subdiagram of order at most17. Then z contacts only
> the two marked roots.

There are two explicit counterexample families for every N≥18. They do
not assert a finite-volume ambient realization. They show that the listed
local diagram hypotheses do not prove the proposed reduction; additional
conditions from an actual unpaid quadrilateral must be used.

Both families also satisfy a useful negative continuation test: if a new
root y is right to every node of A and has label3 contact with z, at least
one of the two required spherical deletions becomes nonspherical. This
does not yet prove nonextendibility for every possible white-contact
profile.

The [source](../scripts/fresh_large_affine_white_contact_counterexamples_20261017.py)
and [certificate](fresh_large_affine_white_contact_counterexamples_20261017.json)
use separate exact quadratic-field arithmetic, with no project arithmetic
imports. They check instances N=18,19,23,31 of each family, including
positive kernels, every displayed congruence pivot, literal finite forests,
negative witnesses, and the continuation failures. The all-rank arguments
below explain why the conclusion is not a finite-rank extrapolation.

## Common notation and a minimal-affine criterion

Let A be the N-node affine D diagram (conventionally tilde D_(N−1)). It
has a long central path with two leaves at each end. At the left fork let
p be the center, a,c its two leaves, and t the next node on the long path.
All bonds have label3. Its norm2 Gram is rational, with positive kernel
vector equal to1 on the central path and1/2 on all four leaves.

For any symmetric Z-matrix G with this connected nonzero graph and a
positive kernel vector nu, the matrix D_nu G D_nu is a connected weighted
graph Laplacian: its off-diagonal entries are−w_ij with w_ij>0 and its
diagonal entries are the incident-weight sums. Therefore G is positive
semidefinite of corank1, and every proper principal block is positive
definite. This proves the needed minimal-affine assertion directly.

In the source numbering, the central path is0,...,N−5, p=0,t=1,
a=N−4,c=N−3, the right leaves are N−2,N−1, and z=N.

## Family 1: a label5 leaf at a white fork leaf

Mark v1=p and v2=t. Add only the bond z–a with label5. Thus z has a
nonzero contact with the white node a and no contact with either marked
root.

Deleting p leaves

`I2(5) + A1 + D_(N−3)`.

Deleting t leaves

`H4 + D_(N−4)`,

where the H4 path is `z --5-- a --3-- p --3-- c`.
The two long D tails are proper principal blocks of the old affine D
matrix, hence positive definite by the Laplacian criterion. The I2(5)
and H4 blocks are positive definite by their exact leading minors. Thus
both entire marked deletions are spherical for every N≥18.

There is no connected affine subdiagram containing z. Indeed, any such
subdiagram must also contain a. If its old-node part B is a proper
principal block of A, then B is rational positive definite. Eliminating B
leaves the scalar

\[
 2-\tau\,(G_B^{-1})_{aa},\qquad
 \tau=4\cos^2(\pi/5)=(3+\sqrt5)/2.
\]

The inverse entry is a positive rational number, while tau is irrational,
so this scalar is nonzero. The whole block is nonsingular and therefore
not affine. If the old-node part is all A, its positive null vector has
nonzero product with z. The span of that null vector and z has negative
determinant, so the block is indefinite. Consequently the only connected
affine subdiagram of A+z is the original N-node A. In particular none has
order at most17.

This proof uses no assumption that a label5 bond occurs in a classified
affine diagram; the irrational Schur coefficient supplies a direct
exclusion.

For the continuation test, add y–z with label3 and y right to all A.
The deletion of t contains the four-node path

`y --3-- z --5-- a --3-- p`.

Its determinant is exactly `3−2sqrt5<0`. Hence that deletion is no longer
spherical. A stronger nonright contact to z also cannot restore positive
definiteness, by the nonnegative test-vector criterion for acute Grams.
Only the label3 test is needed in this package.

## Family 2: a black4/white3 fork switch

This family was proposed independently by the parent agent. Mark v1=a
and v2=p. Add z–a with label4 and z–c with label3. The contact to c is
again a white contact.

Deleting a produces a finite D_N path-and-fork diagram: the new long
endpoint runs z–c–p toward the right fork. This is positive definite,
for example as a proper principal block of the corresponding(N+1)-node
affine D diagram. Deleting p produces

`B3 + D_(N−3)`,

where the B3 path is `a --4-- z --3-- c`. The B3 leading principal minors
are positive, and the D tail is an old proper principal block. Thus both
marked deletions are spherical.

The induced cycle z,a,p,c has three label3 bonds and one label4 bond.
The all-ones vector has square `2−2sqrt2<0`. A positive-semidefinite
subdiagram containing z must therefore omit at least one of a,p,c.
The full a-deletion and p-deletion are spherical, as just shown.
The full c-deletion is an N-node affine B diagram: its long path has the
doubled terminal z–a at the left and the unchanged fork at the right.
A positive kernel vector has value1 along the main path including a,
sqrt2/2 at z, and1/2 at the two right leaves. The connected weighted
Laplacian identity proves that this block is minimal affine.

It follows that the only connected affine subdiagrams of A+z are the
original N-node D diagram and this N-node B diagram. Again none has
order at most17. This argument checks the entire possible affine
subdiagram domain, rather than testing a bounded list of small subsets.

For the continuation y–3–z, the p-deletion contains the four-node star
with center z and leaves a,c,y, of labels4,3,3. It is affine, with positive
kernel `(1,sqrt2/2,1/2,1/2)` in order(z,a,c,y). Thus the required deletion
cannot remain spherical.

A superficially similar single label4 leaf at a is not a counterexample:
it creates a small affine B diagram using p,a,c,t,z. The two-contact
switch above changes that induced subdiagram, and its no-small-affine
proof depends on the resulting negative cycle.

## Scope for the finite-quadrilateral transfer

Both full(N+1)-root Grams are Lorentzian of inertia(N,1,0). A checked
positive-definite N-root deletion and a negative determinant certify
this in the saved instances; in every rank, adjoining z to the old
positive-semidefinite corank1 matrix with nonzero null-vector coupling
has exactly one negative direction and no radical. The connected acute
Gram also has a strict timelike center. These are legitimate local
Coxeter configurations, not numerical diagram artifacts.

They are not claimed to complete to finite-volume chambers having only
large affine cusp factors. The original full-volume hypothesis may
impose further incidence restrictions. In particular, the label5 family
creates an H4 finite deletion involving a marked root; a zero-payment
assumption on the full quadrilateral can potentially exclude it. The
fork switch uses a contact to a marked root and cannot model an external
opposite black root required to be right to both marked roots.

The [private-tree lemma](fresh_private_affine_tree_exclusion_20261017.md)
assumes contact-only-marked-roots at its outset and is not contradicted.
What fails is deriving that assumption from the two spherical deletions
and the small-affine exclusion alone. The more appropriate remaining
question is whether every surviving white-contact profile is terminal
along a shared white path. Only the two families above have been checked
for that continuation in this package.

**FINAL/FROZEN:** generation and fresh saved replay both passed. The source
hash, local report links, and source syntax were checked. The displayed
finite instances verify the formulas; the written minimal-affine and
continuation arguments establish the two families in every N≥18.
