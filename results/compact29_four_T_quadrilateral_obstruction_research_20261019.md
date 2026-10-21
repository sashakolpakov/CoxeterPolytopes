# Compact29 frontier: an exact four-T quadrilateral obstruction

19 October 2026. **New written local model, awaiting independent review.**
No numerical search, finite census, optimizer, checker, Lean or build ran.
This does not construct a compact or finite-volume H29 chamber or exclude
one. It supplies a concrete countermodel to removing T-to-T quadrilateral
consumers from the remaining refund problem using only finite corners,
native Coxeter products, Lorentz signature and absence of affine subsets.

## 1. Exact remaining frontier

The reviewed all-rank4 refund leaves98 typed rank29 forests:

* 80 of rank pattern7+7+7+8, with ranks7/8 types A,B,D,E and original
  allowance87;
* 18 of rank pattern9+10+10, types A,B,D, with original allowance85.

After the compatible E7 facewise refund,79 have allowance87, the special
T=E7+E7+E7+E8 has allowance171/2 (at most81 when active), and18
have allowance85. The threshold remains84. These facts are inherited
from `compact29_all_rank4_leaf_refund_20261018.md`; they do not exclude
the98 types. No arithmetic-field assumption is permitted in this lane.

The inactive T refund's only half-unit loss is the possibility that a
quadrilateral's paid corners are all eligible T angles. The following
model shows that even FOUR T corners survive the local tests above.

## 2. Native nine-wall diagram

Use old E7 roots r_0,...,r_6, all of norm two, with edges

    0--1--2--3--4--5,   2--6,

all simple. Put a=r_3, b=r_5, s=r_4 and xi=sqrt(2). Add norm-two
roots c,d with the following COMPLETE nonzero new products:

    c.r_1=-1,    c.s=-1,    c.a=-xi,
    d.s=-1,      d.b=-xi.

All other products of a new root with an old root are zero, and c.d=0.
Every pair is a finite native Coxeter pair: labels2,3,4 only.
Let W={r_0,r_1,r_2,r_6,s}=A4+A1. The cyclic black order is a,b,c,d.
All four sets W plus adjacent black roots are E7:

* (a,b): the old E7;
* (b,c): fork r_1 has arms r_0; r_2--r_6; c--s--b,
  of lengths1,2,3;
* (c,d): the same fork has arms r_0; r_2--r_6; c--s--d;
* (d,a): fork r_2 has arms r_6; r_1--r_0; a--s--d.

In each corner the marked black pair lies in this connected E7, so
its original Vinberg weight is one. Appending untouched2E7+E8 makes
all four full finite corners have type T. The old whole-component
assumptions are retained exactly.

## 3. Exact white projection and compact quadrilateral

In the order (0,1,2,6), the inverse A4 Gram matrix is

    (1/5) [[4,3,2,1],
           [3,6,4,2],
           [2,4,6,3],
           [1,2,3,4]],

and the inverse Gram entry at s is1/2. Eliminating W gives the
following four-side Schur matrix:

    R = [[3/10, -1/2, -(xi+13/10), -1/2],
         [-1/2, 3/2, -1/2, -(xi+1/2)],
         [-(xi+13/10), -1/2, 3/10, -1/2],
         [-1/2, -(xi+1/2), -1/2, 3/2]].

For example a's white projection has square6/5+1/2=17/10,
as does c's; their cross projection is4/5+1/2=13/10.
The b,d projections each have square1/2. Each adjacent cross
projection is1/2. These give all sixteen entries directly.

Here is an explicit realization in Lorentz signature(-,+,+). Set

    t²=(xi+1)/2,      A²=(5xi+8)/10,
    u²=(xi-1)/2,      B²=(xi+2)/2,

using positive square roots, and choose projected outward normals

    a_bar=(t,A,0),    c_bar=(t,-A,0),
    b_bar=(u,0,B),    d_bar=(u,0,-B).

Their Gram matrix is exactly R: t*u=1/2, their squares are3/10
and3/2, and their opposite products have the displayed magnitudes.
Thus R has inertia(2,1,1), listing positive, negative and zero counts,
and rank three. The full nine-wall Gram matrix has inertia(7,1,1)
and rank eight. No numerical eigenvalue inference is involved.

In the Klein chart with time coordinate one their intersection is

    |X|<=t/A,    |Y|<=u/B.

Its four vertices lie STRICTLY inside the unit disk, because

    t²*u²=1/4 < (3/10)*(3/2)=9/20,

which is equivalent to (t/A)²+(u/B)²<1. Hence this is an actual
compact hyperbolic quadrilateral, not just a corank condition. All
four side inequalities are nonredundant and each nonincident side
inequality is strict at a vertex. The center (1,0,0) is strict on
all four black halfspaces. The induced angles of this intrinsic
quadrilateral need not be Coxeter angles.

Lift the normals by their stated W projections. W is positive
definite, so this realizes the full native Gram matrix in signature
(7,1). The face lies on W=0. A sufficiently small negative white-dual
perturbation of the strict face center makes every W inequality
strict while preserving all four black inequalities. Keeping one
white equality, and perturbing the others, similarly realizes each
white wall. The finite corner points already realize each black
wall. Thus the nine walls are actual facets of a nonempty acute
local chamber and the quadrilateral is an actual two-face. The
untouched22-dimensional spherical summand gives the H29 local model.

This construction also fixes the time orientation of all four finite
corners explicitly. It does not rely on a determinant branch with
incorrectly oriented new roots.

## 4. No affine principal subdiagram

This can be checked by cases, without a subset census. Any subset
with neither opposite black pair {a,c} nor {b,d} lies in one of
the positive-definite E7 corners. If it contains s and an opposite
pair, it has a negative principal three-node form: the triangle on
that pair and s has products -xi,-1,-1, so its all-ones vector
has square2-2xi<0. Such a subset is not positive semidefinite.

Suppose s is absent. The {b,d} component is just I2(4), possibly
with one root removed, and is separate from all other roots.
For a subset containing a,c, if it includes r_1,r_2, the cycle
a--c--r_1--r_2--a has one4-bond and three3-bonds. Its all-ones
vector again has square2-2xi<0. If r_1 is absent, the component
containing a,c is a subdiagram of the B4 chain c--a--r_2--r_6;
r_0 is isolated. If r_2 is absent, it is a subdiagram of the B4
chain a--c--r_1--r_0; r_6 is isolated. All these are positive
definite. This covers every subset, including ones containing both
opposite pairs. Every positive-semidefinite principal block is
therefore positive definite. No affine principal diagram occurs.
Appending finite orthogonal summands does not alter this conclusion.

## 5. Consequence and precise limit

All four old payments on this quad equal one and are eligible at T.
The existing symmetric facewise refund sets each to1/4, refunding3/4
per corner; this is not the worst half-refund case, but it proves
that eligible T-to-T quad incidence cannot simply be excluded. The
face simultaneously consumes payments at four full heavy-type corners.
An argument that deletes their units independently would fail the
quad inequality. Directed redistribution or a global incidence proof
would need to account for this configuration.

The model does not establish that this nine-wall local chamber has
finite volume, that it can be compactly completed while preserving
the face, or that several such faces coexist around one T vertex.
No claim about an actual compact H29 polytope follows. In particular
the current98-type frontier and compact dimension ceiling29 remain
unchanged. Its useful content is a concrete obstruction to an overly
strong local no-T-to-T-quad lemma, together with a small exact diagram
on which any proposed global payment improvement can be tested.
