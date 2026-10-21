# A uniform ideal-four transfer for affine A cycles of order at least26

21 October 2026. **New written proof, held for independent review.**
Source/manual mathematics only; no enumeration, scientific, symbolic,
numerical, Lean or TeX job, manuscript or canonical edit. No new dimension
ceiling or all-family ideal-four envelope is claimed.

**Theorem.** In the current finite-volume57 actual-face framework, every
distance-five unit in a whole varying affine A cycle of order m>=26 can
be replaced by zero. Each affected quadrilateral has a different actual
corner paid by an unchanged full unit: either an adjacent black bond,
or a whole finite component with marked distance at most5. An affected
triangle is impossible. This theorem allows arbitrary supporting whites,
spectator ranks and native contact magnitudes.

Consequently the A-family factor itself can use precisely distance<=4
at every order m>=26, with total4m rather than5m. This decreases its
actual fiber charge by V for each such factor. The other affine families
and orders10..25 are not changed, so no all-family envelope or ambient
bound follows. The current negative-margin low-rank cusp profiles are
not removed by this theorem.

## 1. Actual profiles, separate edge budgets and native heights

Use square-two normals and the actual quadrilateral star W+a+b+c+d,
with opposite pairs {a,c},{b,d}. Its chosen ideal factor is a whole
cycle F, and the remaining supporting forest Z is positive definite
and orthogonal to F. Set a=0,b=5. The two old arcs have lengths
P=5 and Q=m-5>=21. Each W+one-black edge Gram is strictly positive
definite; each corner Gram is positive semidefinite; the complete star
has one negative direction and corank one.

Write x,y for the nonnegative contact rows of c,d on Z. Their actual
separate edge conditions are retained, individually:

    C_F(c)+x^t G_Z^-1 x<2,
    C_F(d)+y^t G_Z^-1 y<2.

There is no replacement by a shared budget or by a bounded spectator
list. The proof below is stronger than these two inequalities: after
their actual geometric role is established, it needs only the
nonnegativity of the normalized spectator difference energy.

If b.c, a.d or c.d is negative, that adjacent marked pair already gives
an unchanged unit. Suppose all three are zero. Let nu be the positive
old cycle kernel and h=-c.nu>0,k=-d.nu>0. Deleting a gives A_(m-1).
Its inverse diagonal at position j is j(m-j)/m. At m>=26 only the
first two positions at either end can have a native contact of cost
at most2. An endpoint allows magnitude1 or sqrt2; a second node
allows only1. Two contacts allow only the two unit endpoints; their
cost is exactly2. All other pairs and all triples exceed2. This is
also exactly the large-path finite/affine one-root classification:
isolated, end3, end4, second-node3 or both-end3. An omitted-mark
contact remains unrestricted by this principal test.

A single port on the short P arc gives a marked path from c to b of
length5 at an endpoint, or4 at its second node. At the adjacent corner
F-a+c it lies in a whole finite component. Indeed the opposite old
endpoint is a fixed simple leaf with a fixed degree-two neighbor;
Z attaches only through c and cannot change those two degrees. No
classical affine diagram has that leaf configuration, while the old
A_(m-1) path is too large for an exceptional affine factor. Every
remaining component is a positive-definite principal part of Z.
Thus this whole actual corner is finite, and its distance<=5 unit is
unaffected by all current deep deletions, protected-leaf changes and
context changes. The same argument applies to d.

It remains only to exclude profiles with no such single short-arc port.
Put u=1/h. Their moments l=(l1,l2), norm constant c0 and native bounds
are as follows; the second coordinate denotes the long Q arc.

| c profile | l | c0 | bound on u |
| --- | --- | ---: | ---: |
| isolated | (0,0) | 2 | u<=1 |
| end3 on Q | (0,1) | 1 | u<=1 |
| end4 on Q | (0,sqrt2) | 0 | u<=1/sqrt2 |
| second-node3 on Q | (0,2) | 0 | u<=1 |
| both-end3 | (1,1) | 0 | u<=1/2 |

These bounds use h equal to the sum of the contacts to F. In the
isolated case its positive opposite-mark contact is a native magnitude
at least1. In the other rows that contact is nonnegative and the shown
ports already contribute1,sqrt2,1 or2. Nothing assumes h=k, bounds
an omitted contact from above or changes its Coxeter alphabet.

## 2. Exact all-rank determinant inequality

Put A=1/P=1/5, B=1/Q<=1/21. For the d-profile write m,d0,v=1/k.
Set L=u*l1+v*m1, M=u*l2+v*m2 and C=c0*u^2+d0*v^2.
The inherited exact Schur identity is

    0=F+(A+B)*(x/h-y/k)^t G_Z^-1*(x/h-y/k),
    F=1-2A*L-2B*M-(A+B)*C-A*B*(L-M)^2.             (1)

This is the actual corank equation, including the separate supporting
forest Z, not a formal sufficient-completion assertion. In each row
|u(l1-l2)|<=2, so |L-M|<=4. Consequently

    (L-M)^2<=4*(|u(l1-l2)|+|v(m1-m2)|).

Thus F>=1-H_c-H_d, where

    H_c=2A*u*l1+2B*u*l2+(A+B)*c0*u^2
        +4AB*|u(l1-l2)|.

Every term increases with u, so the native upper bounds give

| profile | upper bound for H_c |
| --- | ---: |
| isolated | 2(A+B) |
| end3 on Q | A+3B+4AB |
| end4 on Q | 2B+4AB |
| second-node3 on Q | 4B+8AB |
| both-end3 | A+B |

All five are at most2(A+B), because

    A-B-4AB >= 1/5-1/21-4/105=4/35>0.

It follows, uniformly in Q, both native heights and every spectator,

    F>=1-4(A+B)>=1-4(1/5+1/21)=1/105>0.

The other term in(1) is nonnegative, a contradiction. This proves the
quadrilateral transfer without classifying spectator energies or running
a finite profile table. In particular both-end3 changes are included
in this same face inequality; a now-zero cycle corner is not spent
recursively at a different now-zero cycle corner.

For triangles use the accepted analytic cycle lemma: at order>=10 and
marked distance>=5, deleting either mark permits third-side contacts only
within two steps of that mark. The two neighborhoods are disjoint and
the marks are forbidden by the opposite deletion. The third side would
be orthogonal to the positive affine kernel, contradicting its strictly
positive height at the old ideal corner. Thus no changed triangle occurs.

Give this independent theorem priority whenever an actual order>=26
cycle5 corner occurs. Its output units have distance<=5 in a whole
finite component or distance1 at another corner, so none is changed.
If no such actual corner occurs, every old face argument is unchanged.
This also repairs any earlier route whose output used one of the newly
zero entries; there is no recursive dependency on that earlier route.

## 3. Exact boundary: cycle25 forces a new coupon or new global input

The threshold Q>=21 cannot be replaced by Q>=20 in this mechanism while
keeping all finite weights unchanged. Take a square-two25-cycle F with
a=0,b=5, and add c,d with sole contacts c.a=d.b=-1. All other new
products are zero; there is no spectator. Thus P=5,Q=20,h=k=1.
Let W=F-{a,b}=A4+A19. Its inverse elimination has alpha=1/4 and
four-black Schur matrix, in order a,b,c,d,

    [ 1/4 -1/4 -1  0 ]
    [-1/4  1/4  0 -1 ]
    [-1    0    2  0 ]
    [ 0   -1    0  2 ].

Eliminating the bottom2I leaves the rank-one negative matrix
(-1/4)*ones(2,2). Hence the full27-root Gram has inertia(25,1,1),
the required rank26 and corank1. Each W+one-black edge is strictly
positive definite. Its four corners are exactly

| corner | complete selected-corner diagram and marked roles |
| --- | --- |
| (a,b) | affine cycle25, marked distance5 |
| (b,c) | whole A24 plus A1; b at zero-based A24 position4, c in A1 |
| (c,d) | whole A4+A19+2A1; c,d in the two different A1 factors |
| (d,a) | whole A24 plus A1; reflected position4 against A1 |

Every contact is native3. The exact determinant is zero because
1-4(1/5+1/20)=0. Both separate retained costs are zero, strictly
below2. All finite cross roles in this table have weight zero in the
current manuscript; no same-component or protected-leaf entry can pay
these distinct-component marks. Thus removing only the cycle5 unit
would make the whole quadrilateral total zero, even if ALL present
finite same-component payments were retained.

This is a genuine locally realizable acute Lorentz star. After W
projection let U=a+b,V=c+d. Then U^2=0,V^2=4,U.V=-2.
The timelike point U+V/4 has square-3/4 and strictly negative products
with all four projected sides. The adjacent corner Grams give one
ideal and three finite right-angle corners: an ideal Lambert
quadrilateral. A small correction off W yields a strict ambient chamber
point, preserving these four side inequalities. It is not a singular
Gram mistaken for an actual hyperbolic face.

Orthogonal square-two roots z_j can pad its dimension. At the old cusp
put beta_j=nu-z_j, with nu=sum_(F) f_i. Then z_j+beta_j=nu,
z_j.beta_j=-2, beta_j.c=beta_j.d=-1, and every other added product
is zero. These complete the padded old cusp by orthogonal affine A1
factors. They are strict on the quadrilateral except at its old ideal
corner, so do not change its sides or the displayed finite corners.
The example therefore survives local cusp-completion requirements too.

No finite-volume completion of the entire ambient chamber is asserted.
The counterprofile proves that the existing actual-star, native-height,
separate-edge-budget and local cusp conditions alone cannot yield the
order25 distance-four transfer with unchanged finite cross weights.
A new coupon (for example a payment at A1/A24 position4), or a genuinely
additional global obstruction, is required at this boundary. It does
not disprove an ideal-four theorem with newly financed finite coupons.

## Source boundary

The unchanged manuscript is `manuscripts/finite_volume57.tex`, SHA256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The actual-star/corank, positive kernel, spectator identity, complete
large-path extension and fixed-leaf arguments are those explicitly
proved in the following existing sources under results/:

* `fresh_cusp_factor_payment_framework_20261017.md`;
* `fresh_simple_cusp_affine_a_quadrilateral_20261017.md`, determinant
  and end-moment formulas only, not its all-simple height exclusions;
* `fresh_distance8_ideal_corner_replacement_20261017.md`, the full
  native end-profile list and positive spectator correction;
* `fresh_order10_cycle_half_transfer_20261018.md` and
  `fresh_cycle11_half_pair_transfer_20261018.md`, actual corner/whole
  component and independent-priority machinery;
* `finite57_structural_reentry_targets_20261021.md`, accepted source
  `822953841f6b551e590496b6d617afb67c5f5eceba37ad1fd9566c26d9ef53f7`,
  for the all-order triangle lemma.

Only the uniform determinant bound and the exact cycle25 boundary above
are new here. The proposed full7/4 rule, its finite budget and ceiling53
remain unproved. This source has not yet received independent review.
