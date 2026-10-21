# A uniform distance-five replacement in long affine A factors

21 October 2026. Root manual proof candidate, pending a full other-author
mathematical review. This is supporting research in the existing
finite-volume paper family. No accepted weight, manuscript, certificate,
dimension bound or execution record is changed. No scientific, symbolic,
enumeration, Lean or TeX process was run.

## Statement

In the current corner-payment rule, set to zero every distance-five pair
in a whole varying affine A cycle of order m>=26. Keep every other weight
unchanged. Every actual triangle still receives at least two and every
actual quadrilateral at least one. No new finite coupon is required.

Thus the intrinsic charge of such an A factor decreases from5m to4m.
This is NOT an improved all-family envelope or global dimension bound:
the B/C/D factors, finite forests and smaller-order cusps still have their
separate obligations. In particular it does not close the ambient57
finite372/curvature364 obstruction. The six delivered papers remain frozen.

The explicit inputs used below are:

| Source under results/ unless stated otherwise | SHA256 | Scope used |
|---|---|---|
| fresh_cusp_factor_payment_framework_20261017.md | c65cfc5a4409a8cc83fcf0a083a6c1a7a29bb7bc7d649c11435773877e70ccab | Actual quadrilateral support, whole varying affine factor, positive-definite spectators, strictly positive opposite heights and the norm-two Lorentz quotient identity |
| fresh_cycle12_ideal5_transfer_20261017.md | 53e69a52b1a72f0b62e37744ea315dab1642404f8c5b52fee9fd6df73877902f | Complete one-root extension profiles of a long deleted A path; fixed simple leaf/degree-two-neighbor recognition of a finite corner |
| fresh_cycle11_half_pair_transfer_20261018.md | ec117a499ddad7cf8e87426e467d54cde295c0956428b2a8bbda5140abcdecd5 | Actual whole-corner interpretation of those profiles, and the short-contact alternatives; no order-eleven scalar bound is transported |
| ../manuscripts/finite_volume57.tex | c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0 | Current simultaneous weights: all marked distances<=2 retain units; finite components of rank>=26 retain distance<=5 units; affine A order>=26 has exactly the distance<=5 units before this proposed change |

These are mathematical interfaces, not claims that hashes prove geometry.
Finite/affine Coxeter classifications and the actual factor-fiber framework
remain explicit inherited inputs. The new estimates below are analytic;
no unexecuted finite enumeration is assumed.

## 1. An electrical gap inequality

Let C be the unit-edge cycle on m nodes. Its Laplacian is the norm-two
affine A Gram matrix. Choose Euclidean point potentials v_i with

    ||v_i-v_j||^2 = d(m-d)/m,  d=(i-j) mod m.

For probability vectors mu and eta supported in disjoint connected cyclic
intervals A and B, suppose the two complementary paths between the
intervals have edge lengths L and M. Then

    ||sum_i mu_i v_i - sum_i eta_i v_i||^2 >= LM/(L+M).       (1)

Here a connected interval can include nodes receiving zero probability.
Consequently enlarging a support to such an interval is legitimate.

To prove (1), put R=LM/(L+M). Define a function f on the cycle to be zero
on A, R on B, and linear along each complementary path. Its Dirichlet
energy is

    sum_edges (f_i-f_j)^2 = R^2(1/L+1/M)=R.

For s=mu-eta, sum s_i=0 and s.f=-R. The squared norm on the left of (1)
is s^t C^+s. To see this identification, ground one point and use the
inverse of the grounded cycle Laplacian; its diagonal/cross polarization
is exactly the displayed resistance formula. Cauchy in the Laplacian
energy inner product gives

    R^2=(s.f)^2 <= (s^t C^+s)(f^t C f)=R(s^t C^+s).

Since R>0, division proves (1). This proves the inequality for arbitrary
continuous contact magnitudes and all cycle orders, not just vertices of
the two convex hulls.

## 2. Actual roots and the complete profile reduction

Write the quadrilateral star as W+a+b+c+d, with opposite pairs a,c and
b,d. At the selected ideal corner, its whole varying affine factor is
the cycle F on0,...,m-1, with a=0 and b=5. Its two arcs have lengths5
and m-5. Write Z for the positive-definite spectator diagram, orthogonal
to all of F. All roots have norm squared two. Native nonzero contacts
have magnitude at least one and every distinct pair has nonpositive
product. These statements hold for actual roots, with no arithmeticity.

An adjacent black contact already gives a retained distance-one unit at
its corner. A common contact of c and d to a node of Z gives a retained
distance-two unit at their corner. We can therefore suppose all adjacent
black products vanish and the two Z-contact supports are disjoint.

Delete a from F. The remaining path has m-1>=25 nodes. Its inverse
diagonals allow c to contact only its two endpoints, with magnitude one
or sqrt2, or their next nodes, with magnitude one. Two contacts are
possible only at both endpoints, both with magnitude one; that gives
the affine-cycle extension. Three contacts are impossible. The same
classification applies to d after deleting b. These are complete
principal-extension statements: adding Z cannot repair a failed principal.

For clarity the inverse diagonal at path node i is i(m-i)/m.
For3<=i<=m-3 this is at least3(m-3)/m>2. Put phi=(1+sqrt5)/2. At an endpoint, magnitude phi
already costs phi^2(m-1)/m>2, while at a second node sqrt2 costs
4(m-2)/m>2. The inverse has positive entries. Both unit endpoint
contacts cost exactly two; every other two-contact combination costs
more than two. Three unit contacts cost at least3(m-1)/m>2.

If c has its single contact at1 or2, the adjacent corner (b,c) has a
marked path of length at most five. This corner is finite: the opposite
old endpoint m-1 is a fixed simple leaf with a fixed degree-two neighbor.
Neither it nor that neighbor can acquire Z contacts. No classical affine
diagram has this leaf configuration, and exceptional affine diagrams
cannot contain the displayed A_(m-1). Thus the whole marked component
is finite, of rank at least26, and its short pair retains a full unit.
Every other component is a positive-definite principal part of Z.
The reflected assertion holds if d contacts3 or4.

This step deliberately does NOT treat a single contact on the long arc
as a short marked path. Such contacts can give very long finite paths
and must instead be kept in the next estimate.

After removing those paid cases, classify c and d independently:

* Type O: no retained contact, or a single contact on the outward side.
  Its full F-contact support is contained in {m-2,m-1,0} for c and
  {5,6,7} for d. This includes any contact to its opposite old mark.
* Type T: both endpoint contacts, each of magnitude one. Its support is
  contained in {m-1,0,1} for c and {4,5,6} for d. Its F-height is at
  least two. No restriction on the opposite-mark magnitude is imposed.

Type O has height at least one. More precisely, its height is either
one, with a SINGLE unit F-contact, or at least sqrt2. In fact a nonzero
native product has magnitude one or at least sqrt2, and the height is
the sum of all F-contact magnitudes. Two nonzero contacts already sum
to at least two. Positivity excludes an empty row. This native gap is
used below; the proof does not replace it by a continuous height range.
The two types exhaust all unpaid profiles, including dotted products.

## 3. The three continuous domains

Let nu be the null vector sum_i alpha_i of F, h=-(c,nu)>0 and
k=-(d,nu)>0. Dividing each full F-contact row by its height gives
probability vectors mu and eta with exactly the support containments
above. The Euclidean quotient at nu splits into the F-space and the
orthogonal spectator space. The norm-two Lorentz identity and (c,d)=0
therefore give

    ||sum mu_i v_i - sum eta_i v_i||^2
      + spectator_difference_energy = 2/h^2+2/k^2.          (2)

The spectator term is nonnegative. We use it only in that direction;
neither row's energy is substituted for a bound on their difference.

In the O/O case, first suppose one height is at least sqrt2. The cyclic
support intervals have complementary gaps5 and m-9. Equation (1) bounds
the left side of (2) below by

    5(m-9)/(m-4) >=85/22>3, because m>=26.

But the right side is at most2+2/2=3. If neither height is at least
sqrt2, the native gap forces h=k=1 and each full F-contact row consists
of a single unit. The contact nodes are then one of {m-2,m-1,0} and
one of {5,6,7}. Their separation through the short gap is an integer
d in {5,6,7,8,9}. Since m>=26, this lies below m/2 and the resistance
d(m-d)/m increases with d. Its minimum is

    5(m-5)/m >4, because m>25.

The right side of (2) is exactly four. Both O/O alternatives are
impossible, without using a finite cutoff on opposite products.

In the T/O case the gaps are4 and m-8. Thus the lower bound is

    4(m-8)/(m-4) > 5/2.

Indeed its difference from5/2 is (3m-44)/(2(m-4))>0 for m>=26. The right side of
(2) is at most2/4+2=5/2, since the T-height is at least two. Reflection
treats O/T with the same estimates.

Finally T/T has gaps3 and m-7, yielding

    3(m-7)/(m-4) > 1,

whose difference from one is (2m-17)/(m-4)>0. Both heights are at least
two, so (2)'s right side is at most one. Contradiction.

All unpaid profiles are impossible. Every affected actual quadrilateral
therefore has one of the explicitly retained units in Section2, with no
new finite coupon or dependence on another deleted distance-five entry.

## 4. Triangles and simultaneous use

In an actual triangle, its third side must give positive-semidefinite
extensions after deleting a and after deleting b. The first allows
contacts only at {1,2,m-2,m-1}; the second only at {3,4,6,7}. These sets
are disjoint. Contacts to a or b themselves are forbidden by the other
deletion, since 5(m-5)/m>2. Hence the third root is orthogonal to F.
Its height at nu would vanish, contradicting the strictly positive
height of a side not containing this ideal vertex. No triangle is
affected, irrespective of spectators.

Rotation and reflection cover every distance-five pair. If an actual
quadrilateral has any changed corner, apply its direct retained-unit
proof above; that unit has distance at most two, or lies in a finite
component of rank>=26 at distance at most five, so it is unchanged by
ALL these deletions. With no changed corner the old face proof is
unchanged. This establishes simultaneous compatibility with the current
rule, rather than a recursive promise to pay a newly deleted corner.

For m>=26 there are exactly m unordered distance-five pairs, so the
proposed A-factor total becomes4m. Other affine families still control
the all-family envelope. A full other-author review, plus any desired
exact corroborating check, is required before changing the accepted
weight tables. No lower global ceiling is claimed by this candidate.

## 5. Why order25 is a genuine boundary of the local argument

This observation concerns local face data, not a finite-volume example.
At m=25 take the two supporting paths A4 and A19, with the old cycle
marks a,b at distance five. Write alpha=1/5+1/20=1/4. Give c only its
unit contact with a, and d only its unit contact with b; no spectator
contacts are present. Every adjacent black product is zero. The original
corner is an affine25-cycle. The two adjacent corners are whole A1+A24,
with the A24 mark at zero-based index4 or19; the opposite corner is
2A1+A4+A19, marked at the two A1 factors.

These are genuine compatible local Lorentz data. After eliminating the
two supporting paths, in coordinates with square u^2+v^2-t^2 the four
residual normals can be taken to be

    a=(1/2,0,0),        b=(-1/2,1,1),
    c=(-2,1/2,3/2),     d=(0,-3/2,-1/2).

Their squares are alpha,alpha,2,2; their adjacent products are
-(alpha),0,0,0 and their opposite products are -1,-1. Restoring the
old a,b positive projections gives precisely the native norm-two root
Gram just described. The strict face center (-1,1,2) has square -2
and product -1/2 with each residual normal. Its consecutive vertices
can be represented by

    (0,1,1), (-4,5,7), (-2,1,3), (0,1,3).

The first is ideal and the other three are timelike; the products
against the unselected sides are strictly negative. Thus the section
is an actual one-ideal, three-right-angle quadrilateral of this local
acute chamber. These formulas are hand identities, not a claimed
executed matrix certificate.

Under the current finite rule its three finite corners have zero
payments at the displayed marks: there is no A1/A24 index4-or19 coupon
and no payment between these separate A1 factors. Removing its ideal
distance-five unit would therefore leave this LOCAL face unpaid.
This does not assert that the local chamber admits a finite-volume
completion. It shows exactly why a proof extending the deletion to
order25 must add a compensating coupon or use more than the present
supporting-star geometry. The strict order26 threshold above avoids
that boundary, rather than silently discarding it.
