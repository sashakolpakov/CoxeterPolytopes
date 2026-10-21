# An essential E7(2,5) quadrilateral is impossible at the compact29 target

19 October 2026. **New written candidate, awaiting independent review.**
No scientific run, finite census, optimizer, Lean or TeX build ran.
Earlier reviewed sources and the primary manuscripts are unchanged.

This note proposes a genuine allowance improvement for
T=E7+E7+E7+E8: the ALREADY valid facewise refund has charge at most84
at every T vertex. It uses no arithmeticity or field-of-definition
assumption. If accepted, it removes this one type from the98-type
frontier;79 types at87 and18 at85 remain, so the compact dimension
ceiling is still29.

## 1. Statement and inherited interfaces

Number each E7 by 0--1--2--3--4--5 and 2--6, with norm-two roots.
Use the original unconditional Vinberg angle weights w0. A connected
marked pair in rank at most seven has positive weight; an E7 angle
has weight one; a connected pair in rank at most four has weight two.
The old face inequality is sum w0>=5-k. Finite Coxeter classification,
compact acute simplicity, actual-face projection and exclusion of
positive-null principal diagrams are the explicit inherited inputs.

**Essential-quad lemma.** At a T vertex the two-face at marked E7
pair(2,5) cannot be a quadrilateral whose other three w0 payments
are zero.

The previously reviewed packet
`compact29_e7_triple_facewise_refund_20261018.md` proves:

* any designated triangle at(2,4),(2,5),(3,5) forces six actual
  triangles with retained weight-two source corners, giving charge81;
* an essential quad at(3,5) forces an actual triangle at(4,5) with
  payments1,0,2; its old unit is fully refunded;
* finite E7-to-E8 extension has projection square3/2, hence reduced
  norm1/2 for a norm-two new root;
* the valid global facewise refund and its protected-source
  composition with the original-heavy rank4 rule are as recalled
  in Section7 below.

The new content is a complete own proof of the essential-quad lemma,
including arbitrary opposite labels/distances, all spectator contacts,
and an unbounded rank-two label branch treated analytically.

## 2. Zero-payment corners and all possible contacts

Write a=r_2, b=r_5 and cyclic side order a,b,c,d. The five local
white roots split as

    L={r_0,r_1}=A2,   R={r_3,r_4}=A2,   h=r_6=A1.

Root a contacts r_1,r_3,h, while b contacts only r_4. The full white
support also contains the untouched2E7+E8. Every white-plus-adjacent-
black-root corner is positive definite.

A zero marked pair at any new corner must be disconnected. Indeed,
the whole old E8 cannot be extended finitely. A new root contacting
a whole E7 makes E8 and admits no further root in that component.
Otherwise a component containing both marks is in the local diagram
of rank at most seven, where its pair has positive weight. Thus an
exceptional omitted E8 weight cannot evade this disconnection claim.

Consequently

    c.b=d.a=c.d=0.

At the b,c corner, b is already connected to every root of R, so
c is orthogonal to R. At the d,a corner, a is connected to ALL
local white roots, so d is orthogonal to L,R,h. Apart from a possible
contact with one untouched E7, c can contact only L+h. A spectator
E7 contact is necessarily its one-root E8 extension, with no other
contacts in that corner. Thus c then has no local white contact;
d likewise has no local contact. The two new roots cannot use the
same outside E7, since the opposite corner contains both of them.
All untouched outside factors are therefore either orthogonal or
explicitly included in the reduced norms below.

## 3. Exact Schur equation and native windows

Let eta denote c's nonnegative local contact magnitudes. Put

    C=2-eta^T(L+h)^(-1)eta,
    s=(L^(-1)eta_L)_r1 + eta_h/2.

If c instead uses an outside E7, C=1/2 and s=0. If it is
orthogonal to all white roots, C=2,s=0. For d put D=2 when it
uses no spectator and D=1/2 when it uses a spectator E7.

Write x=-a.c, y=-b.d and Y=y². These are actual native products:
zero, 2cos(pi/m), parallel2 or strictly dotted>2. Eliminating all
white roots gives

    [[1/6, -1/3, -(x+s), 0],
     [-1/3, 4/3, 0, -y],
     [-(x+s), 0, C, 0],
     [0, -y, 0, D]].

For example a's white projection square is2/3+2/3+1/2=11/6,
b's is2/3, and their cross projection is1/3. Actual opposite
sides are strictly disjoint in the intrinsic compact quadrilateral,
and hence their Schur pairs are strictly indefinite. The four-side
Schur matrix has rank three. Eliminating c,d therefore yields

    (q-1/6)*(Y/D-4/3)=1/9,
    q=(x+s)²/C>1/6,     Y/D>4/3.                 (1)

For D=2, Y>8/3. The smallest native squared bond above8/3 is3
(label6): label5 has value phi²=(3+sqrt(5))/2<8/3. Hence Y>=3,
and (1) gives

    1/6<q<=5/6.                                (2)

For D=1/2, Y>2/3 implies Y>=1, so

    1/6<q<=1/3.                                (3)

Both windows allow finite, parallel and arbitrarily dotted opposite
pairs initially. They are consequences, not numerical cutoffs.

## 4. Complete finite c-profile domain

The inverse L matrix is (1/3)[[2,1],[1,2]], and h's inverse entry
is1/2. Every contact here occurs in the finite b,c corner. The
connected finite extension of L+h by c has the following COMPLETE
possibilities:

* no contact;
* contact at one L endpoint with label3,4 or5, giving A3,B3,H3;
* contact only at h with ANY finite label m>=3, giving I2(m);
* contact one L endpoint and h, giving a path with successive labels
  3,m,n, where (m,n)=(3,3),(3,4),(4,3),(3,5), namely
  A4,B4,F4,H4;
* a spectator E7 extension, with no local white contacts.

Two L contacts would create a finite-diagram cycle. The four-node
connected cases are exactly the finite path classification above;
a fork cannot arise from adding a vertex to an edge plus an isolate.
There is no unproved high-label cutoff on the isolated h branch.

If x>0, native x>=1. Every NONZERO local profile has

    (1+s)²/C>=4/3>5/6.

For an L-only contact, s>=1/3,C<=4/3; for an h-only contact,
s>=1/2,C<=3/2; for joint contacts, s>=5/6,C<=5/6. These prove
the stated lower bound, including all larger labels. A spectator
profile gives x²/C>=2. Thus all nontrivial profiles require x=0.

The no-contact profile C=2,s=0 has x>0. Window(2) requires
1<=x²<=5/3, hence x²=1, since the next native squared value is2.
Equation(1) then gives Y=10/3. Under window(3) it has no solution.
The value10/3 is not a native finite squared product: it is below4,
and a rational finite value4cos²(pi/m)=2+2cos(2pi/m) is an algebraic
integer and hence an integer. This elementary root-of-unity fact
holds for EVERY Coxeter angle; it is not an arithmetic-lattice or
field-of-definition assumption.

For every other profile set x=0 and use q=s²/C. Write phi=(1+sqrt5)/2
and z=phi². The exact filter is:

|Profile|C|s|q=s²/C|
|---|---:|---:|---:|
|L label3 at r0|4/3|1/3|1/12|
|L label3 at r1|4/3|2/3|1/3|
|L label4 at r0|2/3|sqrt2/3|1/3|
|L label4 at r1|2/3|2sqrt2/3|4/3|
|L label5 at r0|(6-2z)/3|phi/3|z/(18-6z)>5/6|
|L label5 at r1|(6-2z)/3|2phi/3|four times previous|
|h label3|3/2|1/2|1/6|
|h label4|1|sqrt2/2|1/2|
|h label m>=5|(4-z_m)/2|sqrt(z_m)/2|z_m/[2(4-z_m)]>5/6|
|joint(3,3), r0|5/6|5/6|5/6|
|joint(3,3), r1|5/6|7/6|49/30|
|joint(3,4), r0 or r1|1/3|1/3+sqrt2/2 or2/3+sqrt2/2|>=11/6+sqrt2|
|joint(4,3), r0 or r1|1/6|sqrt2/3+1/2 or2sqrt2/3+1/2|>=17/6+2sqrt2|
|joint(3,5), r0 or r1|(8-3z)/6|1/3+phi/2 or2/3+phi/2|>6|
|spectator E7|1/2|0|0|

Here z_m=4cos²(pi/m). The h-label m>=5 inequality is equivalent
to z_m>5/2, true already at m=5. The H3 inequality follows from
z>5/2; the H4 bound follows from C<1/6 and s>1. These strict
bounds treat all larger m at once. Values q<=1/6 fail strict
opposite indefiniteness; values above5/6 fail both windows.

## 5. Every candidate left by the windows is excluded

For D=2, the two L-only rows with q=1/3 force Y=4. The native
pair b,d would be parallel, forbidden in a compact chamber. The
h-label4 row with q=1/2 forces Y=10/3, excluded as above. The
joint(3,3) at r0 has q=5/6 and forces Y=3. But its actual roots

    a--r1--r0--c--h--a

form a simple five-cycle with all edge products -1 and no chord:
x=0, c contacts only r0,h. The all-ones vector is positive null.
Compactness excludes this principal diagram. This also excludes the
last nontrivial row; the earlier no-contact case was already removed.

For D=1/2, only the two L-only q=1/3 rows remain. Both force Y=1.
They each have an explicit compactness obstruction, independently
of the spectator used by d:

* Label3 at r1: the roots {r0,c,r1,a,r3,h} form two forks r1,a,
  joined by one edge, with respective leaves r0,c and r3,h.
  Give each fork coefficient2 and the four leaves coefficient1.
  Every row of the norm-two Gram matrix annihilates this positive
  vector. This is the usual affine D5 diagram.
* Label4 at r0: the roots {c,r0,r1,a,r3,h} have path
  c--r0--r1--a with first label4 and remaining labels3, and
  two simple leaves r3,h at a. Give coefficients sqrt2,2,2,2
  along that path, and1,1 on the leaves. Again every Gram row
  is zero. This is a positive-null diagram (usual affine B5).

All listed roots are actual walls. All unlisted pair products in
these principal sets are zero by the full profile and x=0. A
strict interior point pairs negatively with the positive combination,
so its null vector is nonzero and future; acuteness retains every
other halfspace. It would be a feasible ideal point, contrary to
compactness. This is the same full-wall positive-null proof used
in the earlier packets, not an assertion that every subdiagram is
automatically a face. All candidates are excluded, proving the lemma.

## 6. Spectators, labels and consumers that were not omitted

The proof retains both untouched E7 components as possible contacts
for either new root, and the entire old E8 as a sealed factor. Both
spectator residual norms1/2 and2 were treated. It does not assume
that c,d are isolated from outside factors. It handles c's isolated
I2(m) contacts at every label, and handles arbitrary native opposite
x,y before the exact Schur windows. The only integrality used is
the elementary cyclotomic property of a finite angle itself, valid
in arbitrary Coxeter polytopes.

## 7. A full three-unit saving at T for the existing facewise rule

Use the original-heavy rank4 rule as background w1. Its zero set
is the same as w0, every positive value is an integer at least one,
and it changes no T angle. Its protected weight-two triangle sources
have original whole charge80,84 or78 and are untouched. Let eligible
angles be every E7 angle at a T vertex. For each actual k-gon F, let
m count eligible corners and B sum all other background weights.
The established rule replaces each eligible unit by

    q(F)=max(0,(5-k-B)/m),

leaving all other weights fixed. The old inequality gives0<=q<=1,
and the new face sum B+max(0,5-k-B)>=5-k verifies EVERY simultaneous
consumer. No new deletion is made on the assumption that another
changed corner is fixed.

If any designated triangle occurs, the reviewed active argument
still gives charge<=81. Otherwise, in each E7 factor consider the
two distinct marked angles(2,5) and(3,5). Both faces have at least
four sides. The new lemma says(2,5) cannot be an essential quad.
Thus it refunds at least1/2: a long face refunds one; a quad with
B>=1 refunds one; a quad with B=0 has m>=2 and refunds at least1/2.

For(3,5), a long face or a nonessential quad similarly refunds at
least1/2. An essential quad instead forces the DIFFERENT actual
(4,5) triangle with an untouched weight-two source, fully refunding
that old unit. Thus in every factor the combined saving is at least
one: two half refunds, or the(2,5) half refund plus the distinct
(4,5) full refund. Across the three E7 factors all counted angles
are distinct. Consequently the total saving is at least three,
and EVERY T vertex has charge at most87-3=84.

This is a stronger estimate on the already valid rule, compatible
with all rank4 changes. It supplies no saving at the other79+18
remaining heavy types. No dimension ceiling or construction claim
follows before those consumers and the full vertex domain are closed.
