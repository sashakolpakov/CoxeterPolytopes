# Two fork/leaf half refunds at whole D19+A1 vertices

19 October2026. **New written consumer theorem candidate, awaiting
independent review.** Paper and source work only. No scientific job,
enumeration, optimizer, build or Lean process was run. All prior
source-bound A/B/D evidence remains unchanged.

## 1. Statement, whole-forest scope and separation from earlier rules

Let P be a compact arithmetic simplest-type Coxeter polytope in
dimension20 over one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

At a vertex whose ENTIRE diagram is D19+A1, let a be the D19 fork,
b,e its two short leaves and h the A1 spectator. Subtract1/2 at
each of the two fork/short-leaf marks (a,b),(a,e), retaining weight
1/2 at each. Leave the short-leaf-pair mark unchanged. The proposed
saving is exactly1 per such whole forest, from36 to35.

The charge is the eight-field distance-at-most-two charge plus its
entire-I2(m>=6) unit and the already reviewed full-A20, full-D20 and
B reductions. The present source's complete white graph is

    A16+2A1.

It differs from every earlier negative support: A17+A1, A18,
B17+A1, B18 and B16+2A1. The present theorem therefore has no
hidden mixed A20/D20 consumer or B-coupon interface. In particular
we do NOT refund the D19 short-leaf pair, whose white graph would
instead be A17+A1 and would require a different mixed proof.

Triangles retain their unchanged entire high-label I2 weight2.
Any quadrilateral with at least two of the present half-weight
sources retains at least1. If it has any other unchanged positive
angle, that angle alone pays1. Hence the sole remaining consumer
is a quad with exactly one selected source and three ORIGINAL
zero marks. This is what the proof below excludes. Larger polygons
retain nonnegative charge. No claim for D+D, D+E or another
spectator forest is made.

## 2. Coordinates and complete spectator contact alternatives

Write the D19 long arm as a-l1-...-l16 and put

    N=17,       s=19/68.

All old bonds have label3. Use unit outward roots. The complete white
support at (a,b) is

    W={e,h,l1,...,l16}=A1+A1+A16.

Let cyclic black sides be a,b,c,d, with

    a.c=-x,  b.d=-y,  b.c=c.d=d.a=0.              (1)

The latter equalities are forced by the three original-zero marks.
Eliminating the two isolates and the long path gives a_bar^2=s,
b_bar^2=1 and a_bar.b_bar=-1/2. The c,d extensions are both finite
at their adjacent actual vertices. Their native contacts to W
therefore obey the complete finite-classification lists below.

At (a,d), zero weight forces d.e=d.l1=0. In the old A18 path
e-a-l1-...-l16, its only remaining contacts are a simple far endpoint,
a far endpoint4, a simple far next-to-end, or none. Including h,
the COMPLETE d possibilities are

* I: white-isolated;
* J: only d.h=-delta, arbitrary finite delta=cos(pi/ell);
* F: d.l16=-1/2, with d.h=0, -1/2 or -sqrt(2)/2;
* B: d.l16=-sqrt(2)/2, with d.h=0;
* T: d.l15=-1/2, with d.h=0.

For F, the further h bond lies at the end of a large path and can
have label3 or4. For B it would move a4-bond into the interior;
for T it would change arms1,1,long into arms1,2,long. Those large
components are not finite. No other long contact or multiple long
contact is possible.

At (b,c), b is white-isolated, so c has the complete one-root
continuation of A16 plus two isolates. It can meet at most one long
node. Against A16 the choices are endpoint3/4, next-to-end3, or none.
If it meets an endpoint by3, it may also meet one isolate by3/4,
or both isolates by3 (a D19 fork). If it has endpoint4 or a
next-to-end3, it cannot also meet an isolate in these large ranks.
If it has no long contact, its contacts to e,h give an ordinary
finite diagram of rank at most3, including every single dihedral
label. Thus arbitrary finite spectator labels have not been discarded.

Finally the original-zero corner (c,d) gives

    c.h*d.h=0.                                     (2)

Otherwise the mark has distance2 through h. The proof splits into
neither spectator contact, c-only contact, and d-only contact.

## 3. Arithmetic interfaces and one precise finite-c residue lemma

We use cyclic integrality, acute all-wall projection and strict
nonphysical positivity of independent principal blocks. These apply
to proper blocks of the actual ambient20 chamber, not only to
complete vertex diagrams. All scalar conjugate identities below
are determinant-ratio or cyclic identities; unit entries are not
assumed to belong to K without a cycle.

For the old D19 block, omit h for a moment. If c has finite nonzero
x, its old D19 projection cost exceeds1: the formula

    q_c=1-t0+2N(x+kappa)^2

has t0<=1,kappa>=0 and x>=1/2. Thus it is an independent
20-root Lorentz block even when c also touches h. Its conjugate
cost is less than1 by ambient arithmetic admissibility. Adding
the nonnegative squared h cost can only strengthen this inequality.
Consequently the PRE-CORANK finite-c reductions in
`arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`
apply with all their actual contacts preserved. They use no d wall
until the expressly later corank/integrality tests.

Here are the complete relevant outcomes in the eight fields:

|old-white c contacts|finite-x outcome before any quad corank test|
|---|---|
|none|squared-contact norm forces rational negative determinant, hence Q|
|one simple long contact, no e|all labels excluded by the cosine-moment test|
|simple endpoint and simple e|near x=cos(pi/4) or far x=cos(pi/5), with actual ground fields sqrt(2),sqrt(5)|
|long endpoint4, no e|non-native quadratic alternatives, or rational negative determinant and Q|
|simple long endpoint and e bond4|the old cyclotomic conjugation contradiction|
|only e|only (x,e-label)=(1/2,6) or (cos(pi/12),4) in the eight fields|

In the last row the m>=4 proof gives K=Q(x^2); labels5 and8
give the excluded exact fields sqrt(5),sqrt(2). Label12 permits
e-label4 or12, and the equal-label12 case already violates
conjugate positivity. For x=1/2 the actual cycle puts the e cosine
in K; its conjugates force only e-label6 (label3 is affine).
Thus no quad identity or desired refund is hidden in this table.
For long4 the rational negative determinant forces K=Q even before
the older source's later scalar-integrality contradiction.

When c.h=0 and x=0, the same source's positive-null and positive-index
tests leave only c contacting e with label12. In particular a far
simple extension produces a positive D20 block inside h-perp,
whose positive index is19, and is impossible. For an e-only contact,
the old cost19*cos^2(pi/ell)/2 gives finite ell=5 or12; field equality
excludes5 in the eight fields. These three possible zero/finite
residues, all with c.h=0, have the SAME scalar

    R0=X^2/t=(7+4sqrt(3))/4,      K=Q(sqrt(3)),       (3)

where t=c_bar^2 and a_bar.c_bar=-X after eliminating W. They are

    (x,e-label)=(0,12),(1/2,6),(cos(pi/12),4).

All have no long contact. Formula(3) follows by substituting their
native cosines, independently of a d wall. Put P0=R0-s>3.

If c.h!=0, the finite residues disappear. The (1/2,6) case has
cos^2(pi/6)+cos^2(pi/3)>=1 at W+c and is nonfinite. The
(cos(pi/12),4) case permits only the further h label3. Its old
conjugate projection cost after that addition is

    (61-34sqrt(3))/2>1,

because sqrt(3)<59/34. Zero x is also impossible: the e-only
candidate12 cannot coexist with any h bond; the far simple long
extension with an h bond3/4 has total cost1/2+gamma^2<=1,
contradicting positive index or compactness. The other old zero
cases retain their principal affine or conjugate obstruction.
An isolated c with only h contact has total cost<1 as well. Hence

    c.h!=0  implies  x>1.                           (4)

For c.h=0 the complete alternative is x>1 or one of the three
residues(3). Native x=1 is always excluded by compactness.

## 4. No spectator contact: the exact inherited local argument

If c.h=d.h=0, every selected root except h lies in h-perp, a Lorentz
space of signature(19,1). The same21-root local corank equations,
actual four finite corners, arithmetic strict determinant inequalities
and compact positive-null exclusions used by the reviewed D19
fork/leaf ORIGINAL-ZERO proofs apply. The complete isolated and
far F/B/T branches were proved in

* `arithmetic_dn_fork_whisker_quad_attachment_classification_20261019.md`;
* `arithmetic_dn_nonisolated_quad_finite_opposite_obstruction_20261019.md`;
* `arithmetic_dn_far_fork_quad_obstruction_20261019.md`;
* `arithmetic_dn_far_double_quad_obstruction_20261019.md`;
* `arithmetic_dn_far_simple_unpaid_quad_obstruction_20261019.md`.

In the far-simple alternative the retained short-leaf unit belongs
to D18+A1 within the old19 support, so the full ambient forest has
an additional h spectator and is not a present D19+A1 fork source.
It remains paid. We import that positive payment alternative, not
a claim that every geometric far-simple quad is impossible.

This is an algebraic local import with h perpendicular to every
root used. We do NOT assert that an arbitrary Coxeter facet is a
Coxeter polytope or transfer an ambient theorem by such an assertion.

## 5. Only c touches h

Let gamma=-c.h>0. Then t=c_bar^2<=3/4 and x>1 by(4).
Write a_bar.c_bar=-X, so X>=x>1.

### 5.1 d is I

Both b,d are unchanged white-isolated unit opposite walls, so y>1.
The exact corank equation is

    (X^2/t-s)(y^2-1)=1/4.

The first factor exceeds4/3-s=215/204. Consequently

    4<Z=4y^2<1064/215,       0<Z'<8/19.

Here Z is integral primitive in K from the old D19+d determinant
ratio1-19Z/8. Degree at least3 has norm below
(1064/215)(8/19)^2<1. Degree2 has trace5 and norm1 or2, giving
fields sqrt(21) or sqrt(17), outside the eight. Degree1 gives Q.
Thus no I wall occurs.

### 5.2 d has far profile F,B or T, with d.h=0

The reviewed old-D19 far-wall finite-product tests depend only on
that old block and d, not on c or h. They therefore still give y>1.
The old d residual square u and cross r=-a_bar.d_bar satisfy

    u<=9/17,       0<r<=1/17.

Joint finiteness W+c+d, including c's h neighbor, forces any long
contact of c to be a simple contact at the NEAR endpoint. A same
far endpoint would give arms1,2,long; a4-bond with a further h
neighbor or a next-to-end contact is likewise nonfinite. Therefore
the white c/d cross sigma is either0 or r. This is the full
joint-continuation statement used in this estimate.

Project off the Lorentz span of a_bar,c_bar. Its Gram is
[s,-X;-X,t], and the remaining b,d Gram is positive semidefinite.
Put D=X^2-st. Then

    y <= sqrt(BV)+Zc,
    B=1+t/(4D),
    V=u+(t*r^2+2X*r*sigma+s*sigma^2)/D,
    Zc=(t*r+X*sigma)/(2D).                         (5)

All positive correction fractions decrease with X>=1. Since
s<3/10,t<=3/4,D>31/40 and sigma<=r,

    B<5/4,   V<9/17+4/17^2=157/289,
    Zc<35/527<1/15.

The product BV is below25/36. Thus (5) yields
y<5/6+1/15=9/10, contradicting y>1. This covers every d when
only c touches the spectator.

## 6. Only d touches h and d has no long contact

Write d.h=-delta, u=1-delta^2<=3/4. Here c.h=0 and d has no other
white contact. With R=X^2/t, P=R-s>0 the exact rank equation is

    P*(y^2/u-1)=1/4.                              (6)

Projected opposition gives y^2>u and y>0. The old D19+d block
has cost19y^2/2>1 even for any finite native y>=1/2. Consequently

    K=Q(y^2),       0<(y^2)'<2/19.                (7)

For finite y the usual exact cosine-square trace/norm step leaves
labels5 or12. Label5 gives the excluded exact ground field sqrt(5).
Thus the sole finite possibility in the eight fields is

    y=cos(pi/12),       K=Q(sqrt(3)).

Since delta^2 belongs to this K, its native finite labels can only
be3,4,6 or12. The full old D19+A1+d cost is19y^2/2+delta^2.
Its other conjugate is less than1 only for delta-label3 or12;
labels4,6 already exceed1. These two cases are treated below.

If y>1, consider c's alternatives from Section3. For a residue(3),
P=P0>3 and (6) gives y^2<(3/4)(1+1/12)=13/16<1. For x>1,
P>1-s=49/68 gives

    4<Z=4y^2<198/49,       0<Z'<8/19.

Z is integral primitive by(7). Degree at least3 has norm below1.
In degree2 its trace is strictly between4 and5, impossible. Degree1
is Q. Therefore y cannot be dotted either; native y=1 is forbidden
by compactness. Only the two finite cases remain.

### 6.1 Spectator bond12

Put eta=sqrt(3). Substitution into(6) gives

    R=21/136+eta/12<1.

If x>1, X>=x and t<=1 give R=X^2/t>1. A zero/finite residue
instead has R=R0>3. Both contradict the required value. Thus this
branch is excluded.

### 6.2 Spectator bond3: the exact quadratic scalar

Now substitution into(6) gives

    R=R1:=89/136+(3/8)eta.                        (8)

This is strictly between1 and2. It differs from R0>3, so c must
have x>1. We must nevertheless exclude EVERY finite white
continuation of that dotted c, rather than assuming an isolated c.
The ground field is exactly Q(eta), whose integer ring is Z[eta].

**No long contact and no e contact.** Then x^2=R1. The required
integral two-cycle4x^2 has trace89/17, not an integer.

**Only an e contact beta.** Native beta^2 in K permits precisely
1/4,1/2,3/4,(2+eta)/4. The two-cycles4x^2 and the actual a-e-c
triangle cycle4x*beta are integral. Since

    (x+beta/2)^2=(1-beta^2)R1,

their sum has to equal4(1-beta^2)R1-beta^2. Its field traces at
the four listed values are respectively

    233/68, 55/34, -13/68, -43/68.

None is an integer. Thus every e-only dihedral label is excluded.

**A single simple long contact, no e contact.** There are four
positions. The native cycle makes A=2x=m+ell*eta integral, and
(A+2kappa)^2=4tR1. Compare the eta coefficient to get

    m+2kappa=3t/(4ell).

Here ell!=0. If ell<0 the left side is negative and m<0, making
A<0, impossible. Thus ell>0. Since0<2kappa<1 and the right side
is at most3t/4<1, m=0. The resulting required ell values are

|long position|t|kappa|required ell|
|---|---:|---:|---:|
|l1|9/17|8/17|27/64|
|l2|2/17|15/34|1/10|
|l15|2/17|1/17|3/4|
|l16|9/17|1/34|27/4|

None is integral.

**A single long4, no e contact.** Its native cycle makes
A=2sqrt(2)x=m+ell*eta integral. Put kappa0=16/17 near and1/17
far. The equation is(A+2kappa0)^2=(8/17)R1, so

    m+2kappa0=3/(34ell).

Again ell>0 by positivity of A. The smallest positive left side
is15/17 near or2/17 far, both larger than3/34. No solution exists.

**Simple long endpoint plus simple e.** Here t=19/68 and kappa is
49/68 near or19/68 far. The same eta comparison with integral
A=2x gives

    0<m+2kappa=3t/(4ell)<=57/272.

Its smallest positive left value is15/34 near or19/34 far, both
larger. There is no solution.

**Simple long endpoint plus e bond4.** The cycle through the old
path and e puts sqrt(2) in K, impossible for K=Q(sqrt(3)).

Finite classification allows no long4 together with e, and no
next-to-end contact together with e. The list is therefore complete.
It excludes the last finite-y branch, including every long/e
attachment of c. This is the arithmetic step that an isolated-c
argument would miss.

## 7. Only d touches h and d also has a far long contact

The d inventory allows only a simple l16 contact, with spectator
bond3 or4. Put r=1/34 and u=9/17-delta^2, hence u<=19/68.
The old full vertex cost is

    q_d=19y^2/2+y+1/2+delta^2.

The long cycle makes2y integral and K=Q(y). At other embeddings,
its absolute value is less than9/20 in either spectator case:
for bond3 use(2+sqrt(42))/19<9/20; for bond4 its conjugates
lie in(-4/19,0). Finite2y<2 gives impossible norm below9/10 in
nontrivial degree. Zero y has q_d<1 or q_d=1, contradicting
positive index or compactness. Therefore y>1.

The original-zero (c,d) mark forces c.l16=0. Joint finiteness of
the large path ending in d-h then allows c only near endpoint3/4,
near next-to-end3, or e-only/isolated old-white contacts, with the
usual finite e additions. In particular, its white c/d cross obeys

    0<=sigma<=2r.

This retains the near next-to-end case; it is not a claim that
every unspecified contact is zero.

If x>1, use the positive-complement inequality(5), now with
t<=1, s<3/10, X>=x>1, D>7/10, sigma<=2r. It gives

    B<7/5,
    V<19/68+(62/7)/34^2<3/10,
    Zc<15/238<1/15.

Thus y<sqrt(21/50)+1/15<2/3+1/15<1, contradicting y>1.

If c is one of the zero/finite residues(3), it is e-only and
sigma=0. In that case P=P0>3 and the same projection gives

    B<13/12,       V<u+r^2/3<3/10,
    Zc<r/6=1/204.

Again y<3/5+1/204<1. Thus all d-long-plus-spectator cases are
excluded as well.

## 8. Conclusion and remaining scope

Equation(2) exhausts the spectator possibilities. Neither contact
uses the exact inherited19-dimensional local algebra; c-only
contact is excluded in Section5; d-only contacts, including all
long additions, are excluded in Sections6-7. Every single selected
fork/leaf source with three original-zero corners has been covered.
Two or more half-weight sources retain at least1 automatically.

This completes the written consumer proof for the TWO fork/leaf
half deletions at ENTIRE D19+A1 vertices, subject to independent
review. It saves1, so that forest's cost becomes35. The distinct
short-leaf pair is unmodified. No spectator generalization to D+D,
D+E or D18+I2, no global allowance35, no compact arithmetic ceiling19,
and no executable or formal verification is claimed here.
