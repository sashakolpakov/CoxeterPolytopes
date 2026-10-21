# Both A19/A20 end-distance-two units: a complete consumer proof

19 October 2026. **New written theorem candidate, awaiting independent
review.** No scientific job, proof script, optimizer, Lean process or build
was run. Previously stable evidence and weight sources are unchanged.

The assertion is restricted to whole A19/A20 vertices in the eight
high-label fields. If the written proof passes review, both simple-end
distance-two units can be deleted, reducing a full A20 cost from 37
to35. This does not establish a global allowance35: the B quarter-gap,
D20 and outside-factor cases remain separate obligations.

## 1. Exact rule, field scope and inherited inputs

Let P be compact arithmetic simplest-type of dimension n=19 or20,
with ground field one of

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Start with the distance-at-most-two baseline plus one additional
unit on an entire marked I2(m) component, m>=6, as in
`results/fresh_arithmetic_high_label_20261016.md`. Delete the unit
at the first and third roots at BOTH ends of every ENTIRE A_n vertex.
There are exactly two different such marks at n=19,20.

The inherited inputs are compact Coxeter simplicity, finite diagram
classification, acute all-wall positive-definite projection, absence
of affine principal diagrams in a compact chamber, arithmetic cyclic
integrality and nonphysical positivity. Each specified ground field
has degree at most five. The old eight-field theorem pays two on
every triangle at an entire high-label I2 corner and at least one
on every quad.

The proposed deletion leaves that I2 triangle corner unchanged.
Weights remain nonnegative. Consequently it is enough to exclude a
quad whose total becomes zero, including simultaneous deletions
of this same role at its other corners.

## 2. Complete necessary continuation list for an unpaid quad

At an old selected A_n corner write its path

    a-b-c-t1-...-tL,       L=n-3,

and put the quad sides in cyclic order a,c,E,D. The white support is

    W={b}+{t1,...,tL}=A1+A_(n-3).

An unpaid modified corner has adjacent black product zero: a nonright
adjacency retains its baseline, whereas a deleted A_n mark is an
orthogonal distance-two pair. Thus c.E=E.D=D.a=0.

The whole finite corner W+c+E is a one-root continuation of the long
A_(n-1) path b-c-t1-...-tL. Large finite classification and zero
modified charge leave exactly

    E=I: isolated from W+c;
    E=L: one simple contact at b;
    E=F: one simple contact at tL;
    E=B: one 4-bond contact at tL;
    E=T: one simple contact at t_(L-1).

In L, the marked corner c,E is another deleted A_n end mark. A
4-bond at b or a simple contact at c has an undeleted positive
baseline, so those extensions cannot occur in an unpaid quad.

At W+a+D, either the marked a,D pair had original baseline zero,
or it is another selected A_n end mark. In the original-zero case
D is perpendicular to a,b, and its one-root finite continuation
of the long white path is exactly one of

    isolated;
    simple or 4-bond at t1;
    simple at t2;
    simple or 4-bond at tL;
    simple at t_(L-1).                                (1)

In the selected case it is instead a simple bridge from b to t1
or from b to tL. These are J_near,J_far.

These lists retain EVERY endpoint, 4-bond and next-to-end contact
allowed by finite classification. They do not assume that the new
walls miss the long arm. Exceptional finite types are too small;
a finite one-root extension of the large A path is A, B, D or an
isolated node. Entire corners contain exactly W and two black roots,
with no suppressed outside factor.

The complete J_near/J_far case is proved in
`results/arithmetic_an_same_endpoint_quad_reduction_20261019.md`,
SHA256 `9d45659a751fed16d5ae4558d46c0238703f11cbb0fc65904ce1aa4158e8681f`.
That same-endpoint reduction includes all five E types, both long-arm
orientations and all simultaneous marks. We use it as a separate
written lemma, not as a full-refund assumption. It eliminates J_near
and J_far here. The rest of this note treats EVERY row of (1).

## 3. Fixed projected parameters

Set N=n-2 in {17,18}, epsilon=1/(2N), and

    r=1/4+epsilon.

After eliminating W the old norms are a_bar^2=3/4, c_bar^2=r,
and a_bar.c_bar=-1/4. Put a.E=-x and c.D=-y, x,y>=0.
All original-zero D rows have no contact with a,b.

Write t=E_bar^2 and -h=c_bar.E_bar for a far E row. Direct white
A-path inversion gives

|E type|t|h|
|---|---|---|
|F|(N+1)/(2N)|1/(2N)|
|B|1/N|sqrt(2)/(2N)|
|T|2/N|1/N|

In all these rows t<2/3 and h<1/16. The a,E projected product is
-x. For I, E_bar^2=1 and c_bar.E_bar=0; for L, E_bar^2=3/4,
a_bar.E_bar=-x-1/4, c_bar.E_bar=-1/4.

The full old A_n plus E has a negative physical determinant: the
old n roots span a positive space, and E is an actual additional
wall not containing their finite vertex. Its projection cost is
strictly greater than one. At every nonphysical ground embedding
that same independent block is positive definite, so its cost is
strictly less than one. A ground-field embedding fixing the displayed
cost parameter would preserve its negative physical determinant;
this proves the field-generation assertions used below.

## 4. Far E roots force x>1 in these eight fields

This is where arithmetic rules out finite ambient angles of far
continuations. The actual value x=1 is excluded by compactness:
the positive null sum of the two unit normals would be feasible.
We show that finite cosines x<1 are impossible as well.

### E=F

The long cycle through a, the old A_n path, and E has cyclic product
2x up to sign. Its independent projection cost is

    q_E=[2n*x^2+2x+n/2]/(n+1),
    q_E-1=[2n/(n+1)](x-1/2)(x+1/2+1/n).

Thus K=Q(x), the physical x>1/2, and every other conjugate satisfies

    -1/2-1/n < x' < 1/2.                            (2)

For finite x=cos(pi/m), degree at most five lists
m=4,5,6,7,8,9,10,11,12,15. Label 5 gives the excluded ground field
Q(sqrt(5)). Each other listed m has a nonphysical conjugate at most
cos(5pi/7)<-3/5, contrary to (2). Witness indices are respectively
3,5,5,7,7,9,9,11,13 for m=4,6,7,8,9,10,11,12,15. All are coprime
to 2m. The bound cos(2pi/7)>3/5 follows, for example, by evaluating
s^3+s^2-2s-1 at s=6/5: it is -29/125, and this polynomial increases
above 1 to its root 2cos(2pi/7). Since 1/2+1/n<3/5, the contradiction
is uniform. The rational label3 has x=1/2 and zero determinant.

### E=B or T

Put alpha=2sqrt(2)*x for B and alpha=2x for T. The respective long
cycles make alpha an algebraic integer in K. Their costs satisfy

    B: q_E-1=[n*alpha^2+4alpha-4]/[4(n+1)],
    T: q_E-1=[n*alpha^2+4alpha-4]/[2(n+1)].

The field is exactly Q(alpha). At every other embedding the quadratic
inequality gives

    -7/12 < alpha' < 3/8.                           (3)

Indeed its exact endpoints are (-2 plus or minus 2sqrt(n+1))/n,
and the stated looser interval holds at n=19,20.

If x is finite, B has sqrt(2)<=alpha<2sqrt(2)<20/7. A degree at
least three would have nonzero integer norm of absolute value less
than (20/7)(7/12)^2=35/36. Thus its degree is at most two. In degree
two its integer trace is 1,2 or3, and its nonzero norm is +/-1.
The only real positive roots in this physical interval are the
golden ratio, 1+sqrt(2), and (3+sqrt(5))/2. Their fields are excluded.
The other formal pair (trace3,norm-1) exceeds 2sqrt(2); a repeated
root is not a quadratic case. In degree one alpha=2 and K=Q, also
outside the eight fields.

For T, finite x has 1<=alpha<2. The norm argument is even stronger.
In degree two its trace is 1 or2 and norm +/-1; only the golden
ratio lies in [1,2), and its ground field is excluded. The rational
case alpha=1 gives K=Q. This proves for ALL three far E profiles

    x>1.                                             (4)

The use of cyclic parameters rather than individual unit-normal
entries avoids assuming that sqrt(2) or x separately lies in K in B.

## 5. Far E with a far original-zero D is already excluded

Take E in F,B,T and D from the last three far rows of (1). If their
attachment point is the same, their marked pair E,D has distance
two and retains a baseline unit, or the finite corner is impossible.
This is not a whole A_n selected role: b is an isolated white
spectator at that corner.

If one attaches at tL and the other at t_(L-1), the connected tree
in W+E+D has a long fork with short arm lengths 1 and 2. At these
ranks it is neither A nor D nor a finite E diagram. An added 4-bond
also excludes B, which is a path. Thus it is not a finite corner.
Two simple attachments at t_(L-1) either retain the marked distance-two
unit or make a nonfinite higher-degree fork. Consequently no unpaid
far/far case remains.

## 6. Far E with isolated or near D: a uniform Schur bound

It remains to consider D isolated or attached simply/with a 4-bond
at t1, or simply at t2. Let d=D_bar^2, let -H be the projected
c,D product beyond the ambient -y, and let -v=E_bar.D_bar. Then

|D profile|d|H|v|
|---|---|---|---|
|isolated|1|0|0|
|simple at t1|(N+1)/(2N)|(N-1)/(2N)|h|
|4 at t1|1/N|sqrt(2)(N-1)/(2N)|sqrt(2)h|
|simple at t2|2/N|(N-2)/(2N)|2h|

The a,D product is zero. By (4) the a,E block, with entries
3/4,t,-x, is nonsingular Lorentz. Put Delta=x^2-3t/4>1/2.
After eliminating it, the remaining c,D matrix must be positive
semidefinite of rank one: the full local system has n positive,
one negative and one zero direction, and the negative direction has
already been used. Its entries are

    R_cc=r+[t/16+xh/2+3h^2/4]/Delta,
    R_DD=d+3v^2/(4Delta),
    R_cD=-(y+H)+v[x/4+3h/4]/Delta.

Set Z=v[x/4+3h/4]/Delta. Positive semidefiniteness yields the valid
upper bound, independent of the sign of R_cD,

    y <= sqrt(R_cc R_DD)+Z-H.                        (5)

Because 1/Delta and x/Delta decrease for x>=1, the simple bounds
r<1/3, t<2/3, h<1/16 give

    R_cc < 1/3+1/12+1/16+3/512 <1/2.                (6)

The four D rows now close as follows.

* Isolated: R_DD=1 and Z=H=0, so y<sqrt(1/2). But its opposite
  extension requires y^2>r>1/4; an actual ambient cosine or dotted
  product then has y>=cos(pi/4)=sqrt(1/2), contradiction.
* Simple at t1: R_DD<3/4, Z<1/16 and H>7/16. Equations (5),(6)
  imply y<sqrt(3/8)-3/8<1/2. The value y=0 gives the explicit
  affine E8 with arms 1 (D),2 (c,b),5 (t2,...,t6) at t1. Hence
  compactness requires y>=1/2, contradiction.
* 4 at t1: R_DD<1/8, Z<1/8 and H>1/2. Thus (5) gives
  y<1/4+1/8-1/2<0, impossible even for a right opposite product.
* Simple at t2: R_DD<1/4, Z<1/8 and H>7/16. Thus y<1/16<1/2.
  Its y=0 case has affine E8 at t2, with arms 1 (D),2 (t1,c),
  and5 (t3,...,t7). Again compactness forces y>=1/2.

For the displayed bounds use x/Delta<2, 1/Delta<2. In the second
row v<1/16; in the last two v<1/8. The inequalities d<2/3,1/16,1/8
in their respective rows give exactly the R_DD bounds above. All
affine witnesses are principal: extra actual walls do not undo them.
This excludes every far E with every original-zero D.

## 7. The two nonfar E types: finite arithmetic alternatives

We now have E=I or L. Opposite-side indefiniteness gives respectively
x^2>3/4 or x>1/2. If x is finite, the exact old-A_n costs are

    I: q_E=2n*x^2/(n+1),
    L: q_E=[2n*x^2+2(n-1)x+n-1]/(n+1).

They force K=Q(x^2) in I and K=Q(x) in L. In I every nonphysical
cos(2j*pi/m) is less than 1/n. In L, the factorization

    q_E-1=2(x+1)(nx-1)/(n+1)

gives -1<x'<1/n. Cyclic integrality puts 4x^2 in K in I and 2x
in K in L. The complete finite-label argument is Section 7 of the
same-endpoint source cited in Section 2: degree at most five and
explicit positive conjugates leave

    I: labels12 or18;        L: labels6 or9.          (7)

For clarity this uses the changed A_n conjugate threshold 1/n, not
the stricter zero threshold of the B_n calculation. Every excluded
positive cosine witness is at least sin(pi/30)>1/15>1/n. Rational,
Q(sqrt(2)), Q(sqrt(5)) and K7 alternatives are outside the eight fields.

## 8. Isolated D is excluded by an integral norm

For E=I or L and D isolated, put respectively

    delta=4x^2-3,       or       delta=2x-1,
    w=4y^2.

The corank equation in either case is exactly

    delta[Nw-N-2]=N.                                 (8)

Cyclic integrality gives delta,w in O_K. Reducing (8) modulo N
shows 2delta in N O_K, so delta belongs to d0 O_K with

    d0=N/gcd(2,N)=17 or9.

The physical opposite block gives y^2>r>1/4; actual ambient products
therefore have w>=2. Thus (8) gives 0<delta<=N/(N-2)<2.
In I all nonphysical 4x^2 lie between 0 and 2(n+1)/n<3, giving
|delta'|<3. In L, -1<x'<1/n gives the same |delta'|<3. Every
conjugate of the nonzero algebraic integer delta/d0 has absolute
value less than one. Its nonzero integer norm is impossible. This
settles both isolated/isolated and leaf/isolated cases, without a
field-list norm estimate or a spectator assumption.

## 9. All six endpoint/next-to-end D rows for I or L

Write d=D_bar^2 and -H=c_bar.D_bar+y. For the near/far positions
the exact pairs are

|D contact|d|H|
|---|---|---|
|near simple|1/2+epsilon|1/2-epsilon|
|near 4|2epsilon|sqrt(2)(1/2-epsilon)|
|near next simple|4epsilon|1/2-2epsilon|
|far simple|1/2+epsilon|epsilon|
|far 4|2epsilon|sqrt(2)epsilon|
|far next simple|4epsilon|2epsilon|

Set g=x^2-3/4, kappa=1/16 in I, and g=x-1/2, kappa=1/8 in L.
Eliminating E (or splitting the symmetric a,E pair in L) gives

    (y+H)^2=d(r+R),        R=kappa/g>0.              (9)

From (7), for finite x the least g is (sqrt(3)-1)/4 in I and
(sqrt(3)-1)/2 in L. Consequently in BOTH cases R<3/8. For a
dotted x>1, R<1/4; for finite x<1, R>1/4. The value x=1 is
forbidden. These are all possibilities.

Near simple: y=0 has the affine E8 witness from Section 6. Otherwise
y>=1/2 and y+H>15/16, whereas (9) gives

    (y+H)^2 < (2/3)(1/3+3/8)=17/36<1/2,

a contradiction. Near 4 has H>1/2 and d<1/16, so (9) gives a
square below 17/384<1/4 even before adding y, impossible. Near next
has H>7/16 and d<1/8, so its square exceeds 49/256 but is below
17/192, again impossible.

For the far 4 and far next rows, y=0 fails (9), since H^2<dr.
If y is nonzero, y>=1/2 but the right side is below 17/384 or
17/192 respectively, both below 1/4. Thus both rows are excluded.

For far simple put t=1/2+epsilon; here H=epsilon and r+1/4=t.
If x is dotted, R<1/4 gives y+epsilon<t, hence y<1/2. But y=0
fails because epsilon^2<dr, and no positive ambient value is below
1/2. If x is finite, R>1/4 gives y>1/2, whereas (9) and the same
bound 17/36<1/2 give y<sqrt(1/2). This is precisely the impossible
gap between cos(pi/3) and cos(pi/4). All six D rows are excluded.

## 10. Completion, composition and exact limitation

Section 2 exhausted E and D. Its bridge cases J_near,J_far are
excluded by the explicitly linked same-endpoint proof. Sections 4-6
exclude every far E with every original-zero D. Sections 7-9 exclude
the remaining I,L cases with every original-zero D. Hence no newly
zero quad exists, including all simultaneous selected corners.

Triangles retain their unchanged high-label I2 payment of two.
All weights remain nonnegative, so larger polygons retain their
required nonpositive payment. This proves the proposed consumer
theorem, subject to independent written review of this full proof
and its same-endpoint dependency.

Its white graph A1+A_(n-3) differs from the connected A_(n-2)
white graph of the full-D_n whisker deletion and from the 4-bond
white graphs of both B_n reduced roles. A face cannot carry two
different fixed white graphs. Thus no quad mixes these negative
changes. Their nonnegative recipient coupons only increase payments.
This gives the stated composition within the eight fields.

At full A20 the two distinct deleted baseline units reduce 37 to35.
No deletion at A19+A1, D19+A1 or other outside-factor contexts is
being inferred. The reviewed B rules still have costs 141/4 at
B20 and B19+A1; the D20 and other forest obligations remain. There
is no global allowance35, arithmetic-dimension19 theorem, new example,
formalization or computation in this report.
