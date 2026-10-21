# Same-endpoint simultaneous A19/A20 end-mark deletions cannot leave an unpaid quad

19 October 2026. **New scoped written reduction, awaiting independent
review.** No scientific job, proof script, optimizer, Lean process or build
was run. Stable B-rule and survivor sources are unchanged.

This note does not prove an A_n refund. It removes one entire branch
of its simultaneous-consumer problem, retaining every permitted contact
of the new walls with the long white path. The unresolved one-selected
corner and third-root-shared branches remain explicit below.

## 1. Proposed rule and the precise reduction

Let P be a compact arithmetic simplest-type Coxeter polytope of dimension
n=19 or20, with ground field one of the eight high-label fields

    Q(sqrt(3)), Q(sqrt(6)), Q(sqrt(2),sqrt(3)),
    K9,K11,K15,K16,K20,       Km=Q(cos(2pi/m)).

Use the distance-two baseline plus the entire-I2(m>=6) addition from
`results/fresh_arithmetic_high_label_20261016.md`. Consider the PROPOSED
deletion of the two simple-end distance-two units at every entire A_n
vertex: mark its first and third roots at each end of the path.

At one such marked corner label the simple end a-b-c and the rest
of the path t1-...-tL, L=n-3. Thus a is the endpoint and c the third
root. In an actual quadrilateral put the four black roots in cyclic
order a,c,E,D and fix the common white roots

    W={b}+{t1,...,tL}=A1+A_(n-3).

**Same-endpoint reduction.** If the adjacent corner (a,D) is ALSO a
selected A_n end-distance-two role, this quad cannot have total zero
after all proposed A_n deletions. This includes either orientation
of the long white path, every finite continuation for E, and all
simultaneous proposed deletions at the other corners.

The hypotheses do not assert that two such selected corners can never
lie on an actual quad. The conclusion excludes their being an UNPAID
consumer after the proposed deletion.

## 2. Complete finite continuation inventory under zero charge

Suppose for contradiction that this quad is unpaid. Every adjacent
black product other than the original a,c product is zero. A nonright
adjacency would retain a baseline unit, whereas the selected A_n
distance-two marks are orthogonal pairs. Thus c.E=E.D=D.a=0.

Because a is still a leaf at the corner W+a+D, a selected end mark
(a,D) must put D at the third position. Therefore D contacts b by a
simple bond and contacts exactly ONE endpoint of the long white path
by a simple bond. There are precisely two types:

    J_near: D.b=D.t1=-1/2;
    J_far:  D.b=D.tL=-1/2.

No other D-white contact is possible in the entire finite A_n path.

The corner W+c+E contains the A_(n-1) path b-c-t1-...-tL.
A finite one-root extension at these ranks is isolated, a simple
or 4 extension at an endpoint, or a simple extension at a next-to-end
root, giving respectively A,B,D. Exceptional finite types are too
small. The zero-modified-charge condition leaves exactly these E types:

| E type | Contacts with W+c |
|---|---|
|I|none|
|L|simple contact at b|
|F|simple contact at tL|
|B|4-bond contact at tL|
|T|simple contact at t_(L-1)|

The L case is the simultaneous selected corner (c,E), with endpoint E.
A 4-bond at b gives an undeleted distance-two unit in a whole B_n,
and a simple contact at c gives an undeleted adjacent unit at a whole
D_n fork; neither can occur in an unpaid quad. This derives the full
five-row inventory without assuming that the far contacts vanish.

The earlier B_n reductions cannot create additional zero roles here:
their white graphs contain a 4-bond, unlike W. A whole-D_n whisker
deletion has connected white graph A_(n-2), also different from W.
Existing nonnegative coupons may be ignored when deriving this necessary
inventory; retaining them could only rule out more unpaid cases.

## 3. Projected parameters and compactness facts

Put N=n-2 in {17,18}, epsilon=1/(2N), and

    r=1/4+epsilon,       q=1/2-1/N.

After eliminating W, the old norms and product are

    a_bar^2=3/4,       c_bar^2=r,       a_bar.c_bar=-1/4.

Write the opposite ambient products as a.E=-x and c.D=-y, x,y>=0.
Opposite-side indefiniteness gives

    E=I: x^2>3/4;         E=L: x>1/2.

For the far E rows it gives x^2>(3/4)E_bar^2. These follow from actual
opposite sides: a positive-definite opposite extension would produce
an actual intersection, and a singular positive-semidefinite one is
forbidden in a compact chamber. In particular every x used below is
strictly positive.

An ambient actual nonzero product is a finite cosine at least 1/2,
or a dotted product at least one. The value x=1 is impossible in a
compact acute chamber, since the positive null sum of its two actual
wall normals would be a feasible ideal point.

In J_near, y=0 is also impossible. The four actual roots b,c,t1,D
would form the induced all-simple cycle b-c-t1-D-b, with positive
null kernel (1,1,1,1). This is an explicit affine A3 subdiagram.
Thus y>=1/2 in every J_near case.

## 4. J_near with all three far E continuations

The roots c,D have identical white contacts. Splitting into their
symmetric and antisymmetric combinations reduces the corank equation
to

    (y+q)(x^2-3t/4)=3h^2/2+t/8+xh,                (1)

where t=E_bar^2 and -h is E_bar.c_bar=E_bar.D_bar. Direct inversion
of the white A_(n-3) path gives

|E type|t|h|
|---|---|---|
|F|(N+1)/(2N)|1/(2N)|
|B|1/N|sqrt(2)/(2N)|
|T|2/N|1/N|

All shared-white cross terms are included. In particular the 4-bond
and next-to-end rows have different norms and cannot be treated as
white-isolated roots.

For x>sqrt(3t/4), equation (1) makes

    y(x)=(3h^2/2+t/8+xh)/(x^2-3t/4)-q

strictly decreasing: the numerator of its derivative is
-h*x^2-2(3h^2/2+t/8)x-3ht/4<0.

For F, the least permitted ambient x is cos(pi/4), because its
opposite threshold exceeds (1/2)^2. For B and T it is 1/2. At
these least permitted values the exact upper bounds are

|E type|minimal ambient x|upper bound for y|
|---|---|---|
|F|sqrt(2)/2|(3+2sqrt(2))/(N-3)|
|B|1/2|(-N+7+2sqrt(2))/(2(N-3))|
|T|1/2|(14-N)/(2(N-6))|

For N=17,18 the first is strictly below 1/2 and the other two are
negative. Since compactness gave y>=1/2, all three cases are excluded.

## 5. The two remaining J_near rows

For E=I the corank equation is

    (y+q)(x^2-3/4)=1/8.

Since y>=1/2, it forces

    3/4<x^2<=3/4+N/(8(N-1))<=113/128<1.

Thus x is a finite cosine. The value cos^2(pi/9) is strictly larger
than 113/128. To check the close comparison exactly, let
s=2cos(2pi/9)>1. It satisfies s^3-3s+1=0, whose value at 49/32 is
-111/32768. The polynomial is increasing above 1, so s>49/32 and
cos^2(pi/9)=(2+s)/4>113/128. Consequently x has label 7 or 8.

The independent old-A_n-plus-E determinant ratio has projection
cost 2n*x^2/(n+1)>1 physically and less than one at every other
embedding. It therefore forces K=Q(x^2): an embedding fixing x^2
would preserve the negative physical determinant. Labels 7 and 8
give K7 and Q(sqrt(2)), neither one of the eight fields. This excludes I.

For E=L, the corresponding equation is

    (x-1/2)(y+q)=1/4.

It gives 1/2<x<=49/64, so ambient cosine discreteness leaves only
x=cos(pi/4). Substitution gives y=sqrt(2)/2+1/N, strictly between
cos(pi/4) and cos(pi/5). For example sqrt(2)/2<5/7 and 1/N<1/16
give y<87/112<4/5<cos(pi/5). This is not an ambient Coxeter product.
Thus J_near is completely excluded, including four simultaneous marks.

## 6. J_far has only two finite E continuations

The three far rows F,B,T cannot give a finite opposite corner W+E+D.
For F its large connected tree has a fork at tL with arms
1 (E), 2 (D,b), and L-1. It is neither an A path nor D with two
unit arms; its rank is too large for E6,E7,E8. B has the same
branch with a 4-bond and cannot be a finite B path or exceptional
diagram. For T the fork is at t_(L-1), with arms 1 (E), 3 (tL,D,b)
and L-2, again nonfinite. These are actual finite-corner tests,
not conjectural contacts outside the face.

Only E=I or L remains. Their corank equations are respectively

    y(x^2-3/4)=1/8,        y(x-1/2)=1/4.             (2)

They force y>0, hence y>=1/2, and x<=1. The parallel value x=1 is
forbidden, so x is finite. The field alternatives needed for (2)
are recorded next.

## 7. Finite-angle alternatives with the full A_n projection costs

For I, the independent old-A_n-plus-E projection cost is

    q_E=2n*x^2/(n+1).

It forces K=Q(x^2), and every nonphysical squared cosine satisfies
0<(x^2)'<(n+1)/(2n). Hence its conjugate cos(2j*pi/m) is less than
1/n. Since each of the eight fields has degree at most five, the
finite labels m>=7 are drawn from

    7,8,9,10,11,12,14,15,16,18,20,22,24,30.

Labels 7,14,8,10 give the excluded ground fields K7,Q(sqrt(2)),
Q(sqrt(5)). For each of 9,11,15,16,20,22,24,30 there is a positive
nonphysical conjugate cosine at least sin(pi/30)>1/15>1/n. Explicit
conjugation indices j are respectively 2,2,2,3,3,3,5,7. This leaves
only m=12 or18 for I.

For L the exact whole-A_n projection cost is

    q_E=[2n*x^2+2(n-1)x+n-1]/(n+1),
    q_E-1=2(x+1)(nx-1)/(n+1).

The cycle a-b-E puts x in K. The independent determinant ratio
forces K=Q(x) and -1<x'<1/n at every other embedding. Degree at
most five lists m=4,5,6,7,8,9,10,11,12,15. Labels 4,5 give the
excluded quadratic fields. Each of 7,8,10,11,12,15 has a positive
nonphysical cosine at least sin(pi/30)>1/n (indices 3,3,3,3,5,7).
This leaves only m=6 or9 for L. These are the same finite labels
as in the B_n analysis, but the changed A_n projection costs and
the positive allowance 1/n are explicitly included here.

## 8. The final J_far angle gaps

In I, equation (2) and labels 12,18 give

    1/2<y<=1/[2(sqrt(3)-1)]=(sqrt(3)+1)/4<sqrt(2)/2.

The upper value occurs at label 12; label 18 has larger x and a
smaller y, still above 1/2 because x<1. This is the impossible
ambient interval between labels 3 and 4.

In L, the labels 6,9 give exactly the same inequality from the
second equation (2). Label 6 gives y=(sqrt(3)+1)/4; label 9 has
larger x and smaller y, again above 1/2. Thus J_far is also excluded.
This completes the same-endpoint reduction.

## Scope and remaining consumers

No unpaid quad under the proposed full deletion can have an adjacent
selected corner sharing the endpoint a. The proof covers both bridge
orientations and every finite E continuation, including all far
attachments and all resulting simultaneous marks. In particular it
excludes the four-selected-corner configurations in this branch.

A quad with only the old selected corner, or with the adjacent
selected corner (c,E) sharing the THIRD root c while (a,D) is
originally unpaid, is not covered. The original-zero D possibilities
include isolated, endpoint 3/4 and next-to-end 3 continuations of
the long white A path; they have not been set to zero or discarded.
These remaining branches are essential before even a fractional
A20 refund can be asserted.

The two A20 end-distance-two baseline units would lower 37 to35 if
a full consumer theorem existed. This report proves only the stated
branch reduction, and does not assert that theorem, allowance35 or
any dimension exclusion. No charge implementation or formal source
was changed.
