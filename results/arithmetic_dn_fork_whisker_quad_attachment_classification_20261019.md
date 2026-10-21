# Completing the isolated-wall fork/whisker quad classification in D19/D20

19 October 2026. **New written research proof, awaiting independent
review.** No enumeration, numerical solve, proof script, Lean process
or TeX build was run. Earlier reviewed reports are unchanged.

The finite-continuation analysis below closes the remaining finite
ambient-product cases in the isolated-wall branch. Together with the
two reviewed local exclusions, it proves the following restricted
statement. It does not supply a global refund or a new dimension bound.

**Claim.** Let n be 19 or 20 and let an actual finite vertex of a compact
arithmetic simplest-type Coxeter polytope have full diagram D_n. Let
a be its fork and b,e its short whiskers. Suppose an actual quadrilateral
at the angle (a,b) has opposite side pairs (a,c),(b,d), with

    d perpendicular to every old D_n root except b,
    c perpendicular to b and d.

Then this quadrilateral cannot occur. No products of c with a,e or
the long arm are initially set to zero. Their complete finite
continuations are classified below.

The already reviewed dotted/parallel a.c branch uses the published
high-dimensional field restriction through
[the isolated-quad obstruction](arithmetic_dn_isolated_quad_obstruction_20261018.md).
The new finite-product arguments below do not use a published field
list. The zero a.c branch uses, among its cases,
[the reviewed double-whisker lemma](arithmetic_dn_double_whisker_obstruction_20261018.md).
These inherited inputs are distinguished from the new classification.

## 1. Coordinates, arithmetic inputs, and the actual opposite pair

Put N=n−2 and L=n−3=N−1. Thus N is 17 or 18. Write the long arm as

    b
    |
    a--l1--l2-- ... --lL,
    |
    e

All old bonds have label 3. Normalize roots to square 1. The complete
white support W is {e,l1,...,lL}, of type A1+A_L. It has rank n−2.
Every corner of the actual quadrilateral is finite. In particular W+c
is positive definite, as it is a principal submatrix of W+b+c.

We use the ordinary acute-wall Gram property and finite Coxeter
classification. An ambient product between distinct facet normals is
zero, a negative finite Coxeter cosine, or at most −1. We do not infer
the last alternative merely from opposition in the two-face.

Here b and d are perpendicular to W, so their projections to its
Lorentz three-dimensional perpendicular space are unchanged. Those
two projected sides are opposite sides of an acute compact polygon.
If their two-root Gram were positive definite, projection of an
interior polygon point perpendicular to both would remain timelike
and feasible against all other sides: the two-root inverse is
entrywise nonnegative, and the correction is a nonnegative combination
of the two outward normals. This would give a common polygon corner,
contrary to opposition. A positive null sum excludes the parallel
case by compactness. Consequently

    b.d=−y,  y>1.

This argument concerns the unchanged b,d projections. It is not an
assertion that the ambient a,c walls are disjoint.

Arithmetic inputs are cyclic integrality in the totally real ground
field K and positive definiteness of independent root blocks at every
nonphysical embedding. Unit-normal entries need not belong to K.
Squared contacts and closed cycle products do. Equations involving
embeddings below are determinant-ratio identities, with positive
diagonal normalization of arithmetic roots at each embedding.

Whenever we exclude a positive-null induced subdiagram, the reason is
the following compactness argument. If positive coefficients v satisfy
Gv=0 for actual outward normals, their sum z is nonzero because it
pairs negatively with an interior timelike point. It is null and
orthogonal to those normals. Acuteness makes every other wall
nonpositive against z, giving a feasible ideal ray, impossible in a
compact chamber. For a square-two simply laced tree, the displayed
positive-null certificates below satisfy twice each coefficient equals
the sum of its neighbors.

## 2. Complete finite continuations of W+c

Because A_L has at least sixteen vertices, a finite irreducible
component containing it can only be classical. A contact of c with
two long-arm vertices would make a cycle and is impossible in a
finite diagram. The complete possibilities are:

|Long-arm contact|Contact to e|Permitted bond(s)|Reason|
|---|---|---|---|
|None|None|None|c is isolated in W+c|
|None|e|Any finite label ell≥3|One finite dihedral factor|
|l1 or lL|None|3 or 4|An A- or B-series path|
|l2 or l_(L−1)|None|3|A D-series fork|
|l1 or lL|e|Long bond 3, e bond 3 or 4|An A- or B-series path|

There are no other cases. For an internal long-arm attachment, the
three arm lengths at that attachment must have two length-one arms
in the D series; hence the attachment is next to an endpoint. If c
also contacts e, its new arm has length two, which prevents a large
D-series fork; the exceptional finite E diagrams are too small.
For a path, a 4-bond must be at an endpoint in these large ranks.
Thus a bond 4 on the long contact cannot coexist with the further
e neighbor. Larger finite labels cannot occur in a connected finite
component of these ranks. The dihedral-only e contact is deliberately
left unrestricted in the table.

This is a structural classification of all native finite contacts,
not a list obtained by a label cutoff. Its use of the full white
support is essential; outside old components are not being discarded.

## 3. The zero a.c branch is completely excluded

Assume a.c=0. If c has no W contact either, it is perpendicular to
the full positive n-dimensional D_n span. Its perpendicular is a
negative line, which cannot contain a unit spacelike root.

If c contacts only e, the reviewed double-whisker theorem applies
directly. Its finite e contact has any native label ell≥3, and the
other isolated contact is the actual dotted b/d pair proved above.

Suppose c contacts a long vertex l_k and is perpendicular to e.
For a simple bond and k<L, take the induced subtree with joining
path a-l1-...-l_k and leaves b,e at a, c,l_(k+1) at l_k. Assign
coefficient 2 along the joining path and 1 to the four leaves. This
is a positive-null affine D diagram. At k=L, a simple bond instead
makes the complete D_(n+1) Gram positive definite, impossible in
positive index n.

For a 4-bond, the table permits only k=1 or L. Take a,l1,...,l_k,c
and b,e. Give a and all l_i coefficient 2, b,e coefficient 1, and
c coefficient sqrt2 in the square-two gauge. The last bond has
Gram entry −sqrt2. These coefficients give a positive null vector.
This is an induced affine B diagram; unselected later long vertices
do not create extra edges in the selected block.

Now let c contact both e and a long endpoint by a simple long bond.
If c.e is also simple, the induced cycle

    a-e-c-l_k-l_(k−1)-...-l1-a

has all bonds 3 and positive-null vector with all coefficients 1.
Here a.c=0 is needed so that there is no chord.

If c.e has label 4 and k=L, use the induced five-node diagram on
{a,b,l1,e,c}. Its fork a has leaves b,l1 and tail a-e-c ending in
a 4-bond. Coefficients 2 at a,e, 1 at b,l1, and sqrt2 at c give
a positive null vector. The c-lL edge is absent from this selected
block because L≥16.

The only remaining zero-product case is c.e of label 4 and c.l1 of
label 3. Its old-path cycle puts sqrt2 in K. Its D_n projection cost
at the embedding which negates sqrt2 is

    q'=n/4+(n−3)/2−(n−3)/sqrt2.

It is greater than 1 for n≥9: equivalently
(3n−10)^2>8(n−3)^2, whose difference is n²−12n+28>0 in that range.
The relevant quantities on both sides before squaring are positive.
This contradicts arithmetic conjugate positivity. The D_n+c block
is physically independent because its physical projection cost is
greater than 1. Thus the entire a.c=0 branch is excluded.

## 4. Finite nonzero a.c: common formulas

The reviewed isolated-quad lemma excludes a.c≤−1 for n=19,20. It
remains to consider

    a.c=−x,  x=cos(pi/m),  m≥3 finite.

Let alpha be the magnitude of the long contact, if present, and beta
the magnitude of c.e, if present. Missing contacts have magnitude zero.
Schur-eliminating W gives

    s=n/(4N),
    t=1−beta²−alpha²*2k(N−k)/N,
    kappa=beta/2+alpha*(N−k)/N,
    X=x+kappa.

When there is no long contact, omit the terms containing alpha,k.
The finite continuation gives t>0. Eliminating the complete old D_n
from c gives projection cost

    q_c=1−t+2N*(x+kappa)².                         (7)

The complete n+2-root corank equation is

    (X²/t−s)(y²−1)=1/4,
    z:=4y²=4+t/(X²−s*t).                          (8)

The denominator is nonzero because the first equation has right
side 1/4. The two-cycle z is an algebraic integer in K.

Whenever t and kappa are rational and a closed simple-bond path
puts x in K, admissibility of the independent D_n+c block implies

    |sigma(x)+kappa|<r,  r²=t/(2N)                 (9)

at every nonphysical embedding. Independence also follows directly
in all the cases used below from the physical value q_c>1.

## 5. A small cyclotomic moment lemma

Let x=cos(pi/m), m≥3, and d=[Q(x):Q]=phi(2m)/2. Exact root-of-unity
traces give

    Tr(x)=mu(2m)/2,
    Tr(x²)=d/2+mu(m)/4  if m is odd,
    Tr(x²)=d/2+mu(m)/2  if m is even.

For odd m, mu(2m)=−mu(m); for even m, mu(2m)=0. Hence for
0≤kappa≤3/4,

    Tr((x+kappa)²)≥d*(1/2+kappa²)−1/2.

If every nonphysical conjugate satisfies (9), its upper bound is

    Tr((x+kappa)²)<(1+kappa)²+(d−1)*r².

Consequently d≤3 whenever

    B:=1/2+3*kappa²−2*kappa−3*r²>0.               (10)

Indeed at d=4 the lower trace bound exceeds the upper by B; increasing
d increases that difference since 1/2+kappa²−r²>0 in our cases.
The complete list phi(2m)≤6, m≥3, is

    m=3,4,5,6,7,9.

All nonidentity embeddings of Q(x) lift to nonphysical K embeddings
in applications where the actual cycle puts x in K. This is why the
cycle hypothesis in this lemma is explicitly retained.

## 6. Single simple long contact, with or without a simple e contact

First suppose c.e=0 and the long bond is 3. The table permits the
four positions below. The cycle consisting of a,c and the old long
path puts 2x, and hence x, in K. Here are all centering constants:

|k|kappa|t|B from (10)|
|---|---|---|---|
|1|(N−1)/(2N)|(N+1)/(2N)|(N−5)/(4N)|
|2|(N−2)/(2N)|2/N|(N−4)/(4N)|
|L−1=N−2|1/N|2/N|(N−4)/(2N)|
|L=N−1|1/(2N)|(N+1)/(2N)|(2N−7)/(4N)|

Every B is positive. Also r<1/8 for N=17,18; the near kappa values
lie between 7/16 and 1/2, and the far ones lie between 0 and 1/16.
Thus the moment lemma leaves only the six labels above.

For m=3 the induced cycle through a,c and the old long path has all
bonds 3 and is parabolic. For m=4,6 one conjugate of x is at most
−sqrt2/2; the same is true of cos(7pi/9) when m=9. Adding any
kappa<1/2 leaves absolute value greater than 1/8. For m=7 the
conjugate cos(3pi/7)=sin(pi/14)>1/7 is already greater than r;
the sine bound follows from its strict concavity on [0,pi/2].
For m=5 its other conjugate is (1−sqrt5)/4. At a far kappa it is
less than −1/4+1/16=−3/16. At a near kappa it is greater than
−5/16+7/16=1/8, using sqrt5<9/4. All violate (9).
Thus every simple-long-only finite-product case is excluded.

Now let the long bond and c.e both have label 3. There are just the
two endpoint positions. They have

|Position|kappa|t|B|
|---|---|---|---|
|Near k=1|(3N−2)/(4N)|(N+2)/(4N)|(11N−26)/(16N)|
|Far k=L|(N+2)/(4N)|(N+2)/(4N)|(3N−10)/(16N)|

Again d≤3 and r<1/10. Near kappa lies between 7/10 and 3/4; far
kappa lies between 1/4 and 3/10. The same six-label check leaves only
m=4 in the near case and m=5 in the far case. More explicitly m=3
gives the parabolic triangle a,e,c; m=7 has a positive conjugate;
m=6 has a conjugate below −6/7; m=9 has a conjugate cos(5pi/9)>−1/2
which excludes the near case and cos(7pi/9)<−sqrt2/2 which excludes
the far case. The m=5 other conjugate excludes the near case, and
the m=4 negative conjugate excludes the far case.

The last two possibilities fail integrality in (8). Directly using
the two conjugates gives

    Near, m=4: Tr(z)=4+4(3n−8)/(n²−8),
    Far,  m=5: Tr(z)=4+2(7n−8)/(n²+2n−4).

Both traces are strictly between 4 and 5 at n=19,20. Their denominators
exceed the displayed positive numerators by n²−12n+24 and n²−12n+12,
respectively. A quadratic-field trace of an algebraic integer is an
integer, so both cases are impossible.

## 7. A single long 4-bond and no e contact

Put h=sqrt2*x. The actual old-path cycle makes w=2h an algebraic
integer in K. Formula (7), with t=1/N, becomes

    q_c=1−1/N+N*(h+kappa0)²,
    kappa0=(N−k)/N.

Thus every nonphysical conjugate satisfies

    |sigma(h)+kappa0|<1/N.                         (11)

At the far endpoint kappa0=1/N, so all nonphysical conjugates of w
have absolute value less than 4/N. Its physical value is in
[sqrt2,2sqrt2). If its degree were at least 2, its nonzero integer
norm would have absolute value less than 3*(4/17)=12/17<1. Hence
w is a rational integer, necessarily 2. This gives m=4, and (8) gives

    z=4+4/(2n−1),

which is rational and nonintegral.

At the near endpoint kappa0=(N−1)/N, set A=w+2. Its physical value
lies in [2+sqrt2,2+2sqrt2), below 5. Every other conjugate lies in
(0,4/N). Degree at least 3 would give norm less than
5*(4/17)^2<1. In degree 2 its positive integer norm is 1 and its
integer trace is 4 or 5. The two possibilities imply, respectively,

    x²=3/8,  or  x²=(11+sqrt21)/16.

Neither is a native finite cosine-square. In the first case its real
cyclotomic field is rational, whose m≥3 possibilities are 3,4,6 and
give 1/4,1/2,3/4. In the second case it is quadratic with field
Q(sqrt21), whereas the complete quadratic cosine-square labels
5,8,10,12 have fields Q(sqrt5),Q(sqrt2),Q(sqrt5),Q(sqrt3).
This uses only phi(m)≤4, not a classification of arithmetic fields.

In degree 1, A=4, giving m=4 and

    z=4+4/(8n−25),

again rational and nonintegral. Both long-4 cases are therefore excluded.

## 8. A simple long bond and an e 4-bond

Here two actual cycles put sqrt2 and x in K. Formula (7) has

    t=1/(2N),  kappa=kappa0+sqrt2/4,
    kappa0=(N−k)/(2N),

and conjugate positivity requires

    |sigma(x)+kappa0+sigma(sqrt2)/4|<1/(2N).       (12)

If 4 does not divide m, a cyclotomic automorphism can fix x and
negate sqrt2. To see this directly, take the common cyclotomic field
of conductors 2m and 8 and use the Chinese remainder theorem to choose
an exponent congruent to 1 modulo 2m and 5 modulo 8. The congruences
are compatible because the 2-part of 2m is at most 4. Its restriction
to Q(x,sqrt2), and then an extension to K, supplies the required
nonphysical embedding. The expression in (12) is then at least

    1/2−sqrt2/4 > 1/8 > 1/(2N).

This is impossible.

If 4 divides m, the exponent m−1 is coprime to 2m and sends x to −x.
For m congruent to 4 modulo 8 it also negates sqrt2. Since kappa0<1/2
and x≥sqrt2/2, this violates (12) at once. If 8 divides m it fixes
sqrt2, but x≥cos(pi/8)>9/10, while

    kappa0+sqrt2/4 < 1/2+5/14=6/7.

The absolute difference is greater than 3/70>1/34≥1/(2N), again
impossible. Thus no label m remains in this whole case.

## 9. No long-arm contact: the sole-fork case

If c.e=0 as well, its only old contact is −x at a. The two-cycle
v=4x² is an algebraic integer in K. Admissibility gives all its
nonphysical conjugates in (0,2/N). The physical v lies in [1,4).
Degree at least 2 would give a nonzero integer norm below
4*(2/17)<1. Therefore v is the rational integer 1,2 or 3, giving
m=3,4,6. Substitution in (8) gives respectively

    z=8−2n,
    z=8+8/(n−4),
    z=6+2/(n−3).

The first is negative; the others are rational nonintegers at n=19,20.
This excludes the case without assuming x itself belongs to K.

## 10. No long-arm contact: c supported on a and e

Write beta=cos(pi/ell), ell≥3. Define

    u=x², v=beta², p=x*beta.

The two-cycles and triangle a-e-c give u,v,p in K and p²=u*v.
The projection cost is

    q_c=2N*u+2N*p+(n/2)*v
       =2N*(x+beta/2)²+beta².                    (13)

At any nonphysical embedding, completing the beta-square gives

    q_c >= (4N/n)*sigma(u).

This is a statement about squared contacts and their product; it
does not require separate embedded values of x,beta in K. One may
choose real square roots with the sign of sigma(p), since
sigma(p)²=sigma(u)*sigma(v). Hence

    0<sigma(u)<n/(4N)≤2/7.                        (14)

For m≥4, physical u≥1/2. Any nonidentity K embedding fixing u would
make q_c≥2N/n>1. Total reality therefore forces K=Q(u).

The exact cosine-square trace is d/2+mu(m)/4, where d=phi(m)/2.
If d≥5, its lower bound d/2−1/4 exceeds
1+(d−1)*2/7, contrary to (14) and the physical u<1. Thus d≤4.
The complete labels phi(m)≤8, m≥3, are

    3,4,5,6,7,8,9,10,12,14,15,16,18,20,24,30.

The following labels all have a nonphysical conjugate
cos²(a*pi/m) at least cos²(3pi/10)=(5−sqrt5)/8>2/7:

|m|7|9|10|14|15|16|18|20|24|30|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|a|2|2|3|3|2|3|5|3|5|7|

Each a is coprime to m and a/m≤3/10. Thus, apart from the rational
m=4,6, only m=5,8,12 remain.

For m=4 or6, K=Q. Since v is rational, its finite labels are ell=3,4,6.
The additional condition p in Q permits only ell=m. Formula (8) gives

    m=ell=4: z=4+2N/(4n−9),
    m=ell=6: z=4+2N/(13n−27),

both rational nonintegers.

For m=5, K=Q(sqrt5), and x belongs to K. Hence beta belongs to K.
Its only possible labels are ell=3 or5. Indeed v in K permits
3,4,6,5,10; labels 4,6 have the wrong quadratic root field, and
for ell=10 its v has norm 5/16 and cannot be a square in K.

For m=8, K=Q(sqrt2), and v in K permits ell=3,4,6,8. The norms
of u*v for the first three are 1/128,1/32,9/128, respectively.
None is a rational square, so u*v cannot be a square in K, whereas
p²=u*v. Thus only ell=8 remains.

For m=12, K=Q(sqrt3), and v in K permits ell=3,4,6,12. The identity

    x=sqrt2*(sqrt3+1)/4

puts x in the nontrivial square class sqrt2*K. The condition x*beta
in K therefore permits only ell=4 or12.

Equal labels m=ell=5,8,12 are impossible by conjugate positivity.
At their other embedding, u'>1/16 and (13) becomes

    q_c'=(9n−16)*u'/2 > 1.

The inequality u'>1/16 follows directly from the explicit quadratic
cosine-square values for all three labels.

Only (m,ell)=(5,3) and (12,4) remain. In the first case (8) gives

    Tr(z)=8−36(n−2)(n−3)/(11n²−26n−1),

equal to 4504/869 at n=19 and 2224/431 at n=20. Both are strictly
between 5 and 6. In the second case it gives

    Tr(z)=4+4(7n−15)/(3n²−6n−1),

strictly between 4 and 5 for n=19,20. These traces cannot belong to
an algebraic integer. This completes the m≥4 branch.

Finally take m=3, so x=1/2. The actual triangle cycle puts beta in K.
For every nonphysical conjugate of beta, (13) gives

    q_c−1=(beta+1)*((N+2)*beta+N−2)/2<0.

Every finite cosine conjugate lies in (−1,1), so this forces
beta<−(n−4)/n. The algebraic integer A=2(beta+1) has physical value
in [3,4) and all other conjugates in (0,8/n). Degree at least 3
would give norm less than 4*(8/19)^2<1. In degree 2 its norm is 1
and its trace is 4; hence beta=sqrt3/2 and ell=6. Degree 1 gives
beta=1/2 and ell=3. The latter is the parabolic triangle a,e,c.
For ell=6 the trace of z is the same nonintegral expression

    4+4(7n−15)/(3n²−6n−1).

Thus this last finite-label case is excluded too.

## 11. What has and has not been proved

The table in Section 2 covers every c attachment permitted by the
actual finite W+c corner. Sections 3–10 exclude its zero and finite
nonzero a.c cases without a finite-label cutoff. The previously
reviewed isolated-quad theorem handles a.c≤−1 in dimensions 19 and20,
using its stated published field input. Therefore the claim at the
start follows for the precise isolated-d/orthogonal-c contact pattern.

This does not prove that every unpaid fork/whisker quadrilateral has
an isolated d, or that c.b=c.d=0 under every current payment priority.
Those contact hypotheses remain explicit. A universal refund still
needs all remaining d attachments and all simultaneous triangle and
quadrilateral consumers. The current compact arithmetic ceiling20
and the unconditional compact ceiling29 are unchanged by this report.

No executable certificate or formal verification is claimed. All new
calculations above are written algebra and finite-diagram reasoning;
the report is staged for independent adversarial mathematical review.
