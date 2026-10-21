# A finite upper interval at the opposite A7 edge

20 October 2026. Paper algebra; no scientific execution. This is a successor
to the scalar-facet note, not an amendment to its frozen bytes.

Let s=sqrt3, epsilon=2+s, m=4+s, and use the integral lattice

    M={ (t,z) in O^10 : sum z_i=m*t },
    Q=sum z_i^2-epsilon*t^2.

The scalar-facet note proves that its complete primitive integral-reflector
norm representatives are2 and2epsilon. Its initial chamber corner has A8
walls n_i=-e_i+e_(i+1), i=1,...,8. We now release n8, keeping A7 on the
FIRST eight coordinates. This is the opposite end to the n1 edge initially
described in that note; the orientation matters.

## Root and finite upper endpoint

The root proposed by the coordinating agent is

    T0=5+3s,
    d_z=((4+2s)^3,(3+2s)^5,2+s).

Its sum is29+17s=m*T0. Its squared spatial norm is196+112s, whereas
epsilon*T0^2=194+112s, so Q(d)=2. Thus its reflection is integral on M and
on the original L. Its products with all eight old A8 simple roots are
zero except

    <d,n3>=-1,       <d,n8>=-(1+s).

It is acute to every initial wall. The retained A7 plus d is an E8 spherical
diagram: d attaches at position3 with a3-bond. This proves a finite upper
collision, not that d is the first wall.

Write the outgoing edge in the t=1 slice as

    z1=...=z8=(m+delta)/9,    z9=(m-8delta)/9,    delta>0.

Its square is(1-s)/9+8delta^2/9. For any root of time T and first-eight sum
S8, define Z=9S8-8m*T. Its product with the edge point is

    [T(1-s)+delta*Z]/9.

For d, Z0=11+8s, and

    Z0^2-8(s-1)T0^2=9,
    delta0=(s-1)T0/Z0=(4+10s)/71.

The point is strictly before the light cone because
`delta0^2<(s-1)/8`. Therefore the first mirror collision on this actual
edge exists at a finite point at or before delta0.

## Complete finite profile inventory and direct bounds

The retained A7 permits exactly the following necessary simple extensions:

| Type | q | gap g | high-coordinate count k in the first eight |
|---|---|---|---|
| A7+A1 |2 or2epsilon|0|0|
| A8 |2|1|1 or7|
| D8 |2|1|2 or6|
| E8 |2|1|3 or5|
| B8 |2epsilon|1+s|1 or7|

Here first-eight coordinates are k copies of p+g and8-k copies of p.
Put X=2T. The root equation is

    Z^2-D*X^2=R,       D=2(s-1),
    R=72q-9k(8-k)g^2.

The six profile choices give R144,144epsilon,81,36,9,18epsilon. Required
coordinate congruences are

    T=X/2 in O,
    S8=(Z+8m*T)/9 in O,
    p=(S8-k*g)/8 in O,
    z9=m*T-S8.

The old released wall requires z9<=p. There is NO nonnegativity condition on
p or z9. All eight old wall products must be tested, including n8.

For any positive T,Z, collision position is delta=(s-1)T/Z. The finite upper
interval therefore implies X/Z<=X0/Z0, where X0=2T0=10+6s. Using the norm
equation and Z0^2-D*X0^2=9 gives DIRECT bounds

    Z^2<=R*Z0^2/9,       X^2<=R*X0^2/9.                  (1)

No relative-unit completeness theorem is needed. At the other embedding,

    Z'^2+2(s+1)X'^2=R'.

For all six rows,0<R<540 and0<R'<=144. Since3/2<s<7/4, equation(1) gives

    0<Z<200,       0<X<164,       |Z'|<=12,       |X'|<6.

Consequently Z=A+B*s and X=C+E*s lie in these enlarged INTEGER boxes:

    -5<=A<=105,  -3<=B<=70,
    -2<=C<=84,   -1<=E<=56.

Only8,214 and5,046 raw pairs occur in the two boxes. Indexing X^2 in a
dictionary and matching `(Z^2-R)/D=(1+s)(Z^2-R)/4` covers every possible
first simple endpoint without a four-dimensional product, a height cutoff,
or a unit-orbit search. Exact squared conjugate and interval tests may further
reduce the boxes but may not replace their proof of completeness.

For reference, a relative norm-one unit also exists explicitly:

    (11+6s)+(9+5s)*sqrt(D),
    (11+6s)^2-D*(9+5s)^2=1.

It is unnecessary for this interval. No assertion that it is fundamental is
made.

## First sector and broader scope

At a minimum candidate delta, all simultaneous positive mirrors must be
considered. A candidate simple root must make every other simultaneous root
a nonnegative linear combination of it and the seven retained A7 roots.
With time coordinates T,T_other, this is checked by the seven first-eight
prefix inequalities

    T_other*prefix(root)-T*prefix(other)>=0.

There is no ninth retained ambient rho root in this intrinsic test. Its
orthogonality has already been encoded by M. A projective endpoint can be
stored without divisions as

    t_point=9Z,
    z1=...=z8=mZ+(s-1)T,
    z9=mZ-8(s-1)T,

whose square is exactly `-9(s-1)R<0`. This supplies a complete finite
first-edge computation once implemented and replayed; it is not itself an
execution result.

Independently, the displayed E8 roots have unimodular Gram, so M splits
integrally as E8 plus a rank-one complement. Its norm is1-s up to an O-unit
square. The coordinating/resource agents are developing the explicit
splitting and global reflection-chamber analysis. This does not depend on
whether d proves to be the first wall and may supersede further edge walking.
No compact9 example, complete cone, or nonreflectivity theorem follows merely
from this interval. Any future compactness certificate must first exclude
lineality and bound a future affine cone section, then check ALL extreme rays.
