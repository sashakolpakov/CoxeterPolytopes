# The certified scalar wall and its small eight-dimensional facet lattice

20 October 2026. Paper derivation after the first-edge generation. No new
scientific execution, replay, Lean, or build was run for this note.

The generation of
`scripts/check_compact_h9_sqrt3_diagonal_first_edge_20261020.py` found exactly
one first-edge candidate. Its certificate is
`results/compact_h9_sqrt3_diagonal_first_edge_20261020.json`.
The separately scheduled replay has not started: the resource gate deferred
three attempts before launching a child. This note does not claim a replay.

## The selected wall is the analytic upper wall

Write s=sqrt3, epsilon=2+s, eta=1+s, and alpha=sqrt(s). The selected normal is

    (6230+3597s; (2733+1578s)^9),  square9+5s.

It is epsilon times the simpler primitive normal

    rho=(a;b^9),
    a=1669+964s,   c=423+244s,   b=s*c=732+423s,
    q=rho^2=3-s.

Indeed9b^2-s*a^2=3-s. Its rational norm is6, so this normal cannot have a
nonunit common factor. The identity

    U=((1+alpha)^2)/(s-1),
    (1+alpha)U^4=(1669+964s)+(1269+732s)alpha

identifies its hyperplane with the U^9 midpoint reflection in the route note.
Thus the run found that analytic upper mirror itself, rather than an earlier
mirror. The endpoint is the future finite point

    V=(1;lambda^9),       lambda=s*a/(9b),

with incident roots the eight A8 differences and rho. This is a complete
first-edge statement, not a finite-volume statement.

## A useful refinement of the six-norm necessary pool

The route used the safe primitive norm representatives

    1, epsilon, 2, 2epsilon, 3+s, 9+5s.

These are normalized modulo squares of O-units, where O=Z[s]. They are not
six square classes of K*. The class3+s is actually empty in this literal
lattice L=<-s>+<1>^9. If r=(A;y) had that primitive norm, reflection
integrality at the prime(s) would force every y_i divisible bys. Dividing
the norm equation bys would give

    -A^2+s*sum z_i^2=1+s.

Modulo(s), whose residue field is F3, this says A^2=-1, impossible.
The other norm6 representative9+5s is unit-square equivalent to3-s.

The native finite labels over K=Q(s) are2,3,4,6,12, since their squared
cosines must belong to K. A norm6 root cannot meet a unit-norm or norm2-class
root obliquely: for each nonright native label, the rational norm of
`q*q_other*cos^2(pi/m)` is not a rational square. This quantity would have to
be the square of a K-valued raw Gram product. The same argument allows only
label4 between the TWO norm6 unit classes; the other class was just excluded.
Within the remaining class, labels3 and6 would require q/2 or s*q/2 as an
integral Gram product, which they are not; labels4 and12 require elements
outside K. Thus the only possible intersection between remaining norm6
mirrors would be orthogonal.

In fact they cannot be orthogonal either. Reflection integrality forces
every spatial coordinate of any primitive norm6 root to be divisible bys;
its time coordinate is not divisible bys, since otherwise its norm would
be divisible bys^2. For two such roots their product is

    -s*A*A_other modulo(s^2),

which has valuation one and cannot vanish. Distinct norm6 mirrors are
therefore disjoint. Parallel mirrors are also impossible in this admissible
anisotropic K-form: a degenerate two-root Gram would remain degenerate at
the positive-definite embedding.

Consequently EVERY neighbor of rho is orthogonal to rho, and has primitive
norm in {1,epsilon,2,2epsilon}. The provisional possibility of a4-neighbor
does not survive the local prime test. This is a statement about all
integral reflectors of this literal lattice, not only the saved prefix.

## The entire orthogonal lattice has small coefficients

For (A;y) in L, orthogonality to rho is

    a*A=c*sum y_i.

The elements a,c are coprime: any common prime would divide every coordinate
of rho and have its square divide3-s. Hence A=c*t and sum y_i=a*t for a
unique t in O. Put

    m=4+s,       h=(a-m)/9=185+107s,       y_i=z_i+h*t.

This gives an integral bijection onto

    M={ (t,z) in O^10 : sum z_i=(4+s)t },
    Q_M=sum z_i^2-(2+s)t^2.                              (1)

The identity follows from a^2-9s*c^2=1-s and
`(1-s-m^2)/9=-(2+s)`. In the free integral basis t,z1,...,z8,
with z9=m*t-sum_(i<=8)z_i, the Gram matrix is

    [[17+7s, -(4+s)*1^T],
     [-(4+s)*1, I8+J8]],

and its determinant is1-s. For example det(I8+J8)=9 and
`1^T(I8+J8)^-1*1=8/9`, giving that determinant immediately. Its physical
signature is(8,1); its conjugate is positive definite.

A primitive integral reflector of M has norm dividing2(1-s). The same
Euclidean/Bezout and unit argument as before therefore leaves the necessary
norm pool {1,epsilon,2,2epsilon}. The unit norms are not represented:
reducing the displayed integral Gram modulo2 gives

    Q_M=(1+s)*t^2 modulo2.

For t=A+B*s, t^2=A+B modulo2, because s^2=3=1 there. The possible residues
are therefore0 or1+s, whereas the unit norm representatives1 and epsilon
have residues1 ands. Thus the COMPLETE primitive M-reflector pool reduces to

    {2, 2epsilon}.                                      (2)

Every such reflector of M is also an integral reflector of L: its norm is
an associate of2, so the ambient reflection coefficients are integral.
Primitivity is compatible because M=L intersect rho-perp is saturated.
Conversely every ambient mirror orthogonal to rho preserves M. Together
with the preceding no-oblique-neighbor result, this identifies the rho facet
of the complete L chamber with a full integral-reflection chamber of M.
In particular a compact L chamber would force this M chamber to be compact.
This supplies a genuine smaller global stopping test: nonreflectivity of M
would exclude the L construction. No such nonreflectivity proof or known
lattice identification is asserted here.

## The next genuinely different edge

Release n1 at V while retaining rho and n2,...,n8. Inside (1), this is an A7
edge, with the eight coordinates z2,...,z9 equal. It is different from the
original A8-only edge, and its arithmetic extension field is different.

Finite classification and(2) give precisely these next endpoint profiles:

| Intrinsic type | New norm q | gap g among z2,...,z9 | position k |
|---|---|---|---|
| A7+A1 |2 or2epsilon|0|--|
| A8 |2|1|1 or7|
| D8 |2|1|2 or6|
| E8 |2|1|3 or5|
| B8 |2epsilon|eta|1 or7|

A3-bond requires g^2=q/2, so only q2 supplies its integral gap1. A4-bond
requires g^2=q, so only q2epsilon supplies the integral gapeta. The branching
position4 produces an affine diagram, not another finite case. There is no
bridge to the retained ambient rho component: all these walls are right to it.

If the eight coordinates have k entries p+g and8-k entries p, set
S8=8p+k*g and Z=9S8-8m*t. The remaining coordinate is z1=m*t-S8.
Equation(1) becomes the exact necessary relative norm equation

    Z^2-8(s-1)t^2 = 72q-9k(8-k)g^2.                     (3)

The six profile rows give right sides144,144epsilon,81,36,9,18epsilon.
No first upper mirror or complete orbit interval for(3) has yet been proved.
That is the next mathematical obligation before another worker is justified.
No immediate infinite symmetry follows from the first edge: its midpoint
map is a wall reflection, not a chamber stabilizer, and the A8 reversal
fixing the apex and rho has spatial matrix -reverse+(2/9)J9, not an integral
lattice map. A later symmetry must be established independently.

For a future global compactness certificate one must first exclude lineality
and prove a bounded future affine cone section. Only then does a COMPLETE
extreme-ray census consisting entirely of future timelike rays certify a
compact chamber. A long list of finite corners or closed edges is not such
a certificate. No global compactness claim is made in this batch.
