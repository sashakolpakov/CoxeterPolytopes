# A double-whisker wall obstruction without a field list

18 October2026. New local research proof, written directly from the Gram
matrix and arithmetic conjugate positivity. No enumeration, optimizer,
TeX build or Lean process was run. Independent review is recorded
separately. This is not a global payment rule or a dimension bound.

## Statement and exact scope

Let P be an n-dimensional compact arithmetic Coxeter chamber of simplest
type, with n≥16. Suppose an actual finite vertex has **whole** diagram
D_n. Denote its fork by a and its two whiskers by b,e. Thus a is joined
to b,e and the first node of the long arm, all with label3. Normalize
outward wall normals to square1.

There cannot be two further actual facet walls c,d such that:

* among all the D_n walls, c meets only e nonorthogonally, with
  c.e=−cos(pi/m), for a finite integer m≥3;
* among all the D_n walls, d meets only b nonorthogonally, with
  d.b=−y and y>1;
* c.d=0.

All unlisted products with the D_n vertex walls are zero. In particular
a.c=a.d=b.c=e.d=0. The statement does **not** treat additional contacts
with the long arm, extra outside components, or a nonzero a.c. It is
the isolated-double-whisker branch of the quadrilateral problem, not
all quadrilaterals through a D_n vertex. The proof uses the usual
arithmetic cyclic-integrality and conjugate-positivity criteria; no
published list of ground fields is needed.

Compactness and finite volume are stronger hypotheses than this local
proof needs. They are retained here to match the active arithmetic lane.

## The corank equation

Write G for the unit-normal D_n Gram. In the standard realization
alpha_i=(e_i−e_(i+1))/sqrt2, with the two spinor end roots, the two
whisker fundamental weights have unnormalized squares n/4 and mutual
product(n−2)/4. Unit-normal dual weights are sqrt2 times these weights.
Consequently

    (G^−1)_(bb)=(G^−1)_(ee)=n/2,
    (G^−1)_(be)=(n−2)/2.

Put u=cos²(pi/m). Taking the Schur complement of G in the Gram of the
n+2 actual normals gives

    S = [[1−nu/2,             −sqrt(u)y(n−2)/2],
         [−sqrt(u)y(n−2)/2,    1−ny²/2        ]].

There are n+2 normals in an(n+1)-dimensional Lorentz space, while G is
invertible. Thus det S=0. Expanding gives

    1−n(u+y²)/2+(n−1)u y²=0.

The denominator (2n−2)u−n is nonzero: if it vanished, the same equation
would require nu−2=0, hence n²=4(n−1), or(n−2)²=0, contrary to n≥16.
Therefore the cyclic product z=4y² is exactly

    z = 4(nu−2)/((2n−2)u−n).                         (1)

Arithmetic cyclic integrality says **z is an algebraic integer**. It
does not say that y, sqrt(u), or a unit-normal Gram entry lies in the
ground field. The squared quantities4u and z are the relevant cyclic
products.

## A trace argument limits the dihedral label

The D_n+c principal Gram has determinant ratio1−nu/2. This is strictly
negative at the physical embedding, since u≥1/4 and n≥16, so those
n+1 normals are independent. At every nonphysical embedding of the
totally real arithmetic ground field, admissibility makes the
corresponding rescaled Gram positive definite. Ratios of principal
determinants cancel the positive scaling of the D_n normals, and the
c-normal scaling has positive sign. Hence

    0 < sigma(u) < 2/n ≤ 1/8                        (2)

for every nonphysical conjugate of u. Each nonidentity embedding of
Q(u) extends to such an embedding of the ground field. This argument
does not apply a field embedding directly to a possibly missing
unit-normal square root.

Let d=[Q(u):Q]=phi(m)/2, for m≥3. If zeta=exp(2pi i/m), then

    u = 1/2+(zeta+zeta^−1)/4,
    Tr_Q(u)/Q(u) = d/2+mu(m)/4 ≥ d/2−1/4.           (3)

The trace of zeta+zeta^−1 in the maximal real cyclotomic field is
mu(m), since summing primitive mth roots gives mu(m) and conjugate
pairs give precisely the real-field conjugates. This establishes(3)
also when m is even; replacing m by2m is unnecessary.

The physical u is strictly less than1. If d≥3, equations(2)–(3) imply

    (4d−2)/8 ≤ Tr(u) < 1+(d−1)/8=(d+7)/8,

which is impossible because4d−2≥d+7. Thus d≤2. The elementary
classification phi(m)≤4 leaves exactly

    m in {3,4,5,6,8,10,12}.                         (4)

For completeness, a prime divisor p of such an m must satisfy p−1≤4,
so only2,3,5 occur; the powers and their products under phi(m)≤4 give
the displayed list, together with1,2 which were excluded by m≥3.

## Every surviving label contradicts the corank equation

For m=3 or4, substitution into(1) gives respectively

    z=2(8−n)/(n+1),       z=−2(n−4).

Both are negative for n≥16, whereas z=4y²>4.

For m=6, equation(1) gives z=6+2/(n−3). This is rational and strictly
between6 and7. A rational algebraic integer must be an integer.

For m=8 the other conjugate u'=(2−sqrt2)/4 is greater than1/8;
for m=10 it is u'=(5−sqrt5)/8, also greater than1/8. Both contradict(2).

For m=5, u=(3+sqrt5)/8 has trace3/4 and norm1/16. Set
D=n²−4n−1. Direct substitution into(1), using both quadratic
conjugates, yields

    Tr(z)=8+20(n−2)/D,
    N(z)=−4+20(4n−13)/D,
    4Tr(z)−N(z)=36+100/D.                           (5)

Here z belongs to Q(sqrt5); its two-field trace and norm must be
integers if z is integral, even if the chamber's ground field is
larger. But D≥191 for n≥16, so100/D lies strictly between0 and1.
Equation(5) is impossible.

Finally, for m=12, u=(2+sqrt3)/4 has trace1 and norm1/16. Equation(1)
gives

    Tr(z)=4+4(7n−15)/(3n²−6n−1).                    (6)

The added fraction is strictly between0 and1: the difference between
its denominator and numerator is3n²−34n+59>0 for n≥16. Thus Tr(z)
is not an integer. This eliminates the last label and proves the claim.

## Relation to the previous local obstruction

[The previous full-D_n isolated-quad report](arithmetic_dn_isolated_quad_obstruction_20261018.md)
assumed a.c≤−1, allowed c to contact the old white walls, and used the
published high-dimensional field list. The present proof has a.c=0 and
restricts c to one whisker contact. These are different branches; neither
report subsumes the other. Together they still leave finite nonzero
a.c, contacts with the long arm, and nonisolated replacement walls open.

No global charge has been reduced from these lemmas alone. The current
compact arithmetic ceiling remains20, and this report is research
support rather than a new standalone manuscript.
