# Root review: the complete square-root-of-three diagonal-edge reduction

20 October 2026. Written mathematics only. The worker and its output
require separate review; no enumeration or execution is asserted here.

Source:
`compact_h9_sqrt3_diagonal_construction_route_20261020.md`, SHA256
`1ac6604b8760f5406bfc679fe362b8274c1542165c6e32136cf4f84d85b6b96b`.

PASS for a complete finite first-edge reduction, conditional on the
standard finite Coxeter classification explicitly used in the source.
This is not a completed compact H9 chamber or a dimension bound.

## Norm inventory

For O=Z[sqrt3], rounding the two coefficients makes the absolute norm
of a Euclidean remainder at most3/4. Thus primitive coordinates have
a Bezout relation. Reflection integrality first gives
q | 2sqrt3*a and q | 2y_i, because a linear functional taking the
primitive root to1 recovers each reflection coefficient. The Bezout
relation then gives q | 2sqrt3. The two prime factors are 1+sqrt3
and sqrt3, with exponents at most2 and1. Total positivity forces
even parity of their exponent sum. The elementary unit argument and
unit-square normalization give the six printed representatives.

These are six norm representatives modulo squares of UNITS. They
are not six distinct square classes in K*: for instance
2epsilon=(1+sqrt3)^2. The source uses the correct unit-normalized
inventory and also keeps primitivity/integrality checks explicit.

The omitted roots really are integral reflectors. Therefore the older
norm1/norm2 screen cannot be a completion or nonreflectivity proof for
this full host. The full zero stabilizer is finite, contains the signed
permutations, and cannot properly enlarge B9 by the finite rank-nine
reflection classification. This justifies the actual B9 starting cone.

## All possible finite endpoints

A first finite endpoint of the released A8 edge has spherical rank9
corner. The possible extensions are A8+A1, A9, D9 and B9. In the fixed
coordinate sector, contacts are single positive gaps at an endpoint
or next-to-endpoint of the A8 chain. For a gap g after position k,

    sum y_i^2 = (z^2+k(9-k)g^2)/9.

Together with X=3a this gives z^2-sqrt3*X^2=R exactly. Checking which
of q/2 or q has an integral square root recovers precisely the printed
ten equations and both placements of each nonisolated profile.
The congruences z-kg in9O and X in3O are needed to recover coordinates.
Relations with the ninth zero wall and all other earlier actual walls
must still be checked. A raw solution of the norm equation is not
automatically a simple reflector.

## Upper reflection and bounded orbit search

U=(2+sqrt3)+(1+sqrt3)*sqrt(sqrt3) has relative norm1. Its ninth power
has X in3O and Z=-1 modulo9O. The root a=X/3,
b=(Z+1)/9 satisfies 9b^2-sqrt3*a^2=2b. Substitution in the reflection
formula gives exactly t'=Zt-a*sum(y_i) and
y_i'=y_i+sqrt3*a*t-b*sum(y_j), all integral. Its diagonal intersection
is finite, with lambda^2=(sqrt3/9)*(Z-1)/(Z+1). This bounds the first
collision without assuming this large reflector is itself simple.

For every solution with z,X>0, w=z+sqrt(sqrt3)*X exceeds sqrt(R).
Multiplication by a suitable nonnegative or negative integer power of
the integral unit U gives a representative in
sqrt(R)<=w<U*sqrt(R), still with integral z,X and X>=0. The resulting
physical and conjugate bounds are finite coefficient boxes. It is
irrelevant whether U generates the full relative unit group: every
solution is reduced in its own U-orbit.

The collision ratio increases with w/sqrt(R). For the upper witness
this ratio is U^(9/2). A reduced representative has ratio in[1,U),
so only powers0 through4 can lie before that witness. Exact comparison
is X*(Z_9+1)<=X_9*z. Negative powers would give X<0 and do not describe
the outgoing positive germ. This checks the finite search bound.

Any ensuing first-edge certificate is local. Compactness of a complete
H9 example needs all extreme rays of an actual finite chamber cone
to be future timelike; neither this endpoint nor a large collection
of finite vertices establishes that global condition.
