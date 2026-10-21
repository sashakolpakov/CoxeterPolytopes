# All-height pencil obstructions and finite endpoint profiles in H19

19 October 2026. Independent paper review of the coordinating agent's
local arithmetic arguments, with a further elementary finite-profile
reduction. No computation, search, optimizer, Lean process or build ran.
The constructive and escape reports remain unchanged.

**Result.** The omega2 pencil has no canonical short or long lattice
root, at any height; the omega3 pencil has no canonical long root.
Canonical short and long roots cannot have a nonright finite Coxeter
angle. Consequently only short, nonorthogonal roots can close the
initial omega2 edge. There are exactly 51 remaining labeled finite
white-contact profiles, before the other chamber inequalities and
integral norm equations are imposed. This does not prove that any of
those profiles is realizable, or that none is realizable.

## 1. Normalization and completeness of the orthogonal pencils

Keep K=Q(sqrt(5)), O=Z[phi], phi²=phi+1, b=4+phi,
D=1+5phi=phi*b, and the native Gram

    G=[[-2phi,phi],[phi,2]] + C_E8.

We use the complete primitive reflection norm alphabet {2,2b}, up to
unit rescaling, proved in the unchanged mixed-norm seed report. It is
an explicit dependency. Multiplication by a unit retains integrality
and every orthogonality, so canonical normalization loses no direction.

Let omega_j=C_E8^{-1}e_j. A vector orthogonal to old0 and to all negative
E8 simples except -e_j has exactly the form

    q=(2u,-phi*u,t*omega_j),
    Gq=(-D*u,0;t*e_j),
    q²=2(k_j*t²-D*u²),  k_j=omega_j²/2.             (1)

Here u,t are integral: old0-orthogonality gives u=-phi^{-1}y in O,
and the native integral covector gives t in O. In particular there is
no hidden denominator at 2. For a canonical long reflector, full-lattice
reflection integrality further implies Gq in b O^10 and therefore
t=b*tau with tau in O. Thus its equation, after division by b, is

    k_j*b*tau²-phi*u²=1.                            (2)

The native fundamental weights in question are

    omega2=(10,20,30,24,18,12,6,15), omega2²=30,
    omega3=(8,16,24,20,15,10,5,12), omega3²=20.

Their displayed Cartan images are e2,e3; thus k2=15 and k3=10.

## 2. The entire omega2 pencil is obstructed modulo 3

The polynomial T²-T-1 is irreducible over F3. Hence O/(3)=F9, and
its norm to F3 is the reduction of the rational field norm. The latter
satisfies

    N(D)=-19=2 mod3,   N(phi)=-1=2 mod3.

For a canonical short root, (1) would give

    15t²-Du²=1,
    u²=-D^{-1} in F9.

The right side has F9/F3 norm 2, while a nonzero square has norm 1.
This is impossible. Neither side can vanish.

For a canonical long root, (2) instead gives

    15b*tau²-phi*u²=1,
    u²=-phi^{-1} in F9.

Its right side again has norm 2 and is not a square. Thus no short or
long lattice reflection direction exists in this entire pencil, not
merely in a bounded height range.

## 3. The omega3 long pencil is obstructed at sqrt(5)

The ramified prime ideal (2phi-1) has residue field F5, with phi=3,
b=2 and D=1. Reducing the canonical long equation

    10b*tau²-phi*u²=1

gives -3u²=1 and therefore u²=3 in F5, a nonsquare. Equivalently one
can reduce 10t²-Du²=b directly to -u²=2. No long direction in this
pencil survives at any height. This argument does not exclude its
canonical short equation 10t²-Du²=1.

## 4. No nonright finite angle between the two norm classes

For roots a²=2 and l²=2b with finite angle pi/m, m>=3, their integral
product would satisfy

    (a,l)²=4b*cos²(pi/m).

Since the square is in K, cos(2pi/m) is in K. The elementary real
cyclotomic degree formula first restricts m to

    {3,4,5,6,8,10,12}.

The cases 8 and 12 generate Q(sqrt(2)) and Q(sqrt(3)), respectively,
and are outside K. For the remaining five values the required square
and its rational norm are

|m|4b cos²(pi/m)|Rational norm|
|---:|---:|---:|
|3|b|19|
|4|2b|76|
|5|b*phi²|19|
|6|3b|171|
|10|b*(2+phi)|95|

Every displayed rational norm is a nonsquare. But the rational norm
of a square in K must be a rational square. This excludes all these
angles. Only orthogonality remains at a finite common point. The
same argument even excludes a parallel mixed pair, since N(4b)=304
is a nonsquare; that stronger observation is not needed below.

For two canonical short roots the finite contacts are precisely
0,-1,-phi. The same real-cyclotomic restriction applies; m=4,6
would require sqrt(2),sqrt(3) in K, and m=10 would require
sqrt(2+phi), whose square has nonsquare rational norm 5. This short
contact alphabet is complete, not a cutoff at an arbitrary bond size.

## 5. Consequence for the initial omega2 edge

Its eight incident white normals are old0 and the seven negative E8
simples other than old3=-e2. Their positive-definite Gram type is

    W=A1(old0) + A2(nodes0,1) + A4(nodes3,4,5,6) + A1(node7).

If another wall closes this geodesic edge at a finite point, its normal
together with W has a positive-definite Gram matrix. A canonical long
wall must be orthogonal to every white normal by Section4; it is then
in the omega2 pencil, contradicted by Section2. A short wall orthogonal
to all of W is excluded by the same section. Therefore a closing wall
must be canonical short and must have nonzero finite contact with W.

This is a condition on an actual finite endpoint, not merely on an
abstract attachment. It does not assume that an arbitrary candidate
whose Gram is positive-definite passes the other retained walls.

## 6. Exactly 51 nonzero finite white-contact profiles remain

Here a contact means -1 or -phi with an actual white normal. We can
enumerate every positive-definite attachment using the elementary
Schur complement, without a finite Coxeter diagram classification.
If alpha records the nonnegative absolute contact values, positivity
of the enlarged normalized Gram is exactly

    alpha^T W^{-1} alpha < 2.                       (3)

Within A2 or A4 two distinct nonzero contacts already contribute at
least 2, since all inverse entries are nonnegative. For A2 the value
with both contacts 1 is 2. For A4 the minimum over two distinct contacts
is also 2, attained at its two endpoints. Increasing any contact to
phi increases this value. Thus at most one node per white component
can be contacted.

A single contact 1 has the following contributions to (3):

|White component and position|Contribution|Number of labeled choices|
|---|---:|---:|
|A4 endpoint|4/5|2|
|A4 internal node|6/5|2|
|A2 node|2/3|2|
|Either labeled A1|1/2|2|

A contact phi multiplies its contribution by phi². These values follow
directly from (A_l^{-1})_ij=min(i,j)(l+1-max(i,j))/(l+1).

With only contacts 1 there are:

* 1 profile with no contact;
* 8 with one contacted component;
* 21 with two: A4+A2 gives 8, A4+A1 gives 8, A2+A1 gives 4,
  and the two A1 give 1;
* 12 with three: A4(endpoint)+A2+A1 gives 8,
  A4(endpoint)+both A1 gives 2, and A2+both A1 gives 2.

The other triples fail (3), and all four components have minimum sum
4/5+2/3+1>2. Hence there are exactly 42 zero/1 profiles.

For a contact phi:

* A4 is impossible because (4/5)phi²>2.
* A2 alone is allowed, giving 2 labeled choices. Adding any other
  contact gives at least (2/3)phi²+1/2>2.
* Either A1 alone is allowed, giving 2 choices. It may additionally
  contact the other A1 with value 1 (2 choices), or either A2 node
  with value 1 (4 choices). It cannot contact A4 or two other
  components, and two contacts phi contribute at least phi²>2.

These are exactly 10 additional profiles. For example the closest
allowed boundary case is phi²/2+2/3<2, equivalent to phi<5/3.
All strict allowed and excluded inequalities follow from
3/2<phi<5/3. Thus precisely 52 labeled profiles have positive-definite
white-plus-new Gram. The all-zero profile is already obstructed in
Section2, leaving exactly 51 necessary nonzero profiles.

These 51 are finite contact possibilities only. We have not solved
their integral root equations, tested their compatibility with every
other retained wall, or shown that any can close the actual edge.
There is no whole-lattice nonreflectivity, compactness, finite-volume
completion, or new dimension claim in this reduction.
