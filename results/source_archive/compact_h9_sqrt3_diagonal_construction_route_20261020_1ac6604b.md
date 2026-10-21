# A complete diagonal-edge route for a square-root-of-three H9 host

20 October 2026. Source-only construction research. No enumeration, replay,
scientific process, or compilation was run for this note. No compact or
finite-volume H9 example is claimed.

## The host and the part omitted by the old screen

Put s=sqrt(3), O=Z[s], epsilon=2+s, and

    L = <-s> + <1>^9,
    (t,y)^2 = -s t^2 + sum y_i^2.

At the other real embedding the form is positive definite. In particular it
is anisotropic over Q(s), and a finite-volume chamber of integral reflections
would be compact. This is not a rational quadratic form over Q. The old
quadratic diagonal screen in results/new_example_search_2026_10_13.md included
this coefficient but used norms1 and2 only; its infinite-volume prefix does
not exclude the present complete reflection arrangement.

The ring O is norm-Euclidean: rounding both coefficients gives a remainder
a+bs with |a|,|b|<=1/2, whence |a^2-3b^2|<=3/4<1. Its units are
±epsilon^j. For completeness, a positive unit in [1,epsilon) has conjugate
equal to its reciprocal (there is no norm-minus-one unit, by reduction mod3).
Its s coefficient is (u-u')/(2s), strictly between0 and1, so it is rational
and equals1. Reducing by epsilon proves the assertion. Thus its totally
positive units modulo unit squares are1 and epsilon.

If r is a primitive reflector with square q, reflection integrality and a
coordinate Bezout relation imply q divides2 det(L)= -2s. The prime elements
eta=1+s and s have norms -2 and -3, respectively, and eta^2=2epsilon.
The possible prime exponents in q are0<=v_eta(q)<=2 and0<=v_s(q)<=1.
Total positivity forces their sum even. Consequently a complete necessary
norm pool, up to multiplication of r by units, is

    Q = {1, epsilon, 2, 2epsilon, 3+s, (3+s)epsilon}.

This is a safe pool, not a claim that every class occurs. Integral reflection
and primitivity tests remain necessary for every candidate.

Two genuinely omitted reflector directions are

    r_i=(1; s e_i),             r_i^2=3-s,
    r_ij=(1; e_i+e_j),          r_ij^2=2-s.

For r_i, both 2s/(3-s)=1+s and 2s/(3-s) are integral coefficients of the
reflection matrix. For r_ij its square2-s is a unit, so integrality is
immediate. The two squares belong to the last and second pool classes,
respectively, after rescaling by units. They cannot be omitted merely because
their primitive squares were not1 or2 in the earlier screen.

The complete zero-height chamber is B9. Indeed the integral spatial isometry
group is finite: it is bounded at both embeddings. Its reflection subgroup
contains the signed permutations from e_i and e_i±e_j. The finite reflection
classification in rank9 leaves no irreducible enlargement of B9. This proves
the zero chamber without an assumed root-length census.

Use its simple outward roots

    n_i=-e_i+e_(i+1) (1<=i<=8),     n_9=-e_9.

The positive diagonal ray issuing from H=(1;0^9) after releasing n_9 is

    p(lambda)=(1;lambda,...,lambda),  0<=lambda<sqrt(s)/3.

The other eight walls give A8. In this sector the omitted roots r_1 and
r_12 close the coordinate-axis and two-equal-coordinate germs at finite
points: their equations are y_1=t and y_1+y_2=s t. They do not close the
positive diagonal germ. They are mentioned only to show that the omitted
norm classes change actual low-distance geometry.

## Exactly ten relative norm equations suffice at the diagonal

Let r=(a;y_1,...,y_9) be a primitive new simple wall at the first finite
endpoint of this germ, rescaled to a square q in Q. The finite corner contains
the eight actual A8 walls. Complete finite Coxeter classification leaves
only A8+A1, A9, B9, or D9. If there is a contact it is a single gap g>0,
at position k, so

    y_1=...=y_k=b+g,   y_(k+1)=...=y_9=b.

The possibilities are:

| Type | q | g | k | R |
|---|---|---|---|---|
| A8+A1 | any q in Q |0| -- |9q|
| A9 |2|1|1 or8|10|
| D9 |2|1|2 or7|4|
| B9 |1|1|1 or8|1|
| B9 |2epsilon|1+s|1 or8|2epsilon|

To check the length exclusions: a3-bond requires g^2=q/2. Only q=2 gives
an integral g. Although epsilon/2=((1+s)/2)^2, its square root is not in O.
A4-bond requires g^2=q; the two integral possibilities are q=1 and
q=2epsilon=(1+s)^2. The remaining candidates are excluded by their rational
norms or the fact that epsilon and2 are not squares in Q(s).

Set z=sum y_i and X=3a. The norm equation is precisely

    z^2-s X^2 = 9q-k(9-k)g^2 = R.                       (1)

For the isolated profile z=9b. For a nonzero gap z=9b+kg. Thus the required
congruences are X in3O and z-kg in9O (or z in9O when g=0). The ten values
of R in the table give a complete finite list of equations, with two possible
placements for each nonisolated row. This is a necessary endpoint inventory,
not yet an enumerated list of roots or a first-mirror certificate.

## A finite upper endpoint, with no search

Write alpha=sqrt(s), and work in the auxiliary quartic field Q(alpha). Put

    U=(2+s)+(1+s)alpha.

Its relative norm to Q(s) is1, since

    (2+s)^2-s(1+s)^2=1.

Write U^9=Z+alpha X. These are integral Z,X in O, with Z>1 and X>0 at the
physical embedding. Direct modular identities give

    U^3=(98+57s)+(75+43s)alpha,
    Z=-1 (mod9O),       X=0 (mod3O).

For example reducing the displayed cube modulo9 and cubing again gives
Z=8 and X=3s modulo9. Define

    a=X/3,        b=(Z+1)/9,        rho=(a;b,...,b).

Then rho belongs to L and

    rho^2=2b>0.

Its reflection is the explicit integral map

    t'   = Z t-a sum y_i,
    y_i' = y_i+s a t-b sum y_j.                         (2)

The relative norm identity verifies this formula directly; it is an
involution preserving the form. Hence rho gives an actual integral reflector,
whether or not this displayed vector is primitive. Its hyperplane meets the
diagonal at

    lambda_* = s a/(Z+1),
    lambda_*^2 = (s/9)(Z-1)/(Z+1) < s/9.                (3)

Thus the complete integral reflection chamber has a first diagonal collision
at or before a known finite point. Local finiteness makes the first collision
exist. Equation(2) is a useful upper-wall witness; it does not say that rho
itself is the first simple wall, nor that its displayed square is in the
primitive norm pool before dividing a common factor.

## A complete bounded next computation, rather than a height search

No class-number computation, maximal-order basis, or assertion that U is a
fundamental unit is needed. For each positive R in the table, let

    w=z+alpha X,        w_bar=z-alpha X=R/w.

Every solution with z,X>0 can be multiplied by an integer power of U so that

    sqrt(R)<=w<U sqrt(R).

Multiplication by U and U^-1 preserves the O-pair(z,X), equation(1), and
exhaustively moves the entire orbit. In this strip,

    0<z<(2+s)sqrt(R),       0<=X<(1+s)sqrt(R).

At the nonphysical embedding, (1) instead says

    z'^2+s X'^2=R',

so |z'|<=sqrt(R') and |X'|<=sqrt(R'/s). For z=A+Bs these physical and
conjugate intervals give explicit finite integer bounds on A and B; the same
holds for X. Enlarging these bounds rationally is harmless, provided final
tests are exact. Enumerating these ten small boxes is therefore COMPLETE for
orbit representatives, even if the full relative unit group is larger.

In fact only FIVE powers per representative are needed. For a positive
solution of(1), collision position is strictly increasing with w/sqrt(R).
For the explicit upper wall rho this ratio is

    (U^9+1)/sqrt(2(Z+1)) = U^(9/2),

because 2(Z+1)=U^9+U^-9+2. The strip starts at ratio1, so powers U^j with
j>=5 already lie beyond the upper collision. It suffices to test j=0,...,4,
using the exact comparison X(Z+1)<=X_9 z rather than adjoining sqrt(U).
At every retained power
one checks the profile congruences, q, positive orientation, primitivity,
reflection integrality, outgoing direction, and all earlier wall relations.
The smallest permitted collision must be compared with simultaneous incident
mirrors and the actual simple finite sector, rather than selecting an
arbitrary root at that point. This produces a complete first-edge certificate
if implemented. It has not been run.

Even a successful first edge would be local progress only. A compact example
requires a finite acute full root list whose ENTIRE projective cone has only
future timelike extreme rays. Counts of finite corners, many closed germs,
or absence of a sampled escape do not suffice. If a finite list contains a
spacelike generator it has infinite volume; if it retains a null generator
it cannot be the compact target. No such global certification is presently
available for this host.
