# Norm79: a finite A7+c edge upper and nine complete interval models

21 October2026. Manual source/proof continuation, pending full review.
No enumeration, scientific/symbolic import, optimizer, worker, Lean or
TeX process ran. No frozen source, paper, guide or receipt was edited.
This proves an integral finite upper and a finite necessary domain; the
minimum and actual new endpoint have NOT been enumerated or identified.
It supplies neither an H9 example nor a nonreflectivity theorem.

## 1. Literal setting and inherited actual geometry

Use O=Z[phi], phi^2=phi+1, the native E8 basis with diagonal2 and
products-1 on02,13,23,34,45,56,67, and

    a=3phi^2, b=13+18phi, D=phi*b=18+31phi,
    L=E8(O) orthogonal-sum [[-2phi,a],[a,2]].               (1)

The accepted source33f6 proves signatures(9,1)/(10,0), the exact
primitive mirror alphabet2,2b, an ACTUAL compact Bugaenko base, and
the ACTUAL corner

    p1=(2,-a;0), p1^2=-2D,
    full incident germ E8+A1(c), c=(0,-1;0), c^2=2.

No claim about another integral module or the whole K-similarity class
is imported. In particular (1), not a merely similar diagonal form,
is the module used in all integrality calculations below.

Release -e1 and retain A7 on the path0,2,3,4,5,6,7 together with c.
Let Q=omega1=(5,8,10,15,12,9,6,3), Q^2=8. The outgoing edge is

    p(tau)=p1+tau*Q, tau>0, light value tau^2=D/4.        (2)

The primitive norm2 and norm2b profiles orthogonal to ALL retained
walls are impossible at every height. Their equations would be
4Y^2-DX^2=1 and4bV^2-phiX^2=1. Modulo4, nonzero unit squares are
{1,1+phi,2+3phi}, and D=2+3phi mod4. Multiplication by D, or by
phi, never gives -1 in that set. An even X also fails either equation.
This is the accepted all-height congruence from33f6, not a height cutoff.

## 2. A new all-height obstruction: no c-only native5 opposite

A norm2 opposite orthogonal to A7 and joined only to c by5 has

    h=(U/2)p1-(phi/2)c+YQ,
    16Y^2-DU^2=3-phi, phi-aU in2O, U,Y in O.             (3)

The parity means U=phi^-1+2w for some w in O. As
D*phi^-2=5+13phi, substituting into(3) and dividing by4 gives

    D*(w^2+phi^-1*w)=4Y^2-2-3phi.

Modulo4, D^-1=phi^2 and
(-2-3phi)*phi^2=-5-8phi=3 mod4. Reducing once more to
O/2=F4 therefore requires

    w^2+phi^2*w=1 in F4.                                (4)

The map on the left has image {0,phi}: at w=0,phi^2 it is0,
and at w=1,phi it is phi. Since1 is neither, (4) is impossible.
Thus the c-only native5 profile is excluded at ALL heights in this
literal module. No relative-unit completeness or numerical search is
used in this extra elimination.

## 3. An explicit finite integral upper, not a guessed first wall

The elementary identity from the initial construction is

    V=5+9phi, X=27+45phi, bV^2-phiX^2=1.

Put

    Z0=2bV^2-1=8911+14418phi,
    W0=2VX=1080+1746phi.

Direct expansion gives Z0^2-DW0^2=1. Multiplication by
sqrt5=2phi-1, an element of O, then gives

    U*=sqrt5*W0/2=1206+1953phi,
    Z*=sqrt5*Z0=19925+32240phi,
    Z*^2-4D*U*^2=5,
    Y*=(Z*-5)/8=2490+4030phi.                            (5)

Consequently the A9-bridge vector

    h*=(U*/2)p1-c/2+omega0+Y*Q                          (6)

is an integral norm2 reflector. In the ORIGINAL binary basis its
coordinates are

    (1206+1953phi, -4738-7668phi),

since aU*=9477+15336phi; its spatial coordinates are integral as
well. Norm2 in the even lattice makes it primitive, and the reflection
covector is integral automatically. Its retained products are-1 at
old0 and c and0 at the other A7 roots. The retained-plus-h* formal
diagram is A9. Firstness and actual-facet status are NOT inferred
from that positive definite diagram.

Indeed the argument for an upper uses the integral mirror itself:
(p1,h*)=-DU*<0 and (Q,h*)=Z*>0, so it intersects(2) at

    tau*=D U*/Z*,
    D-4tau*^2=5D/Z*^2>0.                               (7)

This is a strictly finite point before the light cone. The full integral
arrangement is locally finite, the initial E8+c germ is actual, and
h* is a mirror oriented negative at p1. Thus the actual A7+c edge has
a first finite cut at or before tau*, even if h* proves redundant or
an earlier different wall is present. No all-wall feasibility or
chamber-stabilizer property of h* is being assumed.

The norm identity Z0^2-DW0^2=1 is a scalar relative-Pell identity.
It does NOT by itself give an integral isometry of L, an actual-germ
map, an infinite chamber stabilizer or a nonreflectivity proof.

## 4. Complete remaining first-SIMPLE inventory

Every norm2b first wall must be orthogonal to the norm2 retained roots:
any nonzero finite angular invariant would have a prime79 factor.
Section1 already excludes that independent profile. For norm2, its
pairings with retained norm2 roots lie in O. The golden angle list
leaves only native3/5: the raw products for4,6 and the other listed
finite labels require square roots outside K. A connected finite
rank8/9 diagram containing A7 cannot have a5-bond, since those finite
components are only I2/H3/H4. Thus the connected extension uses3 and
is of A,D,E type; this also explains why B8/B9 are not omitted cases.
A single attachment gives A8 at an
endpoint, D8 at a next-to-endpoint, or E8 at a third-from-endpoint.
A wall connecting both A7 and c gives A9 only, at either path end.
An internal bridge has three arms of which one has length2; the
possible rank9 D-type arm lengths cannot occur. Other finite diagrams
of rank9 and any branch with a5-bond are unavailable.

Together with sole c bonds3/5 and the two independent profiles, this is
the full twelve-model list. Section1 removes both independent profiles,
giving the ten remaining models of33f6. Section2 removes the c-only5
case, leaving EXACTLY these NINE necessary
model types (not nine actual walls):

| Type | Controls and equation | Extra original-module conditions |
|---|---|---|
| c-only3 | h=(U/2)p1-c/2+YQ; 16Y^2-DU^2=3 | U=phi mod2; U>0,Y>0 |
| A8, j0/7 | h=Xp1+omega_j+YQ; Z^2-16DX^2=9 | Z=8Y+I_j; X>0,Y>=0 |
| D8, j2/6 | same normal; Z^2-16DX^2=4 | same conditions |
| E8, j3/5 | same normal; Z^2-16DX^2=1 | same conditions |
| A9 bridge, j0/7 | h=(U/2)p1-c/2+omega_j+YQ; Z^2-4DU^2=5 | Z=8Y+I_j; U=phi mod2; U>0,Y>=0 |

Here the exact inverse-E8 data are

    j                 0   2   3   5   6   7
    I_j=(Q,omega_j)   5  10  15   9   6   3
    omega_j^2         4  14  30  12   6   2
    I_j^2-8omega_j^2+16
                      9   4   1   1   4   9.

E8 unimodularity makes all listed controls elements of O. The parity
U=phi mod2 is precisely1-aU in2O, not a copied unit-family condition.
Y>=0 in a connected row is actual acuteness with the released -e1.
For an actual opposite, its own product also has to be native: with
norm2 on both sides, Y is0,1,phi or>2. The value2 would be parallel
and is impossible in this anisotropic K-space. Satisfying one norm
equation is not sufficient to make the normal an actual wall.

## 5. Exact complete physical and conjugate boxes

For a c-only row, an extension row of parameter R, or a bridge,
their respective collision parameters are

    tau=DU/(8Y),       tau=2DX/Z,       tau=DU/Z.

In all cases divide tau^2 by D/4 and use the displayed norm equation.
The resulting increasing fractions are DU^2/(DU^2+3),
16DX^2/(16DX^2+R), and4DU^2/(4DU^2+5). Comparison with(7)
gives the EXACT physical domain before or at tau*:

    c-only3: U^2 <= (12/5)*U*^2,
    extension: X^2 <= (R/20)*U*^2,
    bridge: U^2 <= U*^2.                                (8)

The other embedding, D'<0, gives corresponding bounds

    |U'|<=sqrt3/sqrt(-D'),
    |X'|<=sqrt(R)/(4sqrt(-D')),
    |U'|<=sqrt5/(2sqrt(-D')).                            (9)

These hold for EVERY possible first-simple wall before the upper,
without a height search, approximate root ordering or saved survivor
selection. All comparisons in (8)/(9) can be made by positive-square
clearing in Q(phi). Endpoints and ties must be included.

For simple conservative integer loops, 1.618<phi<1.6181 gives
U*<4367, -D'=31phi-49>1.158, sqrt5>2.236 and the following table.

| Row | Physical bound | Conjugate radius | phi-coefficient 0<=B<= | Integer A trials per B |
|---|---:|---:|---:|---:|
| c-only3, one model | U<6800 | 1.61 |3050|4|
| A8, two models | X<3000 | .70 |1350|2|
| D8, two models | X<2000 | .465 |900|2|
| E8, two models | X<1000 | .233 |450|2|
| A9 bridge, two models | U<4400 | 1.04 |1970|4|

For example .70^2*16*1.158>9, .465^2*16*1.158>4,
.233^2*16*1.158>1, 1.61^2*1.158>3 and
1.04^2*4*1.158>5. Physical bounds follow from(8), U*<4367,
sqrt(12/5)<1.55 and sqrt5>2.236. Thus every radius in the table is
strictly conservative. For t=A+Bphi>0, |t'|<r, the identity
B=(t-t')/sqrt5 implies B>=0 (all r<sqrt5) and the displayed maxima;
the chosen integer maxima are deliberately loose.

For fixed B set n=floor(B*(phi-1)). If r<1, A=n,n+1 suffice;
otherwise the four values n-1,n,n+1,n+2 suffice because r<=1.61<2.
The exact floor is (isqrt(5B^2)-B)//2. All strict strip, exact(8)/(9),
parity and original norm conditions are then applied. The total raw
scalar-trial bound, before any parity or strip rejection, is

    4*3051 + 2*2*1351 + 2*2*901 + 2*2*451 + 2*4*1971
       =38784.                                          (10)

Thus a future complete check can be small and bounded. This is a
proved domain size, not a completed runtime or a predicted resource PASS.

## 6. What one future exact certificate must establish

For each control in the complete loops, a checker must solve the field
square exactly, reconstruct the native normal, check its norm and
reflection covector, and enforce all displayed congruences/orientations.
In the c-only row one may square-test DU^2+3 and demand its positive
root be divisible by4; in other rows test R+16DX^2 or5+4DU^2 and
demand Z-I_j in8O. A square rational norm is only a necessary filter:
both possible norm signs, trace/coefficient square conditions, parity
and exact resquaring are needed for a COMPLETE O-square test.

Any actual first wall must also be acute against EVERY known actual
wall in33f6, with legitimate native labels when their hyperplanes
intersect. There is no permission to replace those by finite Gram
positivity alone. Applying such known-wall tests cannot remove an
actual first simple wall. The explicitly constructed h* is an upper
mirror, so its finite cut remains valid even if these stronger actual
facet tests reject it.

The output must retain the exact earliest collision and all ties;
independent/profile-duplicate vectors must not be merged by row order.
To recognize the whole actual germ, exhibit the correct positive
definite simple extension and certify that tied redundant candidates
are in its nonnegative simple cone. Completeness of the first-simple
list, local finiteness and this cone test are the geometric bridge
from the finite result to an actual endpoint. If earlier profiles
exist, h* need not be a chamber wall at all.

No such interval execution or all-tie calculation has been performed
here. The precise completed result is the all-height c-only5 obstruction,
the explicit finite upper, existence of an actual finite next endpoint,
and the complete nine-model domain(8)--(10). Root alone may authorize
ONE source-reviewed, guarded exact census after unchanged disk/RAM,
lock and child-resource gates. No job is authorized or claimed by this
note. Completion of this edge still would not by itself certify finite
volume of the full chamber, an infinite-order chamber period, a new
compact H9 example or a K-form-class obstruction.

## Bound mathematical inputs

| Input | SHA256 |
|---|---|
| results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md | `33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3` |
| results/review_compact_h9_prime79_full_base_and_actual_edges_full_outside_author_20261021.md | `0770adfa417f586d04fb4cc0fb20b65ff29fe0fd33f4f01dfc56f126165405e9` |
| data/known_examples/bugaenko_compact_h8.json | `c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c` |
| data/known_examples/bugaenko_vinberg_transcript.json | `c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a` |

The accepted source's golden ring/primitive-normal facts, finite spherical
classification, complete arrangement-angle list, local finiteness and
finite-corner simplicity retain their explicit scopes. Its root/outside
review status is not evidence for the NEW calculations in this note.
