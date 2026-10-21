# Whole A20 with two strong ends: actual endpoint feet

21 October 2026. New root written proof, pending full outside-author
review and source-bound acceptance. No arithmeticity, scientific job,
enumeration, optimizer, Lean or TeX is used. This excludes ONLY the
zero-weak-end case below, not whole A20 or dimension20.

The direct reference input is
`results/compact_large_a_two_end_small_debit_coupling_20261021.md`,
SHA256 `43951a969029090579637fc98f645dac674271caac0ef3337347fe6a0002ee04`,
Sections1--4 at N21: its COMPLETE .61 inventory, near4 removal under
energy<1, literal-helper nonadjacency and all-hole reference. Standard
compact acute/farthest-point geometry and finite spherical classification
retain their usual explicit roles. No field restriction is an input.

**Theorem.** At an actual farthest vertex with WHOLE spherical diagram
A20, indices5 and16 cannot both satisfy

    h_i>sqrt(2i(21-i)/21)-.61.

Consequently at least one of these indices is weak. The separate
two-weak-end proof bf610786 excludes both being weak; their conjunction,
after independent acceptance, would leave exactly one weak end.

## 1. The available reference and the SAME original vector

Use zero-sum coordinates in R21, old unit normals
n_i=(e_i-e_(i+1))/sqrt2, and duals

    Omega_i=sqrt2*(P_i-(i/21)P21),
    D_i=Omega_i^2=2i(21-i)/21,
    h_i=x.Omega_i, ||x||^2<1.

Here x is the original farthest vector throughout; p is its original
future unit vertex. If5 and16 are both strong, the input gives one
concave reference w, after deleting EVERY actual internal weak knot,
with

    E=||w||^2>.981678-.00190125=.97977675,
    x.w>=E.

Its first and last five coordinates are constant and opposite. Its
height at5 is H5=2.150, so its first mean is

    b0=H5/(5sqrt2)<.305.

Write w=sum beta_i Omega_i, beta_i>=0. Its coefficient at5 satisfies

    beta5 >= (H5/5-(H6-H5))/2
           =(.430-.167)/2=.1315.                       (1)

Deleting a weak6 lowers the next slope and only increases this
coefficient. Reflection gives the same fact at16. All subsequent
arguments use this already hole-adjusted reference.

## 2. COMPLETE finite-white list at an actual endpoint

Release old n1. Its actual opposite c has unit norm and

    c=A p-v, A>0, v=gamma Omega1+beta Omega_j,
    gamma=-(c,n1)>=0.

The second term is absent when c is isolated from the retained chain.
There are no foreign components: the vertex diagram is WHOLE A20.
The actual OTHER finite endpoint retains the full A19 on n2,...,n20.
Finite spherical classification gives exactly these possibilities:

* c is isolated;
* one simple3 contact at j2,3,19 or20, beta=1/2;
* one simple4 contact at j2 or20, beta=1/sqrt2.

Indeed a connected completion of the full A19 has rank20 and is A20,
B20 or D20. Two contacts would give a cycle; exceptional finite types
have rank at most8. This proves the list, rather than imposing a bond
cutoff. Gamma is zero, a native cos(pi/m), a parallel1, or dotted>1.
In particular every positive gamma is at least1/2.

For a nonisolated row put

    v1=(D1 gamma+C_(1,j) beta)/sqrt(D1),
    delta=v1^2-A^2>0, Cosh=cosh L=v1/sqrt(delta).

Positivity of the actual replacement gives delta>0, and the actual
edge comparison gives

    x1>=sqrt(20/21)*sqrt((Cosh-1)/(Cosh+1)).           (2)

For a fixed row, Cosh increases with gamma. The following table gives
its COMPLETE minimum calculation; infeasible gamma0 or own3 values
are removed solely because A^2=v^2-1<=0 and A>0 is required.

| White row | delta | Cosh | Necessary lower gamma |
|---|---|---|---|
| simple3 at2 |21/40|(40gamma+19)/21|gamma>=1/2; gamma0 has A^2<0|
| simple3 at3 |1/10|(40gamma+18)/(2sqrt21)|gamma>=0|
| simple3 at19|1/10|(40gamma+2)/(2sqrt21)|gamma>=1/2; gamma0 has A^2<0|
| simple3 at20|21/40|(40gamma+1)/21|gamma>=1/sqrt2; gamma0 has A^2<0 and gamma1/2 has A^2=0|
| simple4 at2 |1/20|(40gamma+19sqrt2)/sqrt42|gamma>=0|
| simple4 at20|1/20|(40gamma+sqrt2)/sqrt42|gamma>=1/2; gamma0 has A^2<0|

All rows except simple3 at20 give Cosh>3/2 even at the displayed lower
gamma. Equation(2) then gives x1^2>4/21>.396^2. In the remaining row

    Cosh>=(20sqrt2+1)/21,
    x1^2>=200(sqrt2-1)/[21(10sqrt2+11)]>.396^2.       (3)

For the last strict inequality, .396=99/250; cross multiplication
reduces it to10441790sqrt2>14764031. Substituting the rational lower
bound sqrt2>707/500 already leaves positive margin660.06. Thus EVERY
nonisolated actual endpoint gives x1>.396. Reflection gives -x21>.396
for a nonisolated opposite endpoint. No claim that every formal row is
realizable is required.

If BOTH endpoints are nonisolated, legally graft their two endpoint
bounds into the disjoint first/last five constant blocks of w. The
input proves concavity at the junction even at graft coordinate1, hence
also at.396. The single modified reference gains more than

    2*(5/4)*(.396-.305)^2=.0207025.

Its energy therefore exceeds1.00047925, a contradiction. It remains
to consider ONE isolated endpoint; reflection lets it be the left one.

## 3. Isolated endpoint: the complete own-contact alternatives

Now v=gamma Omega1 and

    A^2=(40/21)gamma^2-1>0.                          (4)

Own0/3/4 are excluded by(4). Own6 gives the affine3-node chain on
c,n1,n2 with labels6,3, impossible in a compact chamber. Parallel
gamma1 gives the nonzero feasible null vector c+n1, also
impossible. The remaining alternatives are native5, finite native>=7,
or dotted gamma>1. This exhausts the actual native gap.

For each alternative we construct an actual feasible point giving

    x.U>=H.

The negative-core vectors in Sections5--6 are positive combinations of
actual normals, with nonpositive selected products and positive time
coefficient. Section4 instead uses the finite H4 projected POINT
p+(A/delta)(c+w_old); c+w_old alone is spacelike and is not claimed
feasible. That point adds positive normal corrections to the feasible p,
has zero selected products, and is future timelike. Thus every comparison
uses a future feasible point of the ENTIRE actual polytope, including
unknown walls, and the SAME original p and x.

Average U over the first five coordinates, keeping the remaining
coordinates unchanged. In every case this is

    Pi U=(gamma/5)Omega5, Y=U-Pi U, w.Y=0.

The residual Y is supported in the first five coordinates and has sum
zero there. Define Q=||Y||^2 and Delta=H-.430gamma. For any mu>0 with
mu*gamma/5<=beta5, the SINGLE vector W=w+mu Y satisfies

    x.W>=E+mu Delta, ||W||^2=E+mu^2 Q,
    ||x||^2>=2x.W-||W||^2>=E+2mu Delta-mu^2 Q.        (5)

The coefficient payment is made BEFORE using any strong height floor.
It is not a subtraction of an unknown h5 from a previously established
pressure. Orthogonality and the one completed square prove(5).

## 4. Native5

Here .809<gamma=cos(pi/5)<.81. Set

    w_old=(3/2)gamma n1+gamma n2+(gamma/2)n3,
    delta=1-(3/2)gamma^2>0,
    U=v-w_old=(gamma/4)Omega4.

The actual finite H4 projection is
p+(A/delta)(c+w_old). It is future, feasible, and has square
-1-A^2/delta. With B=sqrt(A^2+delta), its pressure is

    H=A+delta/A-sqrt(delta*(A^2+delta))/A
     =A B/(B+sqrt(delta)).

Positive rational squaring using .809<gamma<.81 gives

    A>.4966, B>.5146, sqrt(delta)<.1352,
    H>.4966*.5146/(.5146+.1352)>.393.

Therefore Delta>.393-.430*.81=.0447. Direct first-block variance gives
Q=gamma^2/10<.06561. Choose mu=.6. Its payment is<.0972<.1315 and
its gain in(5) exceeds

    1.2*.0447-.36*.06561=.0300204.

Adding E>.97977675 already exceeds1.

## 5. Finite native>=7

Here cos(pi/7)<=gamma<1 and gamma>.9. The actual vector

    xi=c+gamma^(-1)n1+2(gamma^(-1)-gamma)n2

has selected products xi.c=xi.n1=0 and
xi.n2=(3-4gamma^2)/(2gamma)<0. Its square is -delta with

    delta=7-4gamma^2-3/gamma^2>0,
    U=v-gamma^(-1)n1-2(gamma^(-1)-gamma)n2,
    H=A-sqrt(delta).

Direct block variance gives

    Q=(28/5)gamma^2-8+3/gamma^2<.6.                 (6)

This expression increases in gamma^2>.81 and tends to.6 at1.
We need the uniform bound Delta>.111, including the native7 endpoint.

For native7 the polynomial8z^3-4z^2-4z+1 is increasing on
[sqrt3/2,1] and vanishes at cos(pi/7). Its values at.9009 and.901
are -.000569738168 and .000257608. Thus .9009<gamma<.901.
Delta's individual terms have the conservative bounds

    A>.7388, sqrt(delta)<.2395, .430gamma<.38743,
    Delta>.11187>.111.                              (7)

The first is positive rational squaring in(4); delta increases on
this short interval, and its upper bound follows at.901. For example
.2395^2*.901^2-[7*.901^2-4*.901^4-3]
=.00004156271425>0.

For native>=8, gamma>.92. The function sqrt(delta) is strictly
concave where delta>0 because delta''=-8-18/gamma^4<0. At.92,
sqrt(delta)>.264 and delta'<.346, giving
(sqrt(delta))'<.346/.528<.66. Hence its derivative remains<.66
for gamma>=.92. Also A'>sqrt(40/21)>1.38, so

    (A-sqrt(delta)-.430gamma)'>1.38-.66-.43>0.

At.92, A>.7824 and sqrt(delta)<.2646 by rational squaring, and
Delta>.7824-.2646-.3956=.1222>.111. Together with(7) this covers
EVERY finite native>=7 without a label cutoff.

Choose mu=.2. Its payment is<.04<.1315. Equations(5)--(7) give gain

    2*.2*.111-.04*.6=.0204.

Consequently ||x||^2>.97977675+.0204=1.00017675, a contradiction.

## 6. Dotted gamma>1

Use the actual vector xi=c+gamma^(-1)n1. It has xi.c=0,
xi.n1=-gamma+gamma^(-1)<0, and square -(1-gamma^(-2)). Thus

    U=v-gamma^(-1)n1,
    H=A-sqrt(1-gamma^(-2)),
    Q=(8/5)gamma^2-2+gamma^(-2).

Put u=gamma^(-2) in(0,1), D=40/21, q=sqrt(D),
b=q-sqrt(D-1). Then Q/gamma^2=8/5-2u+u^2<=8/5. Concavity gives
sqrt(D-u)>=q-bu. The elementary circle/Cauchy bound gives

    bu+sqrt(u(1-u))<(.43+1.089)/2=.7595,

because q>1.38, q<1.381, sqrt(D-1)>.951 imply0<b<.43, and
sqrt(1+.43^2)<1.089. It follows that

    Delta/gamma
      =sqrt(D-u)-sqrt(u(1-u))-.430>.1905.

Choose mu=1/(8gamma). Its payment is exactly1/40<.1315. Its gain
in(5) exceeds

    (1/4)*.1905-(1/64)*(8/5)=.022625.

Again the total energy exceeds1. This treats the ENTIRE unbounded
dotted interval, not a sampled or truncated list.

All actual endpoint types are exhausted. Either the two ordinary grafts
or one actual isolated foot contradicts the original norm bound.
This proves the zero-weak-end theorem at the stated input boundary.
The one-weak-end case, whole B20/D20 and global19 remain open.
