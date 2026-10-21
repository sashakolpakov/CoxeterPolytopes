# Arithmetic large A components: central cap packing and rank-18/19 consumers

21 October 2026. Written proof candidate by root; full outside-author
review is required before acceptance. This is one supporting mathematical
note for the existing compact arithmetic paper, not another manuscript.
No enumeration, optimizer, root search, TeX or Lean job was used.

## Statements and input boundary

Let a vertex of a compact arithmetic first-type Coxeter polytope be
farthest from a strict interior point. Suppose its finite forest contains
A_(N-1), with N in {19,20,21}, and any finite spectator forest. Use UNIT
old normals, dual vectors omega_i, and the common farthest vector x:

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=omega_i.omega_j=2i(N-j)/N (i<=j),
    h_i=x_A.omega_i, h_0=h_N=0,
    e_A=||x_A||^2=(1/2)sum_(i=0..N-1)(h_(i+1)-h_i)^2<1.

An index in 5<=i<=N-5 is weak if h_i<=sqrt(D_i)-61/100.

**Central packing theorem.** There is at most ONE weak index in this
window. Every weak cap is a far helper: its own angle is pi/5, its
helper is on the nearer end of the path, and i is 5,6,7 or a reflection
of one of these. Its actual ground field is Q(sqrt5). The exceptional
near4/simple3 weak profile is impossible by conjugate positivity.

**Dimension-20 consumers.** In a dimension-20 compact arithmetic
first-type chamber, the farthest forest cannot contain A18 or A19.
These conclusions do not exclude whole A20 or the smaller A components,
and do not establish an arithmetic dimension ceiling of 19.

The arithmetic inputs are the usual first-type field-valued normals,
cyclic integrality and positive definiteness at every nonphysical real
embedding. The geometric inputs below are already proved internal
lemmas, not extra classification assumptions. Finite spherical
classification, compact acute simplicity and the native wall gap remain
the stated external inputs.

| Internal input, under results/ | SHA256 |
| --- | --- |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` |
| compact_a19_a20_central_window_actual_six_sevenths_20261020.md | `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` |
| compact_an_large_components_arbitrary_spectator_transfer_20261020.md | `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` |
| compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md | `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb` |
| compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md | `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b` |
| compact_classical_rank21_isolated_i_actual_terminal_obstruction_candidate_20261021.md | `8183924e37e5153a13a141daa6b7b2af44e3998b4008695ba3b1a0eef1d29bf2` |
| compact_a20_small_debit_coupling_exclusion_20261021.md | `6ed1ed3c81f5126c6ebf66c07cf8fc3338db3b64352d8985af09db74c6ee9721` |

Only the universal I/7 height statements of the last two inputs are
used, not their higher-rank component conclusions. The accepted forest
accounting supplies the exhaustive spectator rank-two list: ordinary
rank two; one I plus a nonshort A1; a type-J pair; or either rank-two 7
defect. The respective ordinary/J/7 costs used below are >1/7, >.34,
and >.09; a nonshort A1 costs >.1. There is no false intrinsic 7/32
bound for every ordinary dihedral component.

## 1. The weak inventory at N19,20,21 without a smaller energy premise

We use only Sections 1--2 of the pinned A19/A20 source, before its
adjacent-packet and near4-removal arguments. Their inventory proofs use
e_A<1, not the later contrary premise e_A<=6/7. Here are the rank and
scalar checks needed also for N19.

The two retained paths have rank at least4, and any connected two-white
completion has rank at least18>8. Thus the same finite classification,
zero-own-contact repair and affine principal obstructions apply. The
retained helper-side length is always m>=5. For isolated caps,

    D_i>=140/19, q_i=sqrt(D_i)>27/10,
    sqrt(D_i-1)>5/2, q_i+sqrt(D_i-1)>26/5.

The foreign-cost/u>1 loss is <95/156<.61. The u<1 clearing loss is
at most (5+sqrt701)/52<.61. Isolated native m>=8 has the same
95/156 bound. For native7, use

    sqrt(D_i-100/81)>247/100,
    loss <10000/41877+9/25<61/100.

The squared first inequality follows from 140/19-100/81>(247/100)^2.
The isolated5 seal still uses q_i+sqrt(D_i-25/16)>19/4.
Consequently every isolated central cap is strong.

For the nonisolated finite-white roles, q_i lies in [5/2,31/8]. The
edge target F(q)/q is <761/500 throughout this interval, where

    F(q)=[1+(1-.61/q)^2]/[1-(1-.61/q)^2].

All ordinary role seals in the pinned proof use just m>=5 and this
interval. Clearing a contacted spectator costs at least1/4; its far
helper and near4 estimates in Section4 of the spectator-transfer input
are still >761/500. Hence a weak cap has NO spectator contact.

The complete remaining weak profiles are therefore exactly:

    far helper: v=(t/2)omega_i+omega_e/2,
                e=1 or N-1, 202/125<t<54/25;
    near4:      v=omega_i/sqrt2+omega_(i-1)/2,
                i=5 or6, or the reflected profile.

The near4 formula means OWN4 and WHITE SIMPLE3, in that order. An
own3/white4 profile is strong. To check the stated near4 locations,
the m>=8 role is strong, while at m7 its edge ratio exceeds3/2.
The target F(q)/q increases for q>=5/2: its derivative has the sign
of 2q^2-4(.61)q+(.61)^2. Its value at25/8 is <3/2. This also
handles N19's q7<3, where the earlier source's displayed [3,25/8]
interval was unnecessarily restrictive. The near4 actual mate is at
i-1, hence4 or5, as required by the pinned mate lemma.

No adjacent weak-cap packet estimate is needed in this note; in
particular we do not import a bound C_(i,i+1)>=7 at N19.

## 2. Every central weak far helper is golden

The normalization of the pinned uniform arithmetic input gives old
A roots and c square2 over the ACTUAL ground field K, helper product
-1 and marked product -t. Cyclic integrality makes t an algebraic
integer. Every nonphysical conjugate satisfies the strict Schur bound

    D_i*t_sigma^2+2*C_ie*t_sigma+D_e<4.                (1)

By reversal take 5<=i<=N/2. For negative t the least cost over the
two helpers is

    Q_-(t)=2i*t^2+4t+2-2(it+1)^2/N;

for positive t it is Q_+(t)=2i*t^2+2-2(it-1)^2/N. Each increases
with i on this half-window and with N for the relevant sign. At
N19,i5, Q_-(-193/250)-4=643/59375>0 and
Q_+(12/25)-4=232/2375>0. Since each quadratic's value at0 is
below4, every nonidentity conjugate in Q(t) lies in

    -193/250<t_sigma<12/25.                           (2)

Every embedding of Q(t) lifts to the totally real K, so none is missed.
If d=[Q(t):Q]>=4, the nonzero integer norm has absolute value less
than (54/25)(193/250)^3=388209078/390625000<1. If d3, its absolute
norm is <2. Opposite signs of the two small conjugates give absolute
norm <(54/25)(193/250)(12/25)<1; two positive conjugates do too.
Thus both are negative, say -a,-b. The norm is1, so ab=1/t. Their
integral trace satisfies

    0<202/125-2(193/250)<t-a-b<=t-2/sqrt(t)<1.

The last inequality follows from t(t-1)^2<4 at t<=54/25, and
monotonicity; this is impossible.

In degree2 the integral trace is1 or2. Trace2 gives a nonzero norm
t(2-t) strictly between-1 and1. Trace1 and (2) give t<443/250 and
-2<t(1-t)<0, so the norm must be-1. Thus t=phi=(1+sqrt5)/2.
In degree1 the only physical integer is2, the forbidden compact
parallel pair. We have proved t=phi.

The helper must be on the nearer end. Indeed, still taking i<=N/2,
the OPPOSITE helper N-1 has unit Schur cost at phi'=1-phi

    S=[i(N-i)*phi'^2+2i*phi'+N-1]/(2N).

For fixed N this is concave in i, so its minimum on [5,N/2] occurs
at an endpoint. At i5 it increases with N, since its subtracted term
is (5phi'-1)^2/(2N). At i=N/2 it is
N*phi'^2/8+phi'/2+1/2-1/(2N), also increasing. At N19 the two
values are (84-40phi)/19 and (870-437phi)/152, both >1 using
phi<13/8. This contradicts nonphysical positive definiteness.
The pinned uniform lemma excludes every i>=8 and i<=N-8.
Thus only i5,6,7 with helper1, or their reflections, remain.

## 3. A pure weak cap locks the actual ground field

For either profile the old A component plus c has physically
Lorentzian Gram: c=A*p-v, Q=0 and A>0, so its Schur complement
is negative. In the square-two normalization that Gram has entries
in F=Q(phi) for the far helper, and F=Q(sqrt2) for near4.
The native simple3 helper normalizes c over K, without adjoining
a square root. Thus F is a subfield of K.

If K were a proper extension of F, an additional real embedding of K fixing
F would repeat this Lorentzian Gram at a nonphysical embedding,
contradicting positive definiteness. Such an embedding exists because
K/F is totally real and has degree at least2. Therefore K=F.
In particular the golden and silver weak profiles cannot coexist.

## 4. Two golden weak caps cannot coexist

All normals are now square2 over K=Q(phi). For a left cap at i,
put M_i=2N*A_i^2. Direct substitution of the unit inverse entries gives

    M_i=i(N-i)-(N+1)+(i+2)(N-i)*phi.

The three norms, in order i5,6,7, are

| N | Norm M5 | Norm M6 | Norm M7 |
| --- | ---: | ---: | ---: |
| 19 | -2204 | -1420 | -656 |
| 20 | -2439 | -1519 | -599 |
| 21 | -2684 | -1616 | -524 |

For fixed N the ratio of any two entries is not a rational square.
For N19 use the respective square classes -19*29,-5*71,-41;
for N20 use -271,-31,-599 (271 and31 are prime and divide neither
other class); for N21 use -11*61,-101,-131.

Why is this a necessary arithmetic test? Subtract from each K-valued
normal its K-valued projection onto the old A roots. With Q0 its
residual is along the same one-dimensional time line. Any two such
nonzero K-vectors have a ratio in K, so A_i/A_j is in K. Therefore
M_i/M_j must be a square in K, and its rational norm a rational square.
This argument is unchanged by spectators: the caps are orthogonal to
them and the projections use the rational A Cartan inverse. Consequently
two weak caps must have the SAME i up to reflection.

At a fixed i the only two different possibilities are reflected caps.
Their equal positive time coefficients give, with k=-2*c_L.c_R,

    k=[(N-2i)*(i+(i+2)phi)-(N+2)]/N.                 (3)

The simple3 helper bonds and the old3 path between them make k an
algebraic integer by the cyclic criterion. Since O_K=Z[phi], (3)
would require N to divide BOTH 2(i^2+1) and 2i(i+2). For i5,6,7
these pairs are (52,70),(74,96),(100,126), with gcd2 in every
case. This is impossible for N>=19. Actual opposites at a single
index are unique, so there is at most one golden weak index.

## 5. The weak near4 orientation is impossible

The literal profile is v=omega_i/sqrt2+omega_(i-1)/2, i5 or6,
or its reflection; OWN4 is essential in this formula. Set s=sqrt2.
Now

    M_i=2N*A_i^2
       =2i(N-i)+(i-1)(N-i+1)-2N+2(i-1)(N-i)*s.

At the embedding s maps to -s, conjugate positive definiteness requires
the Schur complement to be positive. But this complement is
-2*A_i'^2=-M_i'/N. In fact M_i'>0 for both i and every N>=19.
For i5,

    M5'=(12-8s)*N-66+40s,

which increases with N. At N19 it is162-112s>0, since162^2-
2*112^2=1156>0. For i6,

    M6'=(15-10s)*N-97+60s,

also increasing with N, with N19 value188-130s>0 because188^2-
2*130^2=1544>0. Both positive leading coefficients follow from
s<3/2. Thus every exceptional weak near4 profile contradicts its
nonphysical positive definite Gram. Reflections give the same result.

The orientation matters: the actual weak profile has the white simple3
contact at i-1. The earlier mate lemma also discusses the opposite
orientation with contact i+1, which can be conjugately admissible but
is already strong in this central range. We have not declared every
two-contact4/3 wall impossible. Sections 3--5 prove the stated theorem.

## 6. A slightly stronger universal I barrier

The isolated-short-I lemma supplies numbers 0<u<1/3, t>0 such that
the I energy is >=u^2, every other actual edge has normalized
projection at least t, and

    (1-u^2)(1-t^2)=8ut.

We will use

    t+(3/2)u>767/1000.                              (4)

To prove it, put t0=.767-1.5u>0. The left-minus-right at t0 is

    P(u)=.411711-3.835u+9.338289u^2-2.301u^3+2.25u^4.

For v=u-.2 in(-.2,2/15), the coefficients are
.00343456,-.3038044,8.497689,-.501,2.25. The cubic is
>=-.0668v^2, so

    P(u)>.0034-.304|v|+8.43v^2>0.

The last quadratic has negative discriminant:
.304^2-4(.0034)(8.43)<0. The original expression decreases
strictly with positive t, proving (4). All decimals in this note
are exact terminating rationals, not floating-point evidence.

## 7. Rational energy certificates and how endpoint gains are counted

Use zero-sum coordinates x_1,...,x_N with
h_i=sqrt2*sum_(k<=i)x_k. A rational concave lower height sequence
H_i defines w_k=(H_k-H_(k-1))/sqrt2. At each strict bend the
corresponding height constraint is available, so summation by parts
gives x.w>=||w||^2 and ||x||^2>=||w||^2. It is enough to use
strictly smaller rational floors at every strong index.

Removing an internal height constraint at j replaces H_j by the
average of its neighbours and lowers energy by
(H_(j-1)-2H_j+H_(j+1))^2/4. All other lower heights are preserved.

On a first block of length m with constant coordinate b, an actual
endpoint bound x1>=z>b raises this dual energy by

    C*(z-b)^2, C=m/(m-1),                           (5)

provided the modified remaining m-1 equal coordinates meet the next
coordinate. This is the direct reference vector with first coordinate
z and the next m-1 coordinates (mb-z)/(m-1). All later prefix sums
are unchanged, the new first-coordinate multiplier is nonnegative,
and its norm increase is exactly (5). Opposite-end modifications on
disjoint blocks add. These are changes in a single common energy
certificate, not independently charged overlapping gains.

For I use z=a*t, a=sqrt((N-1)/N). Equations (4)--(5) and weighted
Cauchy give, with B=.767,

    u^2+C*max(a*t-b,0)^2
       >(a*B-b)^2/[(9/4)*a^2+1/C].                 (6)

The same formula holds with the weaker inherited B=.75. It follows
from max(a*t-b,0)+(3/2)*a*u>a*B-b. Each reference below admits
the endpoint modification for every z<=a<1: even at z1 the
remaining block coordinate is larger than the next coordinate.

## 8. Excluding A19 in dimension20

Here N20 and the spectator has rank1. If there is no weak5 at either
end, the only possible missing index is6 or7, or a reflection. Start
with the symmetric rational heights

    H5=2.128,H6=2.288,H7=2.406,H8=2.488,
    H9=2.536,H10=2.552,H_(20-i)=H_i,

and linear tails from0 to5 and15 to20. They are concave and strictly
below every respective strong floor sqrt(D_i)-.61, by positive square
comparisons. Their energy is .9544848. A missing6 loses .000441;
a missing7 loses .000324. Therefore e_A>.95404. An ordinary
spectator costs >.1. For I, the first5 mean b<.302, a>.974, and
(6) with B=.75 has numerator >(.75*.974-.302)^2 and denominator
<3. Its gain is >.06. Thus I also contradicts e_A+e_I<1.

If5 is weak, reverse if necessary to put it on the left. All indices
6..15 are strong. Keep the preceding rational heights from6 to15
and use a linear left tail from0 to6. Its energy is

    (1/2)[2.288^2/6+2(.118^2+.082^2+.048^2+.016^2)
                         +.160^2+2.128^2/5]
    =.9250917333...>.92509.

The ordinary spectator again closes the case. For I the first6
mean b=2.288/(6sqrt2)<.270 and a>.974. Formula (6) has numerator
>(.767*.974-.270)^2=.477058^2 and exact denominator713/240.
Their ratio is >.0765. Hence the SAME actual A-plus-I energy is
>.92509+.0765>1.00159. This excludes A19 with every rank1
spectator. No independent intrinsic positive cost for I was assumed.

## 9. Excluding A18 in dimension20

Here N19 and spectator rank2. The symmetric rational heights are

    H5=2.104,H6=2.255,H7=2.363,H8=2.433,
    H9=H10=2.467,H_(19-i)=H_i,

with linear end tails. They are concave strict lower floors; their
energy is .9258842. If5 is strong at both ends, a possible missing6
loses .00046225 and a missing7 loses .000361. Thus e_A>.9254.
Every allowed rank2 spectator alternative has energy >.09 (even
discarding I in I+ordinary-A1). This already contradicts the norm.

If5 is weak, orient it left. Keep H6..H14 and use a linear first6
tail. The energy is

    (1/2)[2.255^2/6+2(.108^2+.070^2+.034^2)
                          +.151^2+2.104^2/5]
    =.8955541833...>.8955.

Ordinary rank2 costs >1/7 and J costs >.34. For I plus an ordinary
A1, apply (6) with B=.75, a>.97, b<.27 and denominator<3 to
obtain joint A-endpoint/I gain >.06; the A1 additionally costs>.1.
All these cases close.

For either rank-two 7 defect its energy is >.09. The typed height
propagation lemma gives every other normalized edge-height square
theta>(11-4sqrt7)/3>11/80. The A2/own7 alternative has the stronger
theta>.2. The exact bound11/80 follows from sqrt7<847/320; the
positive square difference is609/102400. At each A endpoint,

    x1^2, x19^2 > (18/19)(11/80)=99/760>81/625,

so x1>.36 and -x19>.36. The first6 mean is <.266 and the last5
absolute mean is <.298. Applying (5) to these DISJOINT blocks
gives additional energy greater than

    (6/5)(.36-.266)^2+(5/4)(.36-.298)^2
       =.0106032+.004805=.0154082.

Consequently total energy >.8955+.09+.0154082>1.0009082.
This exhausts the spectator rank2 list and proves the A18 consumer.

## Verification scope and what remains open

The argument is a finite sequence of explicit algebraic identities,
integer-norm obstructions and rational dual-energy inequalities. It
requires no combinatorial enumeration or script-output assumption.
The small tables are displayed in full so that their entries can be
checked by substitution. Full independent source review is pending;
the note is not yet an accepted update to any canonical bound.

Even after acceptance the universal arithmetic compact ceiling remains
20. The remaining smaller A ranks, whole A20 and the open classical
terminal alternatives require separate complete arguments. The six
finished papers and their selected Lean snapshots are unchanged.
