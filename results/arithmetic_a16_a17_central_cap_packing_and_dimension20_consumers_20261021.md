# Arithmetic A16 and A17: central cap packing and dimension20 consumers

21 October 2026. Root's complete written proof candidate, held for
full outside-author review. No enumeration, optimizer, scientific,
symbolic, numerical, Lean or TeX job. This supports the existing compact
arithmetic paper, without a new manuscript or accepted ceiling change.

**Theorem.** At a farthest vertex of a compact arithmetic first-type
Coxeter polytope, an A16 component has actual projection energy>83/100,
and an A17 component has energy>6/7, with arbitrary finite spectators.
Neither component can occur at a farthest vertex in dimension20.

The central packing theorem below extends the accepted A18/A19 proof to
N17 and18. The shorter paths have EXTRA golden types and repeated norm
square classes; these are retained and explicitly resolved. The proof
does not extrapolate the earlier three-type norm table.

## 1. Inputs, field normalization and same-vector geometry

All paths below are under results/.

| Input | SHA256 |
| --- | --- |
| arithmetic_large_a_central_cap_packing_20261021.md: field-valued time-line argument and positive dual certificates, not its higher-rank packing conclusion | `8599baeda0f1cfae41bf3e4f27a89d7d999a49a3284caa8bdd559ea1ba945ee0` |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md: normalization, cyclic integrality and all minimal-field conjugates | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` |
| compact_a19_a20_central_window_actual_six_sevenths_20261020.md: complete finite-white inventory and actual isolated-cap proofs, not its larger-rank conclusions or adjacency estimates | `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` |
| compact_an_large_components_arbitrary_spectator_transfer_20261020.md: zero-own repair and positive clearing of entire spectators | `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` |
| compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md: original finite-white classification and actual H4/A2/singleton point formulas | `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb` |
| compact_dimension21_forest_frontier_and_small_debit_continuation_20261021.md: ordinary cost/defect interfaces only | `2edf213aa7bb25d7947a9d4712fdda28f95e9546a1af91e8f4cf3613c244e5bf` |
| arithmetic_compact_bound_research_reentry_frontier_20261020.md: complete necessary dimension20 dominant-rank list only | `5690e564ed601ff29c901e2dc1b6dfc6075ecabe819cd12ffbba6163f77232f0` |
| compact_large_an_rank14_five_sevenths_uniform_transfer_20261021.md: separate new A13/A14/A15 dimension20 consumers, used only in the final joint-frontier paragraph | `2af6d0032c431610a7ece98c3fc8d590e79b1bd75ade7d3c77486c2d51bb9f27` |

Standard external facts remain finite spherical classification, compact
acute simplicity, the native actual-wall gap and the arithmetic
first-type conjugate/cyclic-integrality criterion. No field census or
parked field-specific bound19 is used.

Put N=r+1 in{17,18}. Use unit A normals, duals and the SAME actual
farthest vector x, with ||x||<1. Its component heights satisfy

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=omega_i.omega_j=2i(N-j)/N (i<=j),
    h_i=x_A.omega_i, h0=hN=0,
    e_A=||x_A||^2=(1/2)sum_(i=0..N-1)(h_(i+1)-h_i)^2.

Call5<=i<=N-5 weak when h_i<=sqrt(D_i)-61/100. We first prove
that there is at most ONE such index, including every foreign-contact
and short-rank branch.

## 2. Complete weak inventory at N17 and18

The retained A paths both have rank at least4. Their one-node finite
extensions give exactly no contact, one admissible end/D/E simple3,
one end4, or two end3 contacts. The two-path connected rank is16 or17,
still>8. Hence the same finite tree/bond arguments apply; no exceptional
H/F/E completion appears at these ranks. All local affine principals
are present. The zero-own-contact repair uses only m>=5: the far-end3
height square is<1/(2m), the node2 foreign cost is Q<2/m, and the
end4 cost is Q<1/m; the full contacted spectator, if present, forces
the actual opposite to close a direction in an orthogonal component.
These formulas and the orthogonal-opposite contradiction do not
require N>=19. The gamma=0 affine arm alternatives likewise remain
present on the whole old A path. Thus gamma>0.

Here are fresh isolated-cap seals, needed because the old N19 lower
bound q>2.7 fails. At all present indices

    D_i>=120/17,
    q_i>2.6568, sqrt(D_i-1)>2.4614,
    a_i=1/[q_i+sqrt(D_i-1)]<500/2559<49/250.       (1)

These are positive square comparisons at120/17; their sum exceeds
5.1182>2559/500. For0<u<1 the actual cleared-point pressure gives

    h_i>sqrt(D_i-u)-sqrt(u(1-u)).

Its loss is<a_i u+sqrt(u(1-u))<=
(a_i+sqrt(1+a_i^2))/2<76/125<.61, using a_i<.196 and
sqrt(1+.196^2)<1.02. This covers isolated dotted caps too.

For foreign cost Q>0 and u=(1-Q)/gamma^2>1, positive clearing
gives h_i>sqrt(D_i-1)-sqrt(u-1). Since Q>=1/4 and gamma>=phi/2,
u-1<=5-3phi<73/500<.383^2. Equation(1) makes the loss<.579<.61.
The case u=1 is an actual feasible future null ray, excluded by
compactness. No cleared vector is assumed to be an actual facet.

For Q0 and finite native m>=8, the actual singleton finite projection
gives h_i>sqrt(D_i-1)-tan(pi/m). Using
sqrt2<70711/50000 makes its loss less than

    500/2559+20711/50000<61/100;                   (2)

the last integer cross margin is50051. The radical upper bound has
square numerator70711^2=5000045521>5000000000. For isolated5 the
actual H4 proof is unchanged: q_i+sqrt(D_i-25/16)>19/4 and
delta5<1/36 give loss<25/76+5/24<.61. For native7, the actual A2
core,9/10<gamma<10/11, gives loss

    <(100/81)/(q_i+sqrt(D_i-100/81))+9/25
    <5000/20493+9/25<1/4+9/25=61/100,

since q_i>2.65 and sqrt(D_i-100/81)>2.41. Isolated3/4/6 and
parallel contacts retain their explicit affine/null exclusions.
Therefore EVERY isolated cap is strong at both new ranks.

For all nonisolated finite-white roles,5/2<q_i<=3 and the edge target

    R(q)=100/61-(q-.61)/[q(2q-.61)]<3/2.           (3)

Indeed R increases for q>=5/2, and3/2-R(3)>0 has cleared numerator
1669. The complete ordinary-role seals from the inputs use only
m>=5 and R<761/500, so apply here. Adding spectator cost only
strengthens their edge lower bounds. For the two exceptional roles
the full foreign-clearing estimates are respectively
227/(25sqrt35)>761/500 and squared ratio>121/35. Thus a weak cap
has NO spectator contact. The near4 role at m>=8 is strong by the
old761/500 bound; at m7 its ratio>3/2, so (3) makes it strong too.

The COMPLETE remaining list is therefore

    near4: v=omega_i/sqrt2+omega_(i-1)/2, i5 or6,
           or its reflection;
    far helper: v=gamma*omega_i+omega_e/2,
           e=1 orN-1, Q0, own5/native>=7/dotted.

For a far helper of side parameter m, the exact edge ratio is
(gamma+1/(2m))/sqrt((m+1)/(2m)). Its square minus(4gamma-2) equals
(gamma-1-1/(2m))^2/[(m+1)/(2m)]. Equations(3) and this nonnegative
identity give

    202/125<t=2gamma<17/8.                        (4)

The lower bound is the native5 gap; larger surviving own labels and
dotted contacts only increase it. This is a fresh physical strip,
not the inherited larger-rank54/25 strip.

## 3. A norm/trace seal forces every far helper to be golden

Normalize the old A roots and a simple3-helper cap to square2 over
the ACTUAL ground field K. The helper product is-1 and the own
product is-t. Cyclic integrality makes t an algebraic integer. At
every nonidentity embedding of Q(t), extended to the totally real K,
strict positive Schur complement gives

    D_i t_sigma^2+2C_ie t_sigma+D_e<4.              (5)

Reverse to5<=i<=N/2. For negative t the cheapest helper gives
Q_-(t)=2i t^2+4t+2-2(it+1)^2/N; for positive t use
Q_+(t)=2i t^2+2-2(it-1)^2/N. Each increases with i on this half
and with N for the relevant sign. At the smallest N17,i5,

    Q_-(-98/125)-4=396/53125>0,
    Q_+(12/25)-4=156/2125>0.

Their values at0 are<4. Thus EVERY other conjugate satisfies

    -98/125<t_sigma<12/25.                        (6)

Let d=[Q(t):Q]. At d>=5 the nonzero integer norm has absolute value
less than(17/8)(98/125)^4<1. At d4 its absolute value is<2. If any
small conjugate were positive, the product would have magnitude at
most(17/8)(12/25)(98/125)^2<1. All three must therefore be negative,
say-a,-b,-c, with abc=1/t. The integral trace satisfies

    -1<202/125-3(98/125)<t-a-b-c
       <=t-3t^(-1/3)<0,

where the last inequality follows from t^4<(17/8)^4<27. Impossible.
At d3 the norm likewise has absolute value<2; a positive small
conjugate would give magnitude<(17/8)(12/25)(98/125)<1. Thus both
are negative with ab=1/t, and the integral trace lies strictly in
(0,1): its lower bound is6/125, and its upper bound is
t-2/sqrt(t)<1 because t(t-1)^2<(17/8)(9/8)^2<4.

At d2 the integral trace is1 or2. Trace2 gives nonzero norm
t(2-t) strictly between-1 and1. Trace1 forces t<223/125 and norm
t(1-t) in(-2,0), hence norm-1 and t=phi. At d1 only t2 is possible;
the own unit product-1 is an actual feasible future null ray and
contradicts compactness. We conclude t=phi in every case.

The Q0 old component plus cap has Lorentzian Gram over Q(phi).
A proper totally real extension K/Q(phi) would have a nonphysical
embedding fixing that same Lorentzian principal Gram. Consequently
the ACTUAL field is K=Q(phi). No field list is assumed.

## 4. Near4 is impossible, and all golden orientations are retained

For near4 the same simple3 normalization puts sqrt2 in K. Its
conjugate time numerator M_i'=2N A_i'^2 is

    M5'=(12-8sqrt2)N-66+40sqrt2,
    M6'=(15-10sqrt2)N-97+60sqrt2.

Both increase with N. At17 they are138-96sqrt2 and158-110sqrt2,
strictly positive using sqrt2<17/12. But the conjugate new-root
Schur complement is-M_i'/N and must be positive. This excludes
every near4 weak orientation.

Now reverse a golden cap to5<=i<=N/2. For a helper on the nearer
end set

    M_i=2N A_i^2=i(N-i)-(N+1)+(i+2)(N-i)phi;

for a helper on the opposite end set

    M_i^far=i(N-i)-(N+1)+i(N-i+2)phi.              (7)

At N18 the midpoint i9 has conjugate unit Schur cost
(197-99phi)/36>1. The opposite-end costs on6<=i<=9 are a concave
quadratic in i whose endpoint values are(173-84phi)/36 and the
preceding midpoint value; both>1 using phi<13/8. Thus the complete
necessary list is nearer-end i5,6,7,8 and opposite-end i5, with
their reflections. The new nearer8 and opposite5 are NOT discarded.

At N17 representatives are i5,6,7,8. Only opposite-end i7 is
excluded: its conjugate Schur cost is(170-84phi)/34>1 by
phi<34/21, whose defining polynomial has positive numerator1.
Keep all nearer i5..8 and opposite i5,6,8. Some may be geometrically
strong as well, but none is silently removed from this proof.

## 5. Different golden types cannot coexist

The cap and its orthogonal projection onto the full old vertex forest
are K-valued. Since Q0, their residual lies on the SAME one-dimensional
time line over K. Any two caps therefore have A_i/A_j in K, so the
ratio of their M numerators must be a square in K. Its rational norm
must in particular be a rational square. This remains true with
arbitrary spectators and does not presuppose a whole A vertex.

At N18 the complete table from(7) is

| Type | M | Norm M |
| --- | --- | ---: |
| nearer5 | 46+91phi | -1979 |
| nearer6 | 53+96phi | -1319 |
| nearer7 | 58+99phi | -695 |
| nearer8 | 61+100phi | -179 |
| opposite5 | 46+75phi | -59 |

The five positive magnitudes are pairwise coprime and individually
nonsquares; their pairwise ratios are not rational squares. Thus no
different types coexist.

At N17 the complete table is

| Type | M | Norm M | Norm square class |
| --- | --- | ---: | ---: |
| nearer5 | 42+84phi=42phi^3 | -1764 | -1 |
| nearer6 | 48+88phi | -1216 | -19 |
| nearer7 | 52+90phi | -716 | -179 |
| nearer8 | 54+90phi=18phi^5 | -324 | -1 |
| opposite5 | 42+70phi=14phi^5 | -196 | -1 |
| opposite6 | 48+78phi=6phi^7 | -36 | -1 |
| opposite8 | 54+88phi | -76 | -19 |

Different norm classes cannot mix. Within the-1 class the odd phi
powers differ by squares, leaving rational coefficients42,18,14,6.
Their square classes are42,2,14,6; each pair ratio has an odd3 or7
valuation. A positive rational is a square in Q(sqrt5) only if it
is a rational square or five times a rational square: expand
(a+b sqrt5)^2 and use2ab=0. Neither possibility occurs here.

Within the-19 class the ratio opposite8/nearer6 is

    (25+33phi)/76=(83+33sqrt5)/152,

with norm1/16. If it were(a+b sqrt5)^2, then
a^2=(83+38)/304=121/304 or(83-38)/304=45/304. Both have an odd
denominator valuation at19 and are not rational squares. Hence these
two types cannot mix either. Repeated norm square classes have thus
been checked in K, not dismissed by the insufficient norm test.

## 6. Reflected copies of the same type cannot coexist

At a fixed index the actual opposite is unique by compact simplicity.
For reflected caps of the same type, equal positive time coefficients
and the inverse-A formula give k=-2c_L.c_R as

    nearer: k=[(N-2i)(i+(i+2)phi)-(N+2)]/N;
    opposite: k=[(N-2i)(i+(i-2)phi)-(N+2)]/N.     (8)

The simple3 helpers and the old3 path joining them make k integral by
the cyclic criterion, even if the mutual wall product is dotted.
For nearer i5,6,7,8 the constant numerators in(8) are
16,11,2,-11 at N17 and20,16,8,-4 at N18. None is divisible by
the respective N. For opposite i5,6,8 at17 they are16,11,-11;
for opposite5 at18 the constant numerator is20. Again none is
divisible by N. Since O_K=Z[phi], all these k are nonintegral,
a contradiction. No native angle is assigned to a cleared foot.

Sections2--6 prove the COMPLETE packing theorem: among all indices
5..N-5 there is at most ONE weak direction, for both N17 and18.
No adjacency lemma or unproved smaller-rank C_(i,i+1) bound is used.

## 7. Exact same-vector energy certificates

A concave rational height reference H0=HN=0 defines
w_j=(H_j-H_(j-1))/sqrt2. Its nonnegative dual coefficients are the
half curvatures. They occur only at available strong constraints,
so x_A.w>||w||^2 and e_A>||w||^2. Internal omission at j replaces
H_j by the average of its neighbors and loses exactly
(H_(j-1)-2H_j+H_(j+1))^2/4. With at most one hole there is no
omission-interaction assumption.

For N17 take

    H5=2.046,H6=2.176,H7=2.259,H8=H9=2.300,
    H_(17-i)=H_i,

with linear tails from0 to5 and12 to17. The increments.130,.083,.041
are positive and decreasing after the initial slope.4092. All floors
are strictly below sqrt(D_i)-.61 by positive square comparisons at
120/17,132/17,140/17,144/17. Their energy is.8626932. An internal
missing index loses at most.047^2/4=.00055225. If5 is missing,
reverse if needed and use the linear tail0..6 instead; its remaining
first slope2.176/6 is still greater than the next increment.083.
The resulting energy is EXACTLY

    .8626932-(1/2)[2.046^2/5+.130^2-2.176^2/6]
       =6226597/7500000>83/100.                    (9)

For N18 use

    H5=2.077,H6=2.218,H7=2.314,H8=2.371,H9=2.389,
    H_(18-i)=H_i,

again with linear end tails. The successive positive increments are
.141,.096,.057,.018; concavity and all strict floor squares follow
at65/9,8,77/9,80/9,9. Its energy is.8954558. An internal omission
loses at most.045^2/4=.00050625. A missing5 instead gives EXACTLY

    .8954558-(1/2)[2.077^2/5+.141^2-2.218^2/6]
       =12961241/15000000>6/7.                    (10)

Both references remain concave after the stated endpoint omission;
every positive curvature belongs to a strong index. A single weak
index anywhere else has the smaller internal loss just bounded.
This proves both component energy conclusions. All displayed decimals
are exact terminating rationals, not numerical output.

## 8. Complete dimension20 consumers and joint frontier

An A16 component has four spectator ranks. Ordinary/J alternatives
cost>2/7. One short I leaves ordinary total rank three, cost>17/70.
Either typed7 block costs>.09 and leaves ordinary rank two, cost>1/7.
Each exceeds the remaining17/100 budget. Thus A16 is excluded with
EVERY permitted spectator forest.

An A17 component has three spectator ranks. Ordinary/J cost>17/70.
One I leaves ordinary rank two, cost>1/7; with(10) this strictly
exceeds1 without charging I. A typed7 block leaves one nonshort A1,
cost>.1, and contributes>.09;6/7+.19>1. All are excluded. The
ordinary total-rank-three cost is17/70, not the false universal1/4;
no intrinsic dihedral7/32 rate is used.

For the joint dimension20 frontier, the accepted reentry proof first
leaves A13..A20 or the classical rows. The separate rank14 five-sevenths
candidate excludes A13/A14/A15 by its explicit consumers; this note
excludes A16/A17; accepted packing8599baed excludes A18/A19; and the
accepted whole-A20 summed-foot proof457debf9 excludes the whole A20
case. Accepted classical continuation/B20/D19+I results already leave
only whole D20 in the B/D sector. Therefore, AFTER complete reviews
and acceptance of both new candidates, the entire arithmetic dimension20
frontier reduces to ONLY WHOLE D20. This is not yet a ceiling19 theorem:
the D20 case must be excluded separately and the global composition
then reviewed. The current accepted ceilings stay57/21/20.

No primary manuscript, PDF, frozen packet or formal source is changed.
