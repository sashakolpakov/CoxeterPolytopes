# Whole A21 farthest-vertex exclusion: proof candidate

20 October 2026. PAPER ONLY. Awaiting independent mathematical review.
No scientific execution, numerical optimization, Lean, TeX, new uniform
dimension bound, or extension to arbitrary spectator components is claimed.

**Candidate theorem.** A vertex maximizing distance from a fixed strict
interior point of a compact hyperbolic Coxeter21-polytope cannot have
WHOLE finite diagram A21.

This note supplies complete central endpoint coverage and the subsequent
energy assembly. The following already written geometric inputs are used:

* Actual positive acute projection and farthest comparison:
  `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
  SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`.
* The actual near4 mate and the geometric far5 distance-one prohibition:
  `compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`,
  SHA256 `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`.
* Complete adjacent weak far-m/simple3 packet, INCLUDING both far-end
  assignments and all native cap products:
  `compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md`,
  SHA256 `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7`.
  Its outside source review is
  `review_compact_a21_adjacent_weak_far_packet_batch_guard_20261020.md`,
  SHA256 `a64d90c7b01fda3de511c9230f8f85646e51701400e197d378e7daf341f58b42`.

All paths in this note are under results/. No arithmeticity is assumed.

## 1. Common vector, exact coordinates, and actual-point pressures

Use UNIT outward old normals n_1,...,n_21, with adjacent product-1/2.
Let omega_i be the spatial duals, h_i=x.omega_i, h_0=h_22=0, where x
is the spatial farthest vector. Then

    ||x||<1, x.n_i>0,
    D_i=omega_i^2=2i(22-i)/22,
    C_ij=omega_i.omega_j=2i(22-j)/22, i<j,
    x.n_i=h_i-(h_(i-1)+h_(i+1))/2,
    ||x||^2=(1/2)sum_(i=0)^21(h_(i+1)-h_i)^2.      (1)

Thus the h_i are strictly concave. Set q_i=sqrt(D_i), c0=61/100,
g_i=q_i-c0. An index is weak if h_i<=g_i; all other indices are
called strong. A strong central index has h_i/q_i>39/50, because
q_i>=sqrt(85/11)>61/22 for5<=i<=17.

For an actual outside unit wall write c=A p-v, A>0, v^2=A^2+1.
For any finite old set J, let beta_j=-c.n_j, w=B_J^-1 beta in its
old-normal coordinates, and delta=1-beta^T B_J^-1 beta. If delta>0,
actual finite-face projection gives

    x.(v-w)>=H(A,delta),
    H(A,delta)=A/[1+sqrt(delta/(A^2+delta))].        (2)

If delta<0, the vector xi=c+w is future timelike and ACTUAL feasible:
xi.J=0, xi.c=delta<0, and every other wall has nonpositive pairing
with this positive combination of actual normals. Its norm is delta.
Farthest comparison gives instead

    x.(v-w)>=A-sqrt(-delta).                      (3)

Minimality of a Lanner core is unnecessary for (3). We also use the
actual dotted-pair foot xi=c+n_i/gamma when c.n_i=-gamma<-1:

    x.(v-n_i/gamma)>=A-sqrt(1-gamma^-2).          (4)

Its norm is -(1-gamma^-2), its pairing with c is0, with n_i is
-gamma+gamma^-1<0, and with all other walls is nonpositive. These
are actual feasible points, not merely inverse-Gram directions.

## 2. Complete central endpoint inventory

Fix5<=i<=17. At its actual opposite c, the retained white forest is
A_(i-1)+A_(21-i), each component having rank at least4. Put
gamma=-c.n_i. Any finite component containing c is a tree, so c
contacts at most one node of each white path. Classification of the
finite Coxeter diagrams gives the COMPLETE possibilities:

1. No white contact (isolated).
2. One white simple3: at an end of its path; at the second node from
   an end (D); or at the E6/E7/E8 admissible internal nodes.
3. One white4, necessarily at a path end.
4. Two white simple3 contacts, each at a path end, joined through c.

Here a path plus c has rank at least5, excluding H and F exceptions.
In case4 its connected rank is21, so only A/B/D could occur. A branch
in an old path has one old arm of length at least2 and the other arm
through c of length at least5; it cannot be D, whose two short arms
both have length1. An internal4 bond cannot be B. Thus case4 is A21.

In all cases gamma>0. If gamma=0, a single old3 contact at a full
path end or next-to-end gives finite A22 or D22, impossible in
signature(21,1). A more interior3 contains an affine star with arms
(1,2,5) or(1,3,3). A single4 gives finite B22 at the old path end,
or the affine star with squared incident cosines1/2,1/4,1/4 internally.
Two white3 contacts give an all3 affine cycle. This also records all
zero-own-contact branches rather than assuming an opposite must have
a nonzero contact with its released normal.

The individual edge has

    cosh L_i=sigma*sqrt(D_i)/sqrt(delta),          (5)

where delta is the white Schur residual and sigma the coefficient
remaining along omega_i. For a white path of length m-1 on one side,
number its nodes j=1,...,m-1 from the FAR end, so m=i or22-i.
For a white3 contact,

    delta=1-j(m-j)/(2m), sigma=gamma+j/(2m).      (6)

For an end white4, delta=1/m and sigma=gamma+j/(sqrt2*m).
For two white3 contacts, delta=1/D_i and sigma>=gamma+1/D_i.

The following elementary comparison is useful: for q in[5/2,10/3],

    [1+(1-c0/q)^2]/[q*(1-(1-c0/q)^2)]<38/25.      (7)

Clearing positive denominators gives the concave polynomial
-36400q^2+163602q-93025, positive at both endpoints. Consequently
sigma/sqrt(delta)>38/25 makes the index strong.

For non-near white3 contacts other than the far end, gamma>=phi/2>
f=101/125: own3 gives an affine cycle; own4 gives the affine star
above; own6 gives an affine6,3 path. The D positions have
delta=2/m and sigma>=gamma+1/m, whose ratio is minimized at m5
and exceeds38/25. The E positions have delta<=1/4 and sigma>=gamma,
again sufficient. A near-end3 with own label at least5 has
delta<=3/5 and sigma>=f+2/5, also sufficient.

A far-end white4 with own3 contains the affine path3,4,3,3; own4/6
are excluded by the same local affine stars/paths. Thus gamma>=f,
delta<=1/5 and sigma>=gamma, which is sufficient. For near-end white4,
gamma>=1/2, sigma>53/50 and delta<=1/5, also sufficient. In case4,
gamma>=1/2 and sigma/sqrt(delta)>=q_i/2+1/q_i>38/25.

### 2a. Isolated caps are all strong

Isolated3,4,6 and parallel contact are excluded by affine subdiagrams:
for3 use old arms(1,3,3), for4 the squared-cosine star above, and
for6 the6,3 path. The remaining finite labels are5 or m>=7.

For5, (2) with the actual H4 segment on either side gives
delta5=(7-3sqrt5)/16<1/36. Since gamma>4/5,

    q_i-h_i <25/76+5/24<61/100.

Indeed H(A,delta)>A-sqrt(delta), A/gamma=sqrt(D_i-gamma^-2),
and q_i+sqrt(D_i-25/16)>19/4.

For7, use (3) with the old A2 adjacent to i. It has
-delta=4gamma^2/3-1, and 9/10<gamma<10/11. Thus

    q_i-h_i <400/1701+9/25<61/100.                (8)

The first summand bounds q_i-sqrt(D_i-gamma^-2), using
q_i>11/4, sqrt(D_i-100/81)>5/2; the second bounds
sqrt(4/3-gamma^-2), since gamma^-2>121/100.

For m>=8, (2) with J={i} gives

    h_i>sqrt(D_i-1)-tan(pi/m),
    q_i-h_i<4/21+5/12=17/28<61/100.              (9)

Here q_i>11/4, sqrt(D_i-1)>5/2 and tan(pi/8)=sqrt2-1<5/12.
For a dotted isolated cap, (4) and x.n_i>0, with u=gamma^-2, give

    q_i-h_i <(4/21)u+sqrt(u(1-u))
      <=2/21+sqrt457/42<127/210<61/100.           (10)

The final seal is sqrt457<107/5.

### 2b. Exactly the near4 packets and far-end3 packets can be weak

The sole near3 exception left by the preceding inventory is the
literal near4/simple3 profile. If its helper-side parameter m>=8,
(6) gives sigma/sqrt(delta)>38/25. For m7 it exceeds3/2, while the
left side of(7) is<3/2 on q in[3,25/8]: after clearing denominators,
-3400q^2+13237q-7442 is positive at the two endpoints. Therefore a
weak near4 can only have i=5,6 with helper i-1, or i=16,17 with
helper i+1. Its ACTUAL mate is respectively4,5 or17,18 and has
normalized height greater than sqrt3/2>43/50, by the pinned mate lemma.

All other weak central opposites have EXACTLY

    v=gamma*omega_i+(1/2)*omega_e, e=1 or21,     (11)

with native finite m5 or m>=7, or a dotted own contact. No additional
white contact is being suppressed.

## 3. Remove near4 weak packets and record nonadjacency

The complete adjacent far-packet input states that two adjacent
central indices of form(11) cannot both be weak. Its native cap
product proof treats both far-end assignments: same-end caps must be
orthogonal; opposite-end caps have product-1/2. The possible4 product
would form the explicit affine F4 path3,3,4,3. Opposite-end product0
is excluded by the spatial determinant, except two5 caps, which the
pinned far5 packing excludes. Positive acute mixed replacement then
forces weighted losses<21/40 or31/80, contradicting two losses>=.61.
This exact statement is imported by the source/hash above, not assumed
for profiles outside(11).

Suppose a LEFT weak near4 packet occurs. It has an actual mate j=4
or5 with h_j/q_j>43/50. If a RIGHT weak near4 occurs, its mate k=17
or18 also exceeds43/50. The dual-direction correlation is<1/3, so
Cauchy gives a contradiction. If there is no right weak near4, at
least one of16,17 is strong: otherwise both would be adjacent weak
far packets. Its normalized height exceeds39/50. For j<=5,k>=16,

    (omega_j/q_j).(omega_k/q_k)<=30/sqrt8160<1/3.

But (43/50+39/50)^2=(41/25)^2>8/3, whereas the squared norm of
the sum of these unit directions is<8/3. This again contradicts
||x||<1. The right-to-left argument is identical. Hence no near4
weak packet exists. ALL weak central indices now have form(11),
and no two are adjacent.

## 4. A positive-dual test excludes a weak endpoint5

Suppose5 is weak. Then6 is strong. Concavity and nonadjacency give

    h_6>=a=sqrt(96/11)-.61,
    h_11>=b=sqrt(120/11)-.61,

because either11 is strong or both10,12 are strong. At least one
of16,17 is strong. Treat16 first, so h_16>=a.

The minimum energy vector x0 for the three old constraints6,11,16
has heights exactly a,b,a. It is linear between those nodes and0,22.
Its squared norm and positive old dual coefficients satisfy

    E=||x0||^2=a^2/6+(b-a)^2/5>47/50,
    c6=c16=(11a-6b)/60>4/25,
    c11=(b-a)/5>0.                               (12)

Useful rational seals, obtained simply by squaring, are

    2.344<a<2.3443, 2.6928<b<2.693,
    2.169<d=sqrt(85/11)-.61<2.17.                (13)

The same box proves (12). Put K=5a/6, L=a/12. For a far helper at21,
all pressure vectors below have old projection

    proj(u)=(5gamma/6)omega_6+(1/12)omega_16,
    x0.u=gamma*K+L.                              (14)

Write T for their actual lower bound x.u>=T, Delta=T-x0.u, and
S=||u-proj(u)||^2. For any mu>0 with nonnegative old coefficients,

    W=x0+mu*(u-proj(u)),
    x.W>=E+mu*Delta, ||W||^2=E+mu^2*S.           (15)

Thus (E+mu*Delta)^2>E+mu^2*S contradicts ||x||<1. The positivity
condition is c6-mu*5gamma/6>=0 and c16-mu/12>=0; c11 is unchanged.
It is checked below. No equality-constrained energy is substituted
for the required inequality problem.

### 4a. Complete pressure table for the helper at21

Here gamma ranges over the actual native finite labels or dotted
contacts allowed in(11), and weakness gives gamma<27/25. The latter
edge estimate, including exact rational seals, is part of the pinned
adjacent-packet proof. All vectors below retain the FULL helper.

|Own contact|Actual pressure vector u|T|S|
|---|---|---|---|
|5|gamma(omega4+omega6)/2+omega19/6|H(A,2/3-gamma^2)|2gamma^2/3+1/12|
|7,8|gamma(omega4+omega6)/2+omega20/4|A-sqrt(gamma^2-3/4)|2gamma^2/3+1/6|
|finite m>=9|gamma(omega4+omega6)/2+omega21/2|H(A,1-gamma^2)|2gamma^2/3+5/12|
|dotted|(gamma-gamma^-1)omega5+(omega4+omega6)/(2gamma)+omega21/2|A-sqrt(1-gamma^-2)|5gamma^2/3+gamma^-2-19/12|

For5, (2) uses old J={n5,n21,n20}: its actual finite diagram with
c is H4, and the old Gram is A1+A2. For7,8, (3) uses the TWO
orthogonal actual contacts J={n5,n21}, giving delta=3/4-gamma^2<0.
For finite m>=9 use(2) with J={n5}; for dotted use(4). Conversion
to the displayed omega-vectors uses x.n_l from(1), not an assumed
new lattice or idealized face. In every row

    A^2=(85/11)gamma^2+(5/11)gamma-23/44.         (16)

The complete lower deficits and chosen multipliers are

|Own contact|Delta lower bound|S upper bound|mu|
|---|---|---|---|
|5|3/10|3/4|1/5|
|7|1/4|3/4|1/5|
|8|9/40|3/4|1/5|
|finite m>=9|13/50|13/12|24/125|
|dotted|3/10|5/4|1/6|

Every inequality for Delta is STRICT. The following are complete
hand seals for the table; they do not use sampled parameter values.

* For5, gamma=(phi/2) in(.809,.81), gamma^2>.6545, so A>2.21,
  sqrt(2/3-gamma^2)<1/9, and gamma*K+L<1.78. Since H>A-sqrt(delta),
  Delta>.3. Also gamma^2<2/3 gives the stated S bound.
* For7, .9<gamma<.902, so A>2.47, sqrt(gamma^2-.75)<.253,
  gamma*K+L<1.96. Thus Delta>.257>.25, and S<.75. The endpoints
  .9,.902 are sealed by u^3-u^2-2u+1 at u=1.8,1.804 and its
  monotonicity above sqrt3.
* For8, 1.414<sqrt2<99/70 gives gamma^2>.8535, gamma>.923,
  gamma<.9239. Equation(16) gives A>2.5475. Also
  sqrt(gamma^2-.75)<.322 and gamma*K+L<2.0003, the last from
  a<2.3443 and(5*.9239/6+1/12)=3413/4000. Consequently
  Delta>.2252>.225. The S bound follows from gamma^2<7/8.
* For finite m>=9, gamma>939/1000: the increasing polynomial
  8gamma^3-6gamma-1 vanishes at cos(pi/9) and is negative at .939.
  The deficit H(A,1-gamma^2)-gamma*K-L is increasing for gamma>=.939.
  Indeed A'>sqrt(85/11)>11/4, H_A>A/(A+sqrt(delta))>20/23,
  H_delta<0 and delta'<0, while K<2. At gamma=.939, exact squaring
  gives A>2.5917, sqrt(A^2+delta)>2.6144, sqrt(delta)<.344.
  Hence H>2.29, since2.5917*2.6144>2.29*(2.6144+.344), whereas
  gamma*K+L<2.03. This proves Delta>.26 for the ENTIRE interval.
* For dotted gamma in(1,27/25), the deficit
  A-sqrt(1-gamma^-2)-gamma*K-L is decreasing: A'<3, K>1.9, while
  the derivative of sqrt(1-gamma^-2) exceeds(3/4)/(2/5)=15/8.
  At gamma=27/25, A>2.99, sqrt(1-gamma^-2)<.38 and
  gamma*K+L<a*(59/60)<2.306. Thus Delta>.304>.3. Its variance S
  is increasing in gamma and its exact value at27/25 is<5/4.

The finite5 radical seals use phi>809/500 and phi<81/50.
The finite8 S and root seals use only the rational sqrt2 interval
already displayed. In the dotted derivative argument, A' decreases
and A'(1)<3; this follows directly from(16).

All dual weights in(15) are positive. For5/7/8, mu*(5gamma/6) is
respectively<.135,.151,.154. For finite m>=9 it is at most4/25,
strictly below c6. For dotted it is<3/20. The right helper correction
is at most1/60, also below c16.

Substitute E=47/50 and the deficit/variance bounds in(15); the
difference is increasing with E. Every row is positive. The two
tight rows are

    (47/50+(1/5)*(9/40))^2-(47/50+(1/25)*(3/4))
      =9/40000>0,
    (47/50+(24/125)*(13/50))^2
       -(47/50+(24/125)^2*(13/12))
      =219/39062500>0.                           (17)

The other three rows have larger margins by direct rational arithmetic.

### 4b. The other strong endpoint and the other helper assignment

If17 is strong instead of16, use old heights a,b,d at6,11,17.
The old energy is

    E=a^2/12+(b-a)^2/10+(b-d)^2/12+d^2/10>24/25.

The old coefficients at6 and17 exceed4/25 and the middle coefficient
is positive. The box(13) proves these assertions by positive partial
derivatives or direct substitution. For a helper at21, its old
projection coefficient becomes1/10 atomega17, and

    x0.u=gamma*K+d/10,
    d/10-a/12<1/40.

Every variance in the table DECREASES by1/60. Therefore retain the
same mu and upper S, replace each deficit by its table lower bound
minus1/40, and replace E by24/25. All rational differences in(15)
are positive (the finite8 row is already1-(.96+.03)=.01).
The helper coefficient is at most.02, so dual positivity persists.

For a helper at1 instead of21, reverse the helper term in each row:
omega19/6 becomesomega3/6, omega20/4 becomesomega2/4, andomega21/2
becomesomega1/2. Its old projection joins theomega6 coefficient:

    proj(u)=(5gamma/6+1/12)omega6.

Compared with the helper21/strong16 reference, x0.u is the SAME,
whereas A increases: its squared value gains12gamma/11. The two
time coefficients have sum<25/4 throughout f<gamma<27/25, since
A_right<3 and A_left<13/4. Thus

    A_left-A_right>192/1375.

All pressure targets increase by more than1/10: for T=A-sqrt(h)
this is immediate; for H(A,delta) use H_A>A/(A+sqrt(delta))>4/5.
The deficit is therefore at least the reference table bound plus.1.
The variance gains gamma/3. All finite rows consequently have
Delta>13/40,S<3/2; choose mu=1/6. The newomega6 coefficient is at
most11/72<4/25. All dotted rows have Delta>2/5,S<5/3; choose
mu=4/25, whose old correction is at most59/375<4/25.
At E=47/50 both rational differences in(15) are strictly positive.
For the strong17 alternative, left-helper x0.u is unchanged and E
is larger, so the same two certificates work. This treats BOTH
helper assignments and BOTH old strong-node alternatives completely.

We have contradicted a weak5 in every finite/native and dotted case.
Chain reversal also excludes a weak17.

## 5. Endpoints strong and the final exact energy loss bound

Now5 and17 are strong; weak central sites inside6,...,16 are
nonadjacent, at most6 in number. For the full set of central floors
g_5,...,g_17, let Efull be the energy of the piecewise linear concave
curve through them and through(0,0),(22,0). It is the minimum energy
under those lower constraints, because its old-normal coefficients
are nonnegative. Exact symmetry gives

    Efull=g_5^2/5+sum_(i=5)^10(q_(i+1)-q_i)^2
      >1007/1000.                                (18)

For a transparent rational seal, successive lower/upper bounds for
q_5,...,q_11 are

    (2.7797,2.7799), (2.9541,2.9543),
    (3.0895,3.0897), (3.1908,3.1910),
    (3.2612,3.2614), (3.3028,3.3030),
    (3.3166,3.3167).

Every endpoint is checked by squaring against the displayed rational
D_i. They lower-bound(18) by exactly1.007192468.

Deleting one isolated internal floor replaces its two adjacent
segments by their chord. The energy decrease is

    [(g_i-g_(i-1))-(g_(i+1)-g_i)]^2/4<1/1600.

Indeed the continuous q(t)=sqrt(2t(22-t)/22) has q''(t)=-1/q(t)^3,
and q>11/4 on[5,17], so its discrete second difference has absolute
value<1/20. Nonadjacent deletions have disjoint segment pairs, so
their losses add. The surviving polygonal curve remains concave and
has nonnegative dual weights; it is still the minimum energy for
the remaining lower constraints. Therefore the ACTUAL vector obeys

    ||x||^2 >= Efull-6/1600
      >1007/1000-3/800=4013/4000>1,

contrary to(1). This proves the candidate theorem, subject to the
independent review requested at the beginning.

## Scope boundary

The conclusion concerns the WHOLE A21 farthest vertex only. It does
not assert an A21 component theorem with arbitrary spectators, an
exclusion of every larger A-rank, or a new compact/general dimension
ceiling. All central finite endpoint roles, zero-own-contact cases,
near4 packets, finite far-m/simple3 cases and dotted cases used in the
assembly are explicit above. No arithmetic field list or numerical
optimizer replaces any of them.
