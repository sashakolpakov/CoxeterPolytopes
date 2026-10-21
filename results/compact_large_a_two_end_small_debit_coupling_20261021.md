# Large A components: two-end energy and small-debit couplings

21 October 2026. Complete supporting written-proof candidate, held for
full root and outside-author review. Manual mathematics only: no
enumeration, optimization, scientific, Lean or TeX job, manuscript or
canonical status edit. **No arithmeticity hypothesis is used.**

## Statements and input boundary

At a vertex farthest from a strict interior point of a compact Coxeter
polytope, use the SAME actual vector x, ||x||^2<1, and its orthogonal
component projections. The following conclusions hold with arbitrary
finite spectators:

1. An A18 component has actual energy greater than6/7.
2. An A19 component cannot coexist with either isolated rank-two7 defect.
3. An A20 component has actual energy greater than .92274 and cannot
   coexist with an isolated short A1 of type I. In particular it cannot
   coexist with any single A1 component, ordinary or I.

Consequently the dimension21 necessary forms A18+I+T2, A19+7, and
A20+A1 are all excluded. The A15/A16 couplings and whole D21 are not
settled here, so the general compact ceiling stays21.

All paths in the table are under results/. Their transitive geometric
and forest-accounting premises are retained.

| Input | SHA256 |
| --- | --- |
| compact_a19_a20_central_window_actual_six_sevenths_20261020.md | `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` |
| compact_an_large_components_arbitrary_spectator_transfer_20261020.md | `8d3b999761be7dd46e7990973b3ec4d4a746d4d682b044d53695f4aae405e3e5` |
| compact_whole_an_geometric_a21_plateau_transfer_20261020.md | `8a58b7747e1a3d8c9affe7df17038781a62983ddd70ccb670502c6f40b1f9cee` |
| compact_a21_adjacent_weak_far_m_simple3_packet_20261020.md | `1c4437c98d916d558081572dd2ecc12bc75eda1328a53633fe111120ec63dcb7` |
| compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md | `8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b` |
| arithmetic_large_a_central_cap_packing_20261021.md, ONLY unconditional Section1 inventory scalars and Section6 I barrier | `8599baeda0f1cfae41bf3e4f27a89d7d999a49a3284caa8bdd559ea1ba945ee0` |
| compact_classical_rank21_isolated_i_actual_terminal_obstruction_candidate_20261021.md, ONLY universal I parameters/height | `8183924e37e5153a13a141daa6b7b2af44e3998b4008695ba3b1a0eef1d29bf2` |
| compact_a20_small_debit_coupling_exclusion_20261021.md, ONLY universal typed7 heights | `6ed1ed3c81f5126c6ebf66c07cf8fc3338db3b64352d8985af09db74c6ee9721` |
| compact_dimension21_forest_frontier_and_small_debit_continuation_20261021.md, ONLY accepted necessary forest/accounting boundary | `2edf213aa7bb25d7947a9d4712fdda28f95e9546a1af91e8f4cf3613c244e5bf` |

The arithmetic source's Sections2--5, its field locking and its
at-most-one-weak-index theorem are NOT inputs. The present proof allows
several nonadjacent weak indices and treats their complete energy loss.
The ordinary nonshort-A1 cost>.1, ordinary T2 cost>1/7, isolated7
cost>.09, I isolation and debit compatibility are the accepted universal
interfaces. There is no intrinsic I2 bound7/32 and no claim that all
ordinary total rank-three forests cost>1/4. Finite spherical
classification, compact acute simplicity and the actual native-wall gap
remain the standard external geometric inputs.

## 1. One weak-window inventory at N19,20,21

For A_(N-1), N in{19,20,21}, use unit old normals and their spatial duals

    D_i=omega_i^2=2i(N-i)/N,
    C_ij=omega_i.omega_j=2i(N-j)/N, i<=j,
    h_i=x_A.omega_i, h_0=h_N=0,
    e_A=(1/2)sum_(i=0..N-1)(h_(i+1)-h_i)^2<1.

The h_i are strictly concave. Index i in5..N-5 is weak when
h_i<=sqrt(D_i)-.61, otherwise strong. The complete weak inventory is:

    far helper: v=gamma*omega_i+omega_e/2,
                e=1 orN-1, f=101/125<gamma<g=27/25;
    near4:      v=omega_i/sqrt2+omega_(i-1)/2,
                i5 or6, or the reflected profile.

Every weak cap has Q=0, meaning no contact to any spectator. The own
contact in the far-helper row is native5, finite>=7, or dotted. The
near4 row is OWN4/WHITE SIMPLE3; its actual mate is respectively4 or5
or the reflection and has normalized demand>sqrt3/2>43/50.

For N20/21 this is exactly Sections1--2 of the first input. Those
inventory arguments use only e_A<1, although that source later works
under e_A<=6/7. For N19 the retained paths still have rank>=4 and
the two-white connected rank is18>8, so the same finite list,
zero-own-contact repair and positive foreign clearing apply. The
needed scalar extension is explicit:

    D_i>=140/19, sqrt(D_i)>2.7, sqrt(D_i-1)>2.5;
    foreign/native>=8 loss<95/156<.61;
    u<1 clearing loss<=(5+sqrt701)/52<.61;
    sqrt(D_i-100/81)>2.47,
    native7 loss<10000/41877+.36<.61;
    sqrt(D_i)+sqrt(D_i-25/16)>19/4 for isolated5.

All nonexceptional finite-white edge estimates exceed761/500; the
target F(q)/q is below this value on2.5<=q<=31/8. Clearing Q>0
costs at least1/4 and makes both exceptional roles strong, exactly as
in the spectator-transfer proof. Thus Q0 is proved, not assumed.
The m>=8 near4 role is strong; m7 has ratio>3/2, while F(q)/q
is increasing for q>=2.5 and is<3/2 at q=25/8. This covers N19's
q7<3 and leaves only the near4 locations just stated. No arithmetic
conjugate argument is used in any of these claims.

## 2. Literal weak helpers are nonadjacent also at N19

For N20/21 this is the actual two-cap proof in Section3 of the first
input, valid before its contrary energy bound is used. Here is the
missing N19 extension, including both opposite-helper orientations.

Let j=i+1 be central. Put A=gamma_i, B=gamma_j and R=c_i.c_j for
the distinct actual caps. The general-N weighted AM--GM identities give

    same helper:
    R>=-AB(1-1/N) +/- (A-B)/N +(N-1)/(2N)
                              +(N+1)(A/B+B/A)/(4N);
    outward opposite helpers:
    R>=1/2-AB+(A-1)(B-1)/N;
    reversed opposite helpers:
    R>=1/2-AB+(A+1)(B+1)/N.

At N19 these imply R>-1/2 in the first case and R>-7/10>
-1/sqrt2 in the other two; for example the latter lower bound is
at least1/2-g^2-(1-f)(g-1)/19>-7/10. The actual native gap
therefore leaves R0 for a common helper, or R0/-1/2 for opposite helpers.

For the outward orientation the retained Schur parameters still obey
eta_i,eta_j<4/5. At N19, C_ij>=130/19, so

    (D_i-4/5)(D_j-4/5)-C_ij^2
      =(2/5)C_ij-24/25+8/(5N)>=884/475>0.

This excludes R0 exactly as in the accepted actual packet proof.

For the reversed orientation put l=i,s=N-j. Then l+s18 and l,s>=5.
The exact determinant is

    NF=4ls A^2 B^2-(l+1)(s+1)(A^2+B^2)-2AB
                         -(l+1)A-(s+1)B+(N+2)/4.       (1)

It increases strictly in A and B on[4/5,g]: the bracket in its
A derivative is at least(39/25)ls-N>=82.4, and the subtracted
term2B+(l+1) is at most16.16. At least one gamma is>=.9 unless
both are native5. In the first case use the actual lower floor
f=.808, rather than the older loose .8. Concavity in l reduces
the minimum to l5 or13, at(A,B)=(f,.9) or its reverse. The two
endpoint values are

    .9601824=300057/312500,   1.6961824.

Hence NF>0, excluding R0. If both are native5, the independently
proved far5 distance-one prohibition, valid for every N>=15 and both
helper assignments, excludes the pair without arithmeticity. Thus
opposite helpers have R=-1/2 in every remaining case.

The mixed finite replacement is now ACTUAL exactly as in the pinned
packet proof. Its two-cap inverse Schur matrix is entrywise positive;
appended spectators form an orthogonal positive block. Positive acute
projection is feasible against every unknown wall. Its displacement
is Z=W_i omega_i+W_j omega_j with W_i,W_j>0. All the accepted
pressure constants remain valid: the helper parameters are>=5,
delta in(1/2,3/5], gamma>.8, and the time coefficients exceed2.
For the last point at N19,

    a_i^2> (16/25)(140/19)+(4/5)(10/19)+9/19-1
           =438/95>4.

Consequently sqrt(1+Dcal)/(W_i+W_j)<7/16 for a common helper,
or<3/10 for opposite helpers. The displacement correction is

    (W_i sqrt(D_i)+W_j sqrt(D_j))-||Z||
       <(W_i+W_j)/[2(sqrt(140/19)+sqrt(130/19))]
       <(W_i+W_j)/10,

using sqrt(140/19)>2.7 and sqrt(130/19)>2.6. The actual mixed
point therefore makes the weighted average of the two losses<43/80
or<2/5, respectively. Both are<.61, contradicting two weak indices.
This completes nonadjacency at all three N, before removing near4.

## 3. Removing near4 with the correct energy or rank hypothesis

If a left near4 occurs, its actual mate has index j<=5 and normalized
demand>.86. At the opposite end either there is a near4 mate, or
helper nonadjacency gives a strong k in{N-6,N-5}. Every such strong
demand is>.77 because q_k>2.7. For normalized dual directions,

    rho=sqrt(j(N-k)/((N-j)k)).                       (2)

At N19 this is at most sqrt(15/91)<1/2. Cauchy excludes a near4
whenever e_A<=6/7, since1.63^2>3*(6/7). At N20 it is at most
1/sqrt7<.4, so Cauchy excludes near4 whenever e_A<=.91, since
1.63^2>.91*(14/5), with margin1089/10000. Opposite near4 mates
have still larger demands and smaller correlation. Reflection treats
a right near4. These arguments use the actual mate, not a formal cap.

At N21 a sharper rank-specific comparison removes near4 under ONLY
e_A<1. If the opposite anchor is strong k15, then q15^2=60/7>
(61/21)^2, giving normalized demand>.79. Also rho<=1/sqrt8<.354.
Thus the demand sum exceeds1.65, whose square2.7225 is greater than
2.708, the upper bound on the squared norm of the sum direction.
If the strong anchor is k16, its demand>.77 and rho<=5/16;
1.63^2>2+2(5/16)=2.625 again contradicts e_A<1. If instead
there is an opposite near4, its mate has k>=16 and demand>.86,
so the latter comparison is stronger. Thus at N21 all central weak
indices are nonadjacent pure helpers with arbitrary spectators.

## 4. One two-end energy certificate, with all weak sites retained

Let H be a concave symmetric rational height sequence below the strong
floors, linear from0 to5 and fromN-5 toN. Its coordinate reference
is w_k=(H_k-H_(k-1))/sqrt2, with zero coordinate sum. At available
strict bends, summation by parts gives x_A.w>=||w||^2; hence the
actual energy is at least that reference energy. Strictly smaller
floors make the final energy bound strict.

If index5 is weak, replace the entire left end by the line0 to6;
then6 is strong by nonadjacency. Make the reflected replacement if
N-5 is weak. For each remaining internal weak j, replace H_j by
the average of its neighbours. Concavity and all other floors are
preserved. The exact energy loss is

    (H_(j-1)-2H_j+H_(j+1))^2/4.                     (3)

Since weak sites are nonadjacent, these losses add: their affected
edge pairs are disjoint. This handles EVERY permitted set of weak
indices, not only a selected or arithmetically packed subset.

Use the following three sequences, reflected as H_(N-i)=H_i:

| N | H5 | H6 | H7 | H8 | H9 | H10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 19 | 2.104 | 2.255 | 2.363 | 2.433 | 2.467 | 2.467 |
| 20 | 2.128 | 2.288 | 2.406 | 2.488 | 2.536 | 2.552 |
| 21 | 2.150 | 2.317 | 2.445 | 2.537 | 2.597 | 2.626 |

At N21 set H11=H10. At N19 H10=H9 is already reflection.
All entries are strictly below sqrt(D_i)-.61 by positive squaring,
and their successive differences decrease. For instance the tight
N21 i7 comparison has28/3-3.055^2=925/3000000>0.

For k=0,1,2 weak end5 indices, the exact initial reference energies
before internal deletions, and upper bounds on ALL internal losses,
are as follows. The repeating decimals merely abbreviate the exact
square formulas for these rational references.

| N | k | Initial energy | Maximum loss from internal weak indices |
| --- | ---: | ---: | ---: |
| 19 | 0 | .9258842 | 4(.00046225)=.001849 |
| 19 | 1 | .8955541833... | 4(.00046225)=.001849 |
| 19 | 2 | .8652241666... | 3(.000361)=.001083 |
| 20 | 0 | .9544848 | 5(.000441)=.002205 |
| 20 | 1 | .9250917333... | 4(.000441)=.001764 |
| 20 | 2 | .8956986666... | 4(.000324)=.001296 |
| 21 | 0 | .981678 | 5(.00038025)=.00190125 |
| 21 | 1 | .9528575833... | 5(.00038025)=.00190125 |
| 21 | 2 | .9240371666... | 4(.000324)=.001296 |

For verification, write a=H5,b=H6 and
T=sum_(j=7..floor(N/2))(H_j-H_(j-1))^2. The three initial
energies are respectively a^2/5+(b-a)^2+T,
(1/2)[b^2/6+a^2/5+(b-a)^2+2T], and b^2/6+T.
This also applies when N is odd because the central edge has zero
difference. In the k2 cases the strong6 and N-6 leave internal
interval7..N-7, with respectively3,4,4 possible nonadjacent holes.
For k0/1 the interval lengths give the counts in the table, and the
largest quarter-squares in(3) give their per-hole bounds.

An endpoint bound x1>=z>b0 on a constant first block of length m,
mean b0, can be added by setting its first reference coordinate z
and all remaining m-1 coordinates to(m*b0-z)/(m-1). It preserves
the block sum and raises the reference energy by

    [m/(m-1)](z-b0)^2.                              (4)

Its new dual at1 is positive. The dual at m stays nonnegative if
the remaining block coordinate is at least the next coordinate.
All references here pass that junction check even at z1: at N20
the first5 remaining mean exceeds .126 while the next coordinate
is below .114; the first6 remaining mean exceeds .123 while the
next is below .084. At N21 these pairs exceed .130 versus .119,
and .127 versus .091. Deleting a neighbouring internal floor only
decreases the next coordinate. Opposite-end changes on disjoint
blocks add, keeping all later prefix sums and the total sum unchanged.
If z<=b0 no modification is needed, giving the squared positive part.

## 5. Consumer A18+I+T2 in dimension21

Assume contrary e_A<=6/7. Sections1--3 give the complete pure-helper
weak window at N19. Every row of its energy table exceeds6/7;
the lowest is

    .8652241666...-.001083>.86414>6/7.

Contradiction. This proves the unconditional component theorem
e_A>6/7 with arbitrary spectators. In particular an ordinary
irreducible rank-two T2 already costs>1/7 in orthogonal coordinates,
so A18+I+T2 is impossible even without using I's positive energy.

## 6. Consumer A19 plus either rank-two7 defect

The actual defect energy>.09 forces e_A<.91. Sections1--3 therefore
give pure-helper coverage and nonadjacency at N20. If k0 or1 end5
indices are weak, the respective table rows give

    e_A>.9522798 or e_A>.9233277333...,

and adding>.09 already contradicts ||x||^2<1.

For k2 both6 and14 are strong. The table gives e_A>.8944026666... .
Both typed7 defects supply every other actual direction the common
normalized projection-square bound

    theta>(11-4sqrt7)/3>11/80;

the OLD A2/NEW7 role in fact gives theta>.2. This universal height
argument retains the entire old7 block at other edges and has no
restriction to an A20 rank. The inequality11/80 follows from
sqrt7<847/320, with positive square difference609/102400.
At the two A19 endpoints,

    x1^2,x20^2>(19/20)(11/80)=209/1600>.36^2.

Thus x1>.36 and -x20>.36. The first and last6 blocks have means
H6/(6sqrt2)<.27 in absolute value. Applying(4) to these DISJOINT
blocks adds more than

    2*(6/5)(.36-.27)^2=.01944.

Consequently the SAME A-plus-defect energy is

    >.8944026666...+.01944+.09>1.0038426.

This closes both typed7 alternatives without assuming any intrinsic
stronger bound for arbitrary ordinary dihedral components.

## 7. Consumer A20 plus an ordinary or isolated-I A1

Section3 removed near4 at N21 using only e_A<1. Hence its energy
table applies with arbitrary spectators and yields

    e_A>.97977675, .9509563333..., or .9227411666...

for k0,1,2. In particular e_A>.92274. Any nonshort A1 has
energy>.1 and is immediately impossible alongside A20.

For I, use its universal parameters0<u<1/3,t>0 with

    e_I>=u^2, (1-u^2)(1-t^2)=8ut,
    t+(3/2)u>.767.

Every OTHER actual edge has normalized projection>=t; this holds
at BOTH A20 endpoints, so x1,-x21>=a*t, a=sqrt(20/21)>.975.
The defining relation implies t<1. Therefore every required endpoint
adjustment has z=a*t<a<1 and is legal by Section4.

Count I only once. For each endpoint the block length is5 or6,
so C=m/(m-1)>=6/5. Let b_bound be a common upper bound for their
absolute means: b_bound=.305 if k0 or1, and .274 if k2. By two
disjoint adjustments, their total energy gain is at least
2*(6/5)*max(a*t-b_bound,0)^2. Weighted Cauchy and the STRICT
barrier give

    e_I + both A-end gains
      >(a*.767-b_bound)^2/[(9/4)*a^2+5/12].          (5)

Indeed y=max(a*t-b_bound,0) obeys
y+(3/2)a*u>a*.767-b_bound, and the minimized cost is
u^2+(12/5)y^2. There is no sum of two separate I charges.
The denominator in(5) is215/84<2.56. When k0/1 its numerator
exceeds(.767*.975-.305)^2=.442825^2; the ratio is>.076.
When k2 the numerator exceeds(.767*.975-.274)^2=.473825^2;
the ratio is>.0875. Thus the worst total energy is greater than

    .9227411666...+.0875>1.01024.

The k0/1 cases are stronger, since even .9509563333...+.076>1.
This excludes A20 with isolated I, also with arbitrary further
spectators. In a dimension21 forest its sole A1 spectator is either
ordinary nonshort or I; a type-J pair needs two A1 components.
Thus every A20+A1 alternative is removed.

## Scope of the continuation

The proof is one actual two-end energy construction, with an N19
geometric packet extension and the full nonadjacent weak-window loss
retained. It uses no arithmetic packing and assumes no desired
component energy or forest exclusion. All endpoint heights concern
one common x, all component charges are orthogonal, and every added
endpoint gain is paid inside a single modified reference.

Subject to complete independent review, the dimension21 frontier
loses A18+I+T2, A19+either7, and A20+A1. Whole B21 was excluded
separately. The remaining general dimension21 alternatives include
the listed A15/A16 debit couplings and whole D21. No ceiling20,
new manuscript, checker execution or end-to-end Lean theorem is
claimed; the finished six-paper snapshots remain unchanged.
