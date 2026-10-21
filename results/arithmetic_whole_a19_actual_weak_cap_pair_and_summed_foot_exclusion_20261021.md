# Whole arithmetic A19: actual weak-cap pairs and summed feet

21 October 2026. Complete written-proof candidate; root and full
outside-author review are pending. Source/manual reasoning only: no
scientific, enumeration, optimization, Lean or TeX job. This is one
supporting lemma for the existing arithmetic paper, not a new manuscript
or an arithmetic dimension-18 claim. The preceding dimension-19 frontier
packet `c5992464` is unchanged.

**Theorem.** A farthest vertex of a compact arithmetic first-type
hyperbolic 19-polytope cannot have whole finite diagram A19. The actual
totally real ground field is unrestricted. No parked field census is
used.

The new geometric step forbids adjacent weak *isolated dotted* caps by
an actual two-cap replacement. Finite eight caps are retained with their
actual corrected pressures. They are not incorrectly declared strong.
The resulting single summed certificate exceeds one by more than
`247852/10^9`. Near-four caps newly admitted by the smaller threshold
are also explicitly treated.

## 1. Inputs, normalization and actual-point boundaries

All paths in this table are under `results/`; their current bytes were
hashed, without executing their mathematical checks.

|Input|SHA-256|
|---|---|
|`arithmetic_large_a_central_cap_packing_20261021.md`|`8599baeda0f1cfae41bf3e4f27a89d7d999a49a3284caa8bdd559ea1ba945ee0`|
|`arithmetic_uniform_classical_weak_helper_obstruction_20261021.md`|`535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905`|
|`arithmetic_full_an_isolated_finite_continuation_packing_20261020.md`|`9ab34f073ea95bed2040bafe6b8573331721f9b1e75d65172f814637da046810`|
|`compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md`|`3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb`|
|`compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`|`3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`|
|`compact_an_projection_pressure_and_weak_cap_packing_review_20261020.md`|`8b0a59493f1fa5dc7081e5e70beff2839ac5918bfcd50e089e05ae9b5546e43b`|
|`whole_a20_summed_isolated_feet_technical_input_20261021.md`|`8d93c23c9bf268813334dd12fef09f37f523eb22208971847cf8b04e1d659ab4`|

Use UNIT old normals `n_i=(e_i-e_(i+1))/sqrt2`, `N=20`, and their
spatial duals `omega_i`. The common farthest vector satisfies

    ||x||²<1, x.n_i>0,
    D_i=omega_i²=2i(20-i)/20,
    C_ij=omega_i.omega_j=2i(20-j)/20  (i<=j),
    h_i=x.omega_i, h_0=h_20=0,
    ||x||²=(1/2)sum_(i=0..19)(h_(i+1)-h_i)².

The whole diagram has no foreign spectator space. An actual opposite is
`c=A p-v`, with `c²=1`, `A>0`. Positive finite projection on selected
actual walls has nonnegative correction coefficients and is feasible
against every unselected actual wall by acuteness. Its full displacement
`Y` obeys

    x.Y>=H(A,delta)=A/[1+sqrt(delta/(A²+delta))]
                  >A-sqrt(delta),                    (1)

where `delta>0` is the selected residual. A feasible future negative core
of square `-d` instead gives `x.Y>=A-sqrt(d)`. These are inequalities
from actual points, including nonvertices. No abstract Gram completion
or intrinsic Coxeter property of an arbitrary face is assumed.

The declared external foundations are finite spherical classification,
compact acute simplicity, the native gap for actual distinct walls,
and the first-type cyclic-integrality/conjugate-positivity criterion.
The cited internal arguments supply the actual projection/farthest
normalization and finite-white inventory. In particular, the old
loss-.61 arithmetic packing theorem at N20 is an input: every such weak
central cap is golden far-helper at 5,6,7 or its reflection, with at most
one weak index. We do not transfer a whole-A20 conclusion to A19.

## 2. Fresh complete loss-.55 inventory

Work on `5<=i<=15`; both retained old paths have rank at least four.
Their finite-white list consists exactly of isolated caps, one admissible
white3 (end, D or E position), one end-white4, or two end-white3 contacts.
H/F completions are too small. A two-path completion has rank19>8, and
the old arm lengths exclude a D branch or internal4; it is A19.

Own contact is positive. With zero own contact, one far/next-to-end3
would give a positive definite full A20/D20 star in positive index19;
an interior3 contains the affine arms `(1,2,5)` or `(1,3,3)`; an end4
gives B20 and an interior4 contains the affine squared-cosine star
`1/2,1/4,1/4`; two end3 contacts give an affine cycle. All are impossible.
These checks use the actual full old A19, not only its retained paths.

For `q=sqrt(D_i)` the sufficient edge ratio is

    R(q)=20/11-(q-.55)/[q(2q-.55)]<419/250=1.676.    (2)

It increases for `q>=sqrt(15/2)`. Since `q<=sqrt10<3.17`, it suffices
to check (2) at `3.17`; positive denominator clearing gives
`26200*2750-391*183543=284687>0`.

For a white3 at node `j` counted from the far end of a retained path
of rank `m-1`, `m>=5`, its Schur data are

    delta=1-j(m-j)/(2m), sigma=gamma+j/(2m).

The old affine own3/4/6 exclusions for non-near roles and own3/6
exclusions for near roles remain principal actual diagrams. Rechecking
all ratios against (2) gives:

* D position: `delta=2/m`, `sigma>=gamma+1/m`. At `m>=6`, its
  minimum `(phi/2+1/6)/sqrt(1/3)` exceeds1.676. At m5 own>=7 is
  strong, as is the nearer j3 position. Only own5/helper j2 at m5
  remains for Section3.
* E positions: `delta<=1/4`, `j/(2m)>=3/16`, and
  `gamma>=phi/2`, so the ratio exceeds1.99.
* Near-end white3: own>=7 has ratio greater than
  `(.9+.4)/sqrt(.6)>1.676`. Own5 is strong by the ratio once `m>=9`;
  its m5..8 roles are treated below. Own4 is treated in Section4.
* Far-white4 has own>=5, `delta<=1/5`, and ratio at least
  `(phi/2)/sqrt(1/5)>1.8`. Near-white4 has ratio greater than
  `(.5+.7*.8)/sqrt(1/5)>2`. Both are strong.
* Two white3 contacts have ratio at least `q/2+1/q>1.73`, hence
  are strong.

These exhaust the nonisolated roles, apart from the far-end white3
helpers treated in Section5. Reflection retains both end orientations.
There is no unclassified spectator or two-component bridge in a WHOLE
A19 diagram.

For isolated finite caps, the all-field input gives exactly labels5,8,12
at these indices (`D_i/2>2`), not a preselected ground-field list.
Label5 has actual H4 loss less than `25/76+5/24<.55`. Label12 has loss
less than `.19+tan(pi/12)<.458`, using
`q+sqrt(q²-1)>1/.19`. Thus only finite8 and dotted isolated caps need
retention. Their precise actual inequalities appear in Section7.

## 3. Own-five near and D-position roles

Write `phi=(1+sqrt5)/2`, `gamma=phi/2<.81`. For the near role
`v=gamma*omega_i+omega_(i-1)/2`, the actual right H4
`{c,n_i,n_(i+1),n_(i+2)}` has positive old correction

    w=gamma*(3n_i/2+n_(i+1)+n_(i+2)/2),
    delta=(5-3phi)/8<.136².

Equation(1) gives `x.v>x.(v-w)>A-.136`. If the index were .55-weak,
`x.v<.81(q_i-.55)+q_(i-1)/2`. The time coefficients and rational
bounds are:

|i|40A²|A lower bound|upper bound on x.v under weakness|A-.136 lower bound|
|---|---|---:|---:|---:|
|5|99+195phi|3.219|3.03809|3.083|
|6|119+224phi|3.469|3.27219|3.333|
|7|135+247phi|3.655|3.44777|3.519|
|8|147+264phi|3.788|3.57319|3.652|

The A bounds follow by squaring and `phi>1.618`. For the upper bounds
use `q_4<2.530` and `q_5,q_6,q_7,q_8<2.739,2.899,3.017,3.099`.
Every row contradicts weakness. This excludes the weak profiles, not
the underlying walls irrespective of their heights. No time-square
packing premise is needed here.

The D-position remainder is `i5,v=gamma*omega5+omega2/2`, or its
reflection. Its simple3 helper normalizes the cap to square2 over the
actual field. The full Lorentzian star locks that field to Q(phi).
Its time value `40A²=71+135phi` has norm `-3599=-59*61`. This has
different square class from all three old-.61 weak-cap time norms
`-2439,-1519,-599` (classes `-271,-31,-599`). Two such K-valued
timelike residuals would have K-square time ratio and therefore a
rational-square norm ratio, impossible. Hence every central index is
old-.61 strong in this case.

Use old rational heights

    L5..10=(2.128,2.288,2.406,2.488,2.536,2.552),
    L_(20-i)=L_i, linear tails;  E_old=.9544848.

All these are strictly below their old floors. Right-H4 projection
gives `Y=(3gamma/4)omega4+(gamma/4)omega8+omega2/2` and demand
`T>A-.136>2.553`, since `A>2.689`. Its reference product is below
`1.963628`, giving deficit `Delta>.589`. Projection onto the available
old dual span gives

    P=(.6gamma+.2)omega5+(gamma/4)omega8,
    S=||Y-P||²=(8/5)a²+(4/5)a+3/5<1.677,
    a=3gamma/4.

Take `mu=.08`. Payments at5 and8 are less than `.05488` and`.0162`,
below their old coefficients `.1328` and`.017`. Thus the single vector
`W=w_old+mu(Y-P)` has nonnegative remaining old coefficients, and

    ||x||²>=E_old+2mu*Delta-mu²*S
            >.9544848+.09424-.0107328>1.

Both orientations of the D-position remainder are impossible.

## 4. Newly admitted near-four profiles

For own4 and the neighboring white3 at `i-1`, let `m=i` be the
helper-side parameter. Its square-two K-normalized star is over
Q(sqrt2), and locks the actual field to this field. At the nonphysical
embedding the time numerator is

    M_m'=60(m-1)-(3m²-2m+1)-2sqrt2*(m-1)(20-m).

It is concave in m, and is positive at m5 and13:
`174-120sqrt2>0`, `238-168sqrt2>0`, by positive squares. Thus the
nonphysical Schur complement `-M_m'/20` has the wrong sign for all
`5<=m<=13`. The only remaining orientations have m14 or15. After
reflection these are precisely own i6/helper7 and own i5/helper6.
Neither is omitted because it was strong at the earlier threshold.

Every old-.61 weak cap locks K to Q(phi), so no such cap can coexist
with either near-four profile. All central old floors L are available.
The actual same-wall mate is j=i+1: replacing old j by c gives the
positive definite full forest `B_(i+1)+A_(18-i)`. Actual-wall projection
makes its vertex feasible, so it is the outgoing endpoint at j. Its
edge formula is

    cosh L_j=[1/2+i/((i+1)sqrt2)]sqrt(j D_j).

For j6 the existing uniform mate calculation gives
`h6²>(3/4)D6=6.3`, hence `h6>2.509`. Use the concave reference
with active values `L5=2.128,h6=2.509,L10=2.552`, old right-hand
values thereafter, linear between active nodes and on the left tail.
Its energy is exactly

    2.128²/10+.381²/2+.043²/8+.9544848/2
      =1.002892425>1.

For j7, the displayed cosh is greater than `1.1*7.9>17/2`.
Consequently `h7²>9.1*15/19>(2.68)²`. Use active values
`L5=2.128,h7=2.68,L11=2.536`, and the old right-hand values from11
onward. Its slopes are `.4256,.276,-.036,-.048,-.082,...`, so all
dual coefficients are nonnegative, and its energy is

    2.128²/10+.276²+2*.036²+(.9544848/2-.016²/2)
      =1.0087208>1.

The first reference is also concave (successive active slopes
`.4256,.381,.01075,-.016,...`). Only active indices are constraints;
interpolated certificate heights are not asserted to be lower bounds
on unselected actual heights. Both are single-vector energy arguments.
Reflection covers the other end.

## 5. All .55-weak far helpers are golden, and are impossible

Here `v=(t/2)omega_i+omega_e/2`, `e1 or19`. A helper3 permits
square-two normalization over K, and cyclic integrality makes t an
algebraic integer. Affine exclusions leave own5, own>=7 or dotted, so
physically `t>202/125`. For `gamma=t/2>=6/5`, the edge ratio squared
is minimized at helper-side m5 and is at least `169/60>(1.676)²`.
Thus weakness forces `t<12/5`.

The nonphysical inverse Schur costs are minimized at N20,i5. For
negative and positive t respectively these are

    Q_minus(t)=7.5t²+3t+1.9,
    Q_plus(t)=7.5t²+t+1.9.

Since `Q_minus(-.8)=4.3>4`, `Q_plus(.5)=4.275>4`, every
nonidentity conjugate lies in `(-4/5,1/2)`. Every embedding of Q(t)
lifts to K. These are bounds for the minimal field, not merely for a
chosen list of embeddings.

* Degree>=5: the nonzero integral norm has absolute value less than
  `(12/5)(4/5)^4=.98304<1`.
* Degree4: the absolute norm is<2; any positive small conjugate
  makes it<`(12/5)(1/2)(4/5)²<1`. Thus the other three conjugates
  are `-a,-b,-c`, norm-1. The integral trace lies in(-1,1), so
  `t=a+b+c`. Norm and AM--GM give `t^4>=27`. The integral second
  coefficient is `-(a²+b²+c²+ab+ac+bc)`, strictly greater than
  `-6(4/5)²>-4` and at most `-(2/3)t²<=-2sqrt3<-3`, impossible.
* Degree3: the same norm bounds force two negative conjugates `-a,-b`
  and norm1. The integral trace is1, so `a+b=t-1`, `ab=1/t`.
  Reality gives `t(t-1)²>=4`, hence `t>9/4`. The integral second
  coefficient `1/t+t-t²` lies strictly between its endpoint values
  `-883/300>-3` and `-341/144<-2`, impossible.
* Degree2: the trace is1 or2. Trace2 gives nonzero integral norm
  `t(2-t)` in(-1,1). Trace1 has `t<9/5` and forces norm-1, hence
  `t=phi`. Degree1 leaves only the compactly forbidden parallel t2.

The old conjugate costs for this exact golden profile, independent of
height, restrict it to i5/6/7 with the nearer helper, or reflection.
Its time numerators are `54+105phi,63+112phi,70+117phi`, with norms
`-2439,-1519,-599`. Distinct i cannot coexist with an old-.61 weak cap
by the time-square norm test. At the same i the reflected pair violates
cyclic integrality: it would require20 to divide both `2(i²+1)` and
`2i(i+2)`, whose gcd is2. At one index the opposite is unique.
Therefore all OTHER central indices are old-.61 strong.

The actual H4 on `{c,n_i,n1,n2}` has correction
`gamma*n_i+(2/3)n1+(1/3)n2` and residual
`delta=2/3-gamma²<.111²`. Dropping its positive old correction in
(1) gives

    x1>sqrt2*(A-.111-gamma*h_i)>.6.                (3)

For i5/6/7 use respectively `A>2.365,2.470,2.546`,
`h_i<2.189,2.349,2.467`, `gamma<.81`, and `sqrt2>1.4`.
These produce lower bounds `.673274,.638834,.611422`.

Take the old reference L, omitting just i's height by its chord. The
energies at i5/6/7 are `.9250917333...,.9540438,.9541608`.
For i6/7 replace its first five equal coordinates by first coordinate
`.6` and four equal coordinates with the same total. This gains more
than `(5/4)(.6-.31)²>.105`, with all old coefficients still nonnegative:
the four equal coordinates exceed `.22` while the next is below`.12`.
For i5 do the same in the first block of six; its old coordinate is
below`.28`, the gain exceeds `(6/5)(.6-.28)²=.12288`, and the five
equal coordinates exceed`.2`, above the next coordinate. The new bend
at1 uses (3), not an assumed first-edge floor. Each resulting SINGLE
positive reference has energy>1. All far-helper weak profiles are
therefore impossible, including their reflected assignments.

Sections2--5 now leave exactly:

    h_i>sqrt(D_i)-.55, OR the actual cap is isolated finite8,
    OR it is isolated dotted.                         (4)

No cap in (4) contacts an undeclared old or foreign wall.

## 6. Adjacent internal weak dotted caps cannot coexist

Suppose adjacent `i,j=i+1` in6..14 are both in the dotted alternative.
Write `u_l=gamma_l^-2`. Their actual dotted feet imply

    h_l>T_l=sqrt(D_l-u_l)-sqrt(u_l(1-u_l)).

Rationalization bounds the loss by `a_l*u+sqrt(u(1-u))`, where
`a_l=1/[sqrt(D_l)+sqrt(D_l-1)]<.18`. If u<=1/3 this is<.532;
if u>=4/5 it decreases and is at most`.544`. Thus

    1/3<u_l<4/5,  5/4<gamma_l²<3,
    3<A_l<11/2,  |A_i-A_j|/(A_i+A_j)<5/17.          (5)

Let `C=C_ij>=7.8`. The adjacent identities are
`D_i+D_j=2C+1.9`, `D_i D_j-C²=2C`. Hence

    (D_i-u_i)(D_j-u_j)-C²>.4C-.88>0.

The actual cap product is negative. Put `k=-c_i.c_j`; native gap gives
`k>=1/2`. Equal caps are impossible because their isolated old supports
differ.

If `k<1`, the actual two-cap Gram is positive definite. Its projection
has positive `lambda=G^-1(A_i,A_j)`, displacement
`Z=lambda_i gamma_i omega_i+lambda_j gamma_j omega_j`, and
`Dcal=A.lambda`. Let `Lambda=lambda_i+lambda_j=Aplus/(1-k)`,
`W=lambda_i gamma_i+lambda_j gamma_j`,
`E0=lambda_i gamma_i sqrt(D_i)+lambda_j gamma_j sqrt(D_j)`.
Then

    Dcal/Lambda² <=1/4+25/3468,
    1/Lambda²<1/144, W²>(5/4)Lambda²,
    (sqrt(1+Dcal)/W)²<10993/52020<529/2500.

Moreover `sqrt(D_i D_j)-C<1`, so

    E0-||Z|| < W/[2(sqrt(8.4)+sqrt(7.8))]<.089W.

Actual finite projection gives `x.Z>=1+Dcal-sqrt(1+Dcal)`, and thus
`E0-x.Z<(.089+.46)W=.549W`. Two weak indices instead give
`E0-x.Z>=.55W`, contradiction.

If `k=1`, `c_i+c_j` is a nonzero feasible future null vector,
contradicting compactness. If `k>1`, it is a feasible future timelike
core of square `-2(k-1)`. With `Y=gamma_i omega_i+gamma_j omega_j`
and `W=gamma_i+gamma_j`, actual comparison gives
`x.Y>=Aplus-sqrt(2(k-1))`. AM--GM gives

    2(k-1)<=1.9gamma_i gamma_j-4
             <=(19/40)W²-4,
    2(k-1)/W²<17/120<.377²,

because W²<12. The same radial averaging estimate now gives loss
less than `(.089+.377)W<.55W`, again impossible. This exhausts ALL
native products, including disjoint caps. Weak dotted internal indices
are therefore nonadjacent. No arithmetic field restriction was used in
this section.

## 7. One summed actual-foot certificate

Use the symmetric rational reference

    H5..10=(2.1886,2.3482,2.4666,2.5483,2.5964,2.6122),
    H_(20-i)=H_i, linear tails to0 and20.

Positive squaring proves `H_i<sqrt(D_i)-.55`. Its old dual coefficients
and energy, by direct rational arithmetic, are

    beta5=.13906, beta6=.0206, beta7=.01835,
    beta8=.0168, beta9=.01615, beta10=.0158,
    beta_(20-i)=beta_i,
    w=sum beta_i omega_i,
    E=w²=sum beta_i H_i=1.006722852.                 (6)

At a dotted index the ACTUAL negative core `c+n_i/gamma` gives

    x.(omega_i-u_i n_i)>=T_i,
    T_i=sqrt(D_i-u_i)-sqrt(u_i(1-u_i)), 0<u_i<1.

At a finite8 index positive projection on the ACTUAL pair `{c,n_i}`
instead gives

    x.(omega_i-n_i)>=T8(q_i),
    T8(q)=r*sqrt((r-s)/(r+s)), r=sqrt(q²-1),
    s=tan(pi/8)=sqrt2-1<83/200.                     (7)

The own-normal coefficient here is ONE, not gamma^-2. Define J to be
exactly the indices with `h_i<=sqrt(D_i)-.55`; by (4) their caps are
finite8 or dotted. Let a_i be u_i at dotted indices and1 at finite8
indices. Set `Z=sum_J beta_i a_i n_i`, `W=w-Z`, and
`A=E+sum_J beta_i(T_i-H_i)`. Combining the actual inequalities gives
`x.W>=A`. Since old adjacent products are-1/2 and all others zero,
`Z²<=sum_J beta_i² a_i²`. Consequently

    ||x||²>=2A-W²
      >=E+sum_J [2beta_i(T_i-H_i)+(2a_i-a_i²)beta_i²]. (8)

All omitted adjacent cross terms are NONNEGATIVE. There is no claim
that the effective dual energy is itself an equality to a squared norm.

### End indices5 and15 have nonnegative error

For finite8, `T8>sqrt(q²-1)-s` gives loss less than
`.19+.415<.61`, and hence a brace greater than
`beta*(beta-.12)>0`. For dotted caps let `a=1/[q+sqrt(q²-1)]<.19`
and `beta=.13906`. The brace divided by beta is greater than

    1.1-2a u-2sqrt(u(1-u))+beta(2u-u²).

For u<=3/5 use `2u-u²>=7u/5` and
`b=a-.7beta<.0927`, `sqrt(1+b²)<1.0043`. The expression exceeds
`1.1-[b+sqrt(1+b²)]>0`. For u>=3/5 its negative part is decreasing:
its derivative is at most `.38-1/sqrt6<0`. Thus the end error is
positive for every dotted u, without an end-specific weak packing rule.

### Exact internal dotted and finite errors

For dotted internal indices, Section6 gives `1/3<u<4/5`. The loss
expression
`2a u+2sqrt(u(1-u))-beta(2u-u²)` decreases for u>=2/3, since its
derivative is at most `.36-1/sqrt2<0`. On[1/3,2/3],
`u(1-u)>=2/9`. With `b=a-beta/2`, its maximum is therefore at most
`b+sqrt(1+b²)-(2/9)beta`. This bounds the negative brace in (8).

|i, reflected|a upper bound|sqrt(1+b²) upper bound|negative brace less than|
|---|---:|---:|---:|
|6,14|.178|1.014|.001589|
|7,13|.1706|1.013|.001291|
|8,12|.1659|1.0124|.001112|
|9,11|.1632|1.012|.001027|
|10|.1623|1.0119|.000993|

The a bounds use lower bounds on q and sqrt(q²-1); all displayed
square-root comparisons are strict positive-square comparisons.

For finite8, `q-T8(q)` decreases with q. This derivative fact is
proved in the pinned singleton-foot input for `q>=sqrt2,s<=1/2`.
The following fully rational seals cover all needed positions:

    i6: T8²>(68/25)²*(461/627)>(2329/1000)²,
        q6<2899/1000, hence loss<.57;
    i7: T8²>(1423/500)²*(2431/3261)>(2457/1000)²,
        q7<3017/1000, hence loss<.56;
    i8: T8²>(1173/400)²*(1007/1339)>(1271/500)²,
        q8<3099/1000, hence loss<.557.

Thus its brace `beta²+2beta(T8-H)` is greater than-.0004 at6/14,
greater than-.00004 at7/13, and nonnegative at8..12 because
`beta>=.0158>.014`. No finite8 spacing theorem is needed.

### Complete nine-position charge

In millionths, first allow finite8 costs400,40,0,0,0,0,0,40,400
at6..14, total880. Replacing an index by a dotted error increases this
upper bound by at most

    1189,1251,1112,1027,993,1027,1112,1251,1189.

Dotted indices are nonadjacent. The elementary maximum-weight
recurrence `F_j=max(F_(j-1),F_(j-2)+weight_j)`, with zero initial
values, has the following nine successive exact values:

    1189,1251,2301,2301,3294,3328,4406,4579,5595.

This is a displayed nine-term hand calculation, not an executed table
or an assumed optimizer result. It allows every finite8 cost at every
remaining index and hence overcounts the possible actual configurations.
The TOTAL error is at most `.000880+.005595=.006475`.

Finally (6)--(8) imply

    ||x||²>1.006722852-.006475=1.000247852>1,

contradicting the common farthest normalization. This proves the stated
whole-A19 exclusion.

## Scope and attribution

The arithmetic ingredients are used for the full isolated finite-label
classification, the near-four conjugate/field gate, the D-position
time-square gate, and the fresh far-helper degree argument. The dotted
pair lemma and final summed certificate themselves are geometric at
their stated actual-cap interfaces. The proof does not exclude arbitrary
A19 components with spectators and does not close wholeD19 or compose
an arithmetic dimension-18 theorem.

The summed-foot identity and finite-singleton derivative were developed
by this author in the earlier shared A20 work; root supplied the original
whole-A19 .55 lead and helped check the near-own5 sign argument. The
fresh weak inventory, near-four consumers, actual isolated-dotted pair
lemma and nine-position certificate are written out here for full
outside-author review. No manuscript, formal module, catalogue, accepted
bound or frozen predecessor was changed.
