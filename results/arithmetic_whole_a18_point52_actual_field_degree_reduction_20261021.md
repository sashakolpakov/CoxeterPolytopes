# Whole A18 at fresh rank18: actual field degree at most three

21 October 2026. Supporting written-proof candidate, pending full root and
outside-author review. Manual mathematics only: no census, scientific job,
optimizer, Lean, TeX or primary-paper edit. The source-audit author derived
the fresh .52 inventory and degree reduction. Root independently checked
the isolated-cap estimates and proposed the quartic Minkowski argument;
the rational discriminant seal below was checked by both. These are shared
inputs, not an assertion of independent review of this source.

**Theorem.** Suppose an actual farthest vertex of a compact arithmetic
first-type hyperbolic18-polytope has WHOLE finite diagram A18. Its ACTUAL
admissible ground field K has degree at most3 over Q. If its degree is3,
every .52-weak index in5..14 has an actual ISOLATED DOTTED opposite, and
at least one such index exists. For each of these opposites, beta=4gamma^2
generates K, is an algebraic integer, and satisfies

    4<beta<40/3,       0<beta_sigma<19/35 (sigma nonphysical).

Its irreducible polynomial is necessarily X^3-T X^2+S X-M with
T in5..14, S in1..14 and M in1..3. These finite coefficient bounds are
necessary conditions, not a claim that any listed triple is realizable.
The remaining quadratic and cubic profiles are NOT excluded here. No
whole-A18 exclusion or dimension17 ceiling follows.

## 1. Exact inputs and external boundary

All local input paths below are under results/. Only the specified
formulas and classification mechanisms are used, not their higher-rank
conclusions or their old weak inventories at a different loss.

| Input | SHA256 | Use here |
| --- | --- | --- |
| arithmetic_large_a_central_cap_packing_20261021.md | `8599baeda0f1cfae41bf3e4f27a89d7d999a49a3284caa8bdd559ea1ba945ee0` | N19 inverse costs and normalization |
| arithmetic_full_an_isolated_finite_continuation_packing_20261020.md | `9ab34f073ea95bed2040bafe6b8573331721f9b1e75d65172f814637da046810` | complete isolated cyclotomic argument and raw normalization |
| arithmetic_uniform_classical_weak_helper_obstruction_20261021.md | `535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905` | native3 normalization and cyclic integrality |
| compact_whole_a21_farthest_vertex_exclusion_candidate_20261020.md | `3abe11bb70ce677069279462a91111ed6b8ab8889e3ec9e656ac71c1094200fb` | finite-white classification proof and role formulas, rechecked at18 below |
| compact_farthest_vertex_mixed_replacement_inequalities_20261020.md | `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84` | actual feasible-point/farthest comparison |
| whole_a20_summed_isolated_feet_technical_input_20261021.md | `8d93c23c9bf268813334dd12fef09f37f523eb22208971847cf8b04e1d659ab4` | isolated dotted actual-foot formula only |

External hypotheses are compact acute simplicity, finite spherical
Coxeter classification, the native finite-angle gap, and the standard
arithmetic first-type model: field-valued normals, cyclic integrality,
and positive definiteness at every nonphysical real embedding. K is
totally real and is the actual field, not a chosen field from a list.

One standard number-theoretic theorem is used in Section7: for a totally
real degree-d field, Minkowski gives sqrt(disc K)>=d^d/d!. See Romyar
Sharifi, *Algebraic Number Theory*, Chapter4, Corollary4.3.6, with its
derivation in Proposition4.3.4 and Definition4.3.5:
[primary lecture notes](https://www.math.ucla.edu/~sharifi/notes/algnum-ch04.html).
We also use the discriminant-index identity
disc(Z[alpha])=disc(K)[O_K:Z[alpha]]^2 for an integral primitive alpha;
it follows by change of integral basis. No field-list or quartic-field
classification is an input, and no D19 field reduction is imported.

## 2. A fresh .52-weak index must exist

Use the sum-zero subspace of R19 and UNIT old normals
n_i=(e_i-e_(i+1))/sqrt2. At the actual farthest vertex p write the
strict interior point in the usual form z=E(p-x), with ||x||<1 and
x.n_i>0. The spatial duals and inverse entries are

    omega_i=sqrt2(P_i-i P19/19),
    D_i=omega_i^2=2i(19-i)/19,
    C_ij=omega_i.omega_j=2min(i,j)(19-max(i,j))/19.

Set h_i=x.omega_i and h_0=h_19=0. A central index5<=i<=14 is weak
when h_i<=sqrt(D_i)-13/25. Its D_i is at least140/19 and less than19/2.

Consider the EXACT rational heights

| i | 5 | 6 | 7 | 8 | 9 | 10 |
| --- | --- | --- | --- | --- | --- | --- |
| H_i | 2.194 | 2.345 | 2.453 | 2.523 | 2.557 | 2.557 |

Reflect by H_(19-i)=H_i and interpolate linearly from0 to5 and14 to19.
The residual numerators of D_i-(H_i+.52)^2, in denominator19,000,000,
for i5..9 are respectively49876,43725,64149,62869,109349, all positive.
The successive left slopes are .4388,.151,.108,.070,.034,0, so this
height function is concave. Its vector w has nonnegative old-dual
coefficients supported on5..14, and

    ||w||^2=sum beta_i H_i
      =2.194^2/5+.151^2+.108^2+.070^2+.034^2
      =1.0032482>1.

If all ten central indices were strong, x.w>||w||^2 and completing
the square would give ||x||^2>=2x.w-||w||^2>1. Hence a weak index exists.
This does not assume any old .61 or .55 prefix floor.

For a weak index, the exact edge ratio R=sigma/sqrt(delta) cannot exceed
the sufficient strong threshold

    T(q)=[1+(1-.52/q)^2]/[q(1-(1-.52/q)^2)]
        =25/13-(q-13/25)/[q(2q-13/25)], q=sqrt(D_i).

T increases on this q-range, and q<31/10 gives
T(q)<101665/57226<89/50. The last cleared difference is9864.
Only the implication R>T(q) => strong is needed.

## 3. Complete white inventory at the changed loss

Deleting i leaves A_(i-1)+A_(18-i), both of rank at least4. The actual
opposite c together with these white roots is finite. A finite diagram
is a tree, so c contacts at most one node of each white path. Its
COMPLETE possibilities are: isolated; one white3 at an end, D-position,
or admissible E-position; one white4 at an end; or two white3 contacts
at path ends. The single-path augmented rank is at least5, excluding
H/F exceptions. A two-path join has connected rank18>8 and is A18:
an internal4 cannot be B and a branch with both long arms cannot be D.
There is no foreign component at this WHOLE rank18 vertex.

Own zero contact is impossible. One endpoint3, next-to-end3 or
endpoint4 would give positive definite A19,D19 orB19 in signature(18,1).
More internal3 gives an affine E principal star; internal4 gives the
affine squared-cosine star(1/2,1/4,1/4). Two white3 contacts give an
all3 affine cycle. With no white contact, c would lie in the timelike
orthogonal complement of the positive old18-space. Own parallel contact
is also impossible: the nonzero future null sum of the two actual unit
normals is feasible against every wall, contradicting compactness.

For one white3, let its side have m-1 nodes, m=i or19-i>=5, numbered
j from the far end. The exact white residual is

    delta=1-j(m-j)/(2m),       sigma=gamma+j/(2m).

Own3 gives an all3 affine cycle; own6 gives an affine6--3 path using
an uncontacted old neighbor of i. The following table covers every
remaining NONISOLATED weak possibility; permissive quadratic entries
need not all be geometrically possible.

| White role | Fresh .52 conclusion |
| --- | --- |
| far endpoint3, j1 | Section4: actual K is Q(sqrt5) or Q(sqrt2) |
| near endpoint3, j=m-1 | dotted gamma>=1 is strong since R>=(1+.4)/sqrt(.6)=7/sqrt15>89/50; finite own>=7 is excluded by Section5; only quadratic own4/5 can remain |
| far D-position3, j2 | dotted is strong: R>=(m+1)/sqrt(2m)>=6/sqrt10>89/50; finite own>=7 is excluded by Section5; own4/5 is at most quadratic |
| near D-position3, j=m-2 | same dotted bound and Section5 finite exclusion; own4/5 is at most quadratic |
| E-position3 | delta<=1/4, so finite own>=7 or dotted gives R>=2gamma>1.8>89/50; remaining own4/5 is quadratic |
| endpoint4 | far own3 is affine F4; near own3 has R>2; own>=4 has R>=(m+1)/sqrt(2m)>=6/sqrt10>89/50 |
| two white endpoint3 | own3 is affine; for own>=4, R>=q/sqrt2+1/q>2>89/50 |

In the far endpoint3 row, own4 is the affine star at i, while own3/6
were just excluded. Thus t=2gamma is at least phi>202/125 before
Section4's upper bound. The own4 near endpoint3 row is deliberately
retained: lowering the loss does not entitle us to discard it.

If there is a native3 helper, rescale c OVER K to square2, as the
old A roots, so its helper product is-1 and marked product is-t,
t=2gamma. Cyclic integrality gives t in O_K. The full old18+c Gram
has entries in Q(t) and is physically Lorentzian, with nonzero time
coefficient. If K properly extended Q(t), an extra real embedding fixing
Q(t) would repeat this Lorentzian Gram in a positive definite conjugate
space. Therefore K=Q(t). In particular own4/5 remaining in the table
locks the ACTUAL field to its quadratic field, with no square-root
extension silently made.

## 4. Far endpoint3: a fresh all-degree norm/trace seal

The exact identity, with delta=(m+1)/(2m), is

    R^2-(4gamma-2)=(gamma-1-1/(2m))^2/delta>=0.

Weakness therefore forces gamma<13/10, since at gamma=13/10 the lower
bound16/5 already exceeds(89/50)^2. Thus202/125<t<13/5.

At every nonphysical embedding the square-two Gram Schur complement
gives D_i t_sigma^2+2 C_ie t_sigma+D_e<4. The N19 inverse-cost
calculation of the pinned input, independent of its old weak threshold,
gives EVERY other conjugate in

    -a<t_sigma<b,       a=193/250, b=12/25.

For completeness the least negative/positive costs on the left half are
Q_-=2it^2+4t+2-2(it+1)^2/19 and
Q_+=2it^2+2-2(it-1)^2/19. They increase with i there; at i5 their
values minus4 at -a and b are643/59375 and232/2375. Every embedding
of the minimal field extends to K, so these bounds omit no conjugate.

If degree d>=5, |Norm t|<(13/5)a^4<1, impossible. In degree4 the
absolute norm is<2; a positive small conjugate would instead give
|Norm t|<(13/5)b a^2<1. Hence the three small conjugates are-a1,-a2,-a3,
all a_j in(0,a), with product1/t. Integral trace is0: it is>-1, and
AMGM bounds it by t-3t^(-1/3)<1 for t<13/5. Thus sum a_j=t, so
t^4>=27. The second integral coefficient is

    S=-(sum a_j^2+sum_(j<k)a_j a_k),
    -4<-6a^2<S<=-(2/3)t^2<=-2sqrt3<-3,

impossible. The upper trace bound follows, for example, from
t(t-1)^3<27 for1<t<13/5.

In degree3 the absolute norm is<2, and any positive small conjugate
would give absolute norm<(13/5)b a<1. Write the two conjugates-a1,-a2,
with product1/t. Integral trace is1 or2. Trace2 would give
0<a1+a2=t-2<.6 but a1 a2=1/t>5/13, impossible. Trace1 gives
t<1+2a<51/20 and t(t-1)^2>=4, hence t>9/4. The second coefficient
S=t-t^2+1/t is strictly between-4 and-2, so S=-3. Its polynomial
X^3-X^2-3X-1=(X+1)(X^2-2X-1) is reducible, again impossible.

In degree2 the trace is1,2 or3. Trace3 forces63/25<t<13/5 and
1<26/25<t(3-t)<756/625<2, impossible. Trace2 has the nonzero
integer norm t(2-t) in(-2,1), hence norm-1 and t=1+sqrt2. Trace1
has t<1+a and norm t(1-t) in(-2,0), hence norm-1 and t=phi.
In degree1 the only physical integer is2, the forbidden parallel pair.
Consequently every weak far helper has actual quadratic field.

## 5. Finite near and D-position labels: no missed higher-degree helper

Use the same t=2cos(pi/l), square-two normalization and K=Q(t).
If the white helper is j, every nonphysical conjugate obeys

    Q_(i,j)(t_sigma)=D_i t_sigma^2+2 C_ij t_sigma+D_j<4.

Reflection allows near positions j=i-1 or i-2 and far D-position j2,
always5<=i<=14. The following boundary tests specify the entire interval
of this convex quadratic containing its minimum:

| j | Left boundary | Right boundary | Interior test |
| --- | --- | --- | --- |
| i-1 | Q(-3/2)=(148+15i-i^2)/38>=81/19>4 | Q(0)=D_(i-1)>=120/19>4 | Q(-1)=36/19<4 |
| i-2 | Q(-3/2)=(288+11i-i^2)/38>=123/19>4 | Q(0)=D_(i-2)>=96/19>4 | Q(-1)=68/19<4 |
| 2 | Q(-3/2)=[(19-i)((9/2)i-12)+68]/19>=215/19>4 | Q(2/5)>=532/95>4 | Q(0)=68/19<4 |

Thus allowed conjugates lie in(-3/2,0) in the near cases and in
(-3/2,2/5) in the far D case. The endpoint bounds in i follow from
concavity or the displayed minimum inverse entries, not an enumeration.

For l7 a conjugate is2cos(3pi/7)>2/5, excluded in every case. The
small positive root of X^3-X^2-2X+1 lies between2/5 and1/2, which
is an exact alternative verification. For even l>=8, the coprime
exponent l-1 modulo2l gives conjugate-2cos(pi/l)<-3/2. For odd l>=9,
the coprime exponent l-2 gives-2cos(2pi/l)<=-2cos(2pi/9)<-3/2;
the final inequality follows from the largest root of X^3-3X+1
being>3/2. These exhaust ALL finite l>=7, including9. Own3/6 were
affine. Therefore no omitted high-degree finite near/D helper remains.

## 6. Isolated caps and their physical bound

For an isolated FINITE cap put beta=4cos^2(pi/l). Over the ground field
it has a raw normal c^2=2beta and marked old product-beta. The full
Lorentzian star locks K=Q(beta); every other conjugate satisfies
0<beta_sigma<4/D_i<=19/35<1. Since the raw Cartan inverse diagonal
D_i/2>=70/19>2, the complete elementary cyclotomic proof in9ab leaves
only finite labels5,8,12 after the affine3/4/6 exclusions. Each has a
quadratic field. No assertion that each is .52-strong is made or needed.

For an isolated DOTTED cap c=A p-gamma omega_i, gamma>1, put
beta=4gamma^2 and u=gamma^(-2)=4/beta in(0,1). The raw normalization
c^2=2beta, c.old_i=-beta is over K; the square cycle gives beta integral.
Again K=Q(beta) by the full Lorentzian-star argument, and every other
conjugate lies in(0,delta), delta=19/35.

The positive combination of ACTUAL normals c+gamma^(-1)n_i is a
feasible future point of norm-(1-u). Farthest comparison gives

    x.(omega_i-u n_i)>=sqrt(D_i-u)-sqrt(u(1-u)).

Since x.n_i>0 and
1/[sqrt(D_i)+sqrt(D_i-1)]<191/1000, its height loss is less than

    .191u+sqrt(u(1-u)).

Indeed sqrt(D_i)>2.714 and sqrt(D_i-1)>2.523, and
5237*191=1000267>1000000. If u<=3/10, the last expression is at most
.0573+sqrt(.21)<.5163<.52. Hence a WEAK isolated dotted cap necessarily
has u>3/10 and therefore

    4<beta<40/3,       0<beta_sigma<delta=19/35.       (1)

This physical upper bound comes from one actual feasible point and the
same original x; no formal wall or hypothetical prefix floor is used.

## 7. Exclude degree at least four without a coefficient census

Write d=[K:Q]=[Q(beta):Q] for an isolated dotted WEAK cap.
If d>=6, its positive nonzero integral norm is less than
(40/3)delta^5=99043960/157565625<1, impossible.

If d5, its norm is<(40/3)delta^4=5212840/4501875<2, so it is1.
Each of the four small conjugates a_j consequently exceeds
3/(40delta^3)=128625/274360>16/35. Thus |1-2a_j|<3/35 and

    0<|Norm(1-2beta)|<(77/3)(3/35)^4<1,

impossible for a nonzero algebraic integer. It is nonzero because the
physical beta>4, and beta is irreducible.

If d4, the norm is<(40/3)delta^3=274360/128625<3, hence1 or2.
Norm2 forces each small conjugate>3/(20delta^2)=3675/7220>16/35;
then the same argument with exponent3 gives
0<|Norm(1-2beta)|<(77/3)(3/35)^3<1. Thus only a quartic UNIT could
remain at this stage.

Write its physical root as b and its other three as a1<=a2<=a3.
The integral trace is<=14, since40/3+3delta<15. Norm1 and AMGM give
sum a_j>=3/b^(1/3)>5/4 because40/3<(12/5)^3. Thus b<51/4, and
every a_j exceeds

    ell=1/[(51/4)delta^2]=4900/18411>1/4.

Put L=delta-ell=178309/644385<277/1000. The ordered small roots obey
prod_(j<k)(a_k-a_j)^2<=L^6/16: their outer gap is at mostL and the
product of the two inner gaps is at mostL^2/4. Also b-a_j<51/4-ell<25/2.
The discriminant of the primitive integral generator beta therefore is

    disc(Z[beta]) < [(25/2)(277/1000)]^6/16
                  =(277/80)^6/16 <12^3/16=108,

because277^2=76729<76800=12*80^2. But Minkowski's stated external
theorem gives disc(K)>=(4^4/4!)^2=1024/9>108, and the discriminant-index
identity makes disc(Z[beta])>=disc(K). Contradiction. This excludes the
quartic unit as well as all larger degrees, without enumerating fields.

## 8. Fresh-rank conclusion and exact remaining interface

Section2 forces an actual .52-weak index. Sections3--5 make every
nonisolated weak cap quadratic in its ACTUAL field; isolated finite
weak caps are also quadratic. Section7 bounds every isolated dotted
weak cap by degree3. The claimed degree<=3 follows for the actual
whole-A18 vertex, over ALL admissible first-type fields.

If K is cubic, every weak central cap is therefore isolated dotted.
For any one of them, write its roots b,a1,a2 as in(1). The positive
integral norm M is<(40/3)delta^2=14440/3675<4, so M in1..3. Its
trace satisfies4<T<40/3+2delta<15, hence T in5..14. Its second
coefficient satisfies0<S<(80/3)delta+delta^2=54283/3675<15, hence
S in1..14. These prove the finite cubic box stated at the outset.

Nothing here counts formal cubic polynomials as actual caps, declares
the cubic box empty, or assumes compatible weak caps can be chosen
independently. Completing a whole-A18 exclusion still requires actual
same-field/time/pair restrictions and one legal same-x energy argument
for the remaining quadratic/cubic configurations. No new enumeration
is proposed or executed by this source, and all frozen inputs remain
unchanged.
