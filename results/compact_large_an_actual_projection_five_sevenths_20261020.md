# A components of rank at least seventeen: actual projection above 5/7

20 October 2026. Complete paper candidate awaiting independent review.
No scientific execution, optimization, enumeration, arithmeticity,
Lean, TeX or predecessor edits.

**Theorem.** At a vertex farthest from a strict interior point of a
compact Coxeter polytope, every A_r component, r>=17, satisfies

    ||x_A||^2>5/7,

with arbitrary finite spectators. The proof uses actual feasible
points and full component vectors. It is not an edge-only estimate.
No further global dimension ceiling is asserted in this draft.

## 1. Complete inventory input and a stronger forcing certificate

Pin the full actual-wall tools, spectator cost facts, ALL rank-uniform
central3/4 and end1/2 inventories, and actual L/H/F feet in
`compact_large_an_actual_half_uniform_transfer_20261020.md`, SHA256
`7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`.
All paths are under `results/`. Its geometric inventories are
independent of an energy threshold. Every previously threshold-based
discard is explicitly reconsidered below at the stronger target.

Use N=r+1>=18, unit A_r normals n_i=(e_i-e_(i+1))/sqrt2,
C_ij=2min(i,j)(N-max(i,j))/N, and u_i=omega_i/sqrt(C_ii).
Write the strict interior point as E(p-x), h_i=x_A.u_i, and assume
||x_A||^2<=5/7. Distinct actual opposites obey
F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2).

Pairs1/2,(r-1)/r have correlation at least sqrt(8/17); pairs3/4,
(r-3)/(r-2) have correlation at least sqrt(7/10). Their distinct
pair sums respectively exceed53/50 and13/10. Exact seals are

    F(53/100)^2=112360000/164070481<137/200<sqrt(8/17),
    F(13/20)^2=270400/323761<209/250<sqrt(7/10).

The cross-multiplication differences are5655897 and66049; the final
radical comparisons are319073<320000 and43681<43750.

Set L=2u1+2u2+u3+u4 and R=2u_r+2u_(r-1)+u_(r-2)+u_(r-3).
The uniform source proves L^2,R^2<28 and

    L.R=(sum_(i=1..4) a_i sqrt(i/(N-i)))^2, a=(2,2,1,1).

This decreases with N. At N18 its positive square root is

    2/sqrt17+1/sqrt2+1/sqrt5+2/sqrt14
      <243/500+177/250+56/125+67/125=1089/500.

Each comparison follows by positive rational squaring. Hence

    (L+R)^2<56+2*(1089/500)^2<6549/100.

If all four pairs had distinct actual opposites, their weighted
sum would give x_A.(L+R)>171/25. Cauchy is impossible because

    (171/25)^2-(5/7)*(6549/100)=123/17500>0.

After reversing the normal numbering a common actual pair is1/2
or3/4. The complete inventories in the pinned source now apply.

## 2. Threshold-sensitive central profiles and the fourth actual foot

All finite-white, affine-principal, zero-own-contact and spectator-
cost exclusions in that inventory remain unchanged. The target
half-tangent square5/7 corresponds to cosh L=6, so the previously
strong central endpoints must be checked at6.

For adjacent contacts at3,4, the native possibilities are p=3/4/5
at3 and q=3/4 at4. Both3 are affine. If q>=4, the actual direction3
has cosine at least

    sqrt(3/(2N))*((N-3)+sqrt2*(N-4))>6,

already at N18, and the expression increases with N. If p=5,q3,
the actual direction4 has cosine

    (3phi+4)phi^2*sqrt((N-4)/N)>6.

The remaining p4,q3 profile has Q=0 because its retained left A3
cost is3/4. It MUST NOT be discarded merely from its old edge
length. Instead use the NEW actual foot, called C:

    xi=c+n3+(2-sqrt2)n4.

All coefficients are positive. Direct multiplication gives

    xi.c=xi.n3=0, xi.n4=1-sqrt2<0,
    xi^2=4-3sqrt2=-k<0.

Every other actual wall has nonpositive product with this positive
combination. Its positive time coefficient A makes it future, so
xi/sqrtk is an ACTUAL point. It has no spectator spatial component.
With B=3sqrt2/2+2 its full-vector data are

    A_N^2=4+3sqrt2-(17+12sqrt2)/N,
    U^2=A_N^2-k=8-(17+12sqrt2)/N,
    U_raw is supported on1,...,5 with raw sum sqrt2 B.

The isolated single5@3 central profile still has actual cosine
3phi^3*sqrt((N-4)/N)>6. The remaining central single5@1 is foot F,
now retained rather than discarded. Thus a surviving central common
pair supplies foot C or F. A surviving end common pair supplies
one of the same L/H/F types from the complete inventory.

For the displayed endpoint inequalities, at N18 one has
sqrt((N-4)/N)>4/5, phi>8/5, phi^2>13/5 and phi^3>4.
The q4 expression exceeds (2/7)(15+14sqrt2)>6. These bounds and
monotonicity justify the6 cutoffs without a numerical search.

## 3. Four actual feet: exact uniform pressure and density bounds

Every retained foot has xi=A p-U, xi^2=-k<0, U purely in A_r.
It is feasible against EVERY actual wall, being a positive
combination of actual normals as proved in the pinned input or
Section2. Define

    g=(A-sqrtk)/sqrt(U^2),  d_U=U^2/B^2.

Actual farthest comparison gives x_A.(U/sqrt(U^2))>=g. Here B>0
is its raw coordinate sum divided by sqrt2. The following bounds
hold for every N>=18:

| Type | Raw support | g^2 | d_U |
| --- | --- | --- | --- |
| L: cleared A2 dual foot | first3 | >49/100 | >5/9 |
| H: 4@1+3@3 | first4 | >29/50 | >4/9 |
| F: isolated end5,Q0 | first5 | >3/5 | >1/3 |
| C: central4@3+3@4 | first5 | >13/20 | >1/3 |

Here are complete scalar seals; increasing N strengthens all of them.

**L.** With B=b+2d and z0 as in the pinned positive cleared foot,

    8/81<z0<1/3, a_N=2/3-2/N>=5/9,
    A^2=B^2(a_N+z0), k=B^2(z0-3z0^2),
    U^2=B^2(a_N+3z0^2).

The equivalent quadratic for g^2>49/100 at a_N=5/9 is
66603z0^2-19600z0+1445>0. Its discriminant is-805340, proving
the claim for the entire continuous z0 domain.

**H.** At N18,

    A^2=(43+30sqrt2)/36>237/100,
    k=sqrt2-5/4<33/200,
    d_U=(88-6sqrt2)/(9(11+6sqrt2))>4/9.

The first two bounds follow from707/500<sqrt2<283/200. Their
ratio exceeds158/11>6241/441, the threshold for g^2>29/50.
The density comparison is22>15sqrt2.

**F.** At N18,

    A^2=(17phi-19)/36>59/250,
    k=34phi-55<167/12500,
    d_U=(1961-1207phi)/(9(phi+1))>1/3.

The first two bounds use809/500<phi<40451/25000; their ratio is
greater than2950/167>16, giving g^2>3/5. The density comparison
is979>605phi, sealed by phi<16181/10000.

**C.** At N18,

    A^2=(55+42sqrt2)/18>569/90,
    k=3sqrt2-4<1/4,
    d_U=(127-12sqrt2)/(9(17+12sqrt2))>1/3.

These use sqrt2>7/5, sqrt2<17/12 and19>12sqrt2. The ratio is
greater than1138/45>1089/49, proving g^2>13/20.

For H,F,C, k and B are fixed while A_N^2 increases with N; therefore
both g^2 and d_U increase. L has the same property for fixed z0,
with the estimates already uniform over every allowed z0.

## 4. Distinct far walls contradict all four feet

Let W=u_(r-1)+u_r. Its squared norm is less than7/2 because
W^2=2+2sqrt((N-2)/(2(N-1)))<2+sqrt2<7/2. For every raw-support
foot above, its support is before the far end and has sum sqrt2 B,
so, with rho=(U.W)/sqrt(U^2),

    rho=[2/sqrt(N(N-2))+sqrt(2/(N(N-1)))]/sqrt(d_U).

The numerator decreases with N. At N18 it is

    C18=sqrt2/12+1/sqrt153<159/800<1/5.

For the sharper first bound use sqrt2<14143/10000 and
1/sqrt153<1617/20000; their sum is4769/24000<159/800.
Both comparisons are positive rational square tests.

Consequently rho_H<3/10 and rho_F,rho_C<7/20. For L, retain its
full z0 minimum to obtain the stronger necessary estimate

    rho_L^2<(159/800)^2/(5/9+64/2187)
            =55289547/818560000<169/2500.

The last cross-product margin is112772500; hence rho_L<13/50.

If far actual opposites are DISTINCT, the pair calibration in
Section1 requires x_A.W>53/50. In fact every foot contradicts it:

* L: for v=x_A.(U/sqrt(U^2))>=7/10, the two-direction expression
  rho*v+sqrt((5/7-v^2)(7/2-rho^2)) decreases with v and increases
  with rho in this domain. Therefore

      x_A.W<91/500+sqrt(1347217/1750000)
             <91/500+439/500=53/50.

  The last square comparison has numerator margin1830 on the
  common denominator1750000.
* H: separate perpendicular Cauchy gives

      x_A.W<(3/10)sqrt(5/7)+sqrt(47/100)
             <13/50+69/100=19/20.

* F: it gives x_A.W<(7/20)sqrt(5/7)+sqrt(2/5)
             <3/10+16/25=47/50.
* C: it gives x_A.W<(7/20)sqrt(5/7)+sqrt(9/40)
             <3/10+12/25=39/50.

All radical upper bounds follow by squaring positive rationals.
Thus no retained foot is compatible with distinct far opposites.

## 5. Common far walls give one uniform four-type contradiction

If far opposites are COMMON, the reflected complete end inventory
supplies another actual foot of type L,H or F. Every raw support
has at most five coordinates. Since N>=18, the left and right raw
supports are disjoint. Correct reversal of the NORMAL basis gives
raw sums sqrt2 B and -sqrt2 B', so their centered spatial product is

    U.U'=2BB'/N>0.

All four density bounds in Section3 exceed1/3. Hence the normalized
correlation is less than (2/N)/(1/3)<=1/3. Every normalized actual
demand exceeds7/10 because every g^2>49/100. Sum-direction Cauchy
therefore gives

    ||x_A||^2>(7/5)^2/(2+2/3)=147/200>5/7.

This contradicts the contrary projection bound and exhausts the
far-common alternative. The theorem follows.

## Scope and contribution

This is one rank-uniform theorem, not a sequence of tuned isolated
cases. The new central C foot is necessary at the stronger target;
the proof does not discard it using an old half-tangent estimate.
All finite inventories, arbitrary spectators, actual point signs and
full-vector pressure are retained. The argument was derived by
finite57_research, using root's general-N formulas and request for a
single stronger large-A pressure mechanism. New global composition
requires independently reviewed B/D rates and small-defect accounting;
no such additional claim is made here. No machine proof or new
manuscript was produced.
