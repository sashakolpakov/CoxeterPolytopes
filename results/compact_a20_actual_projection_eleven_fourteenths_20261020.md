# A20: actual farthest projection above 11/14

20 October 2026. Complete paper candidate awaiting independent review.
No scientific execution, optimization, enumeration, arithmeticity, Lean,
TeX or predecessor edits.

**Theorem.** At a vertex farthest from a strict interior point of a
compact hyperbolic Coxeter polytope, an A20 component satisfies

    ||x_A||^2>11/14,

with arbitrary finite spectators. This controls the ACTUAL component
of the same farthest vector, not an edge-only covering relaxation.
No new ambient dimension ceiling is claimed here.

## 1. Fixed inputs and the contrary projection bound

All paths below are under `results/`. Import the complete arbitrary-
spectator finite inventories, actual reciprocal-shadow/positive-point
interfaces and the positive cleared L foot from
`compact_large_an_actual_half_uniform_transfer_20261020.md`, SHA256
`7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`.
Import the four actual feet L/H/F/C, their full-vector formulas and
their exact pressure/density seals from Sections2--3 of
`compact_large_an_actual_projection_five_sevenths_20261020.md`, SHA256
`7a57b66e1b9b68948c6fe14c493527ff3ac0b963144f904ff512ead723e192b0`.
Its outside audit is
`review_compact_large_an_five_sevenths_batch_guard_20261020.md`, SHA256
`0ee81c7e00a0681c811bb311746e797a98c3e4328f049233ae78006d3d5f986c`.
The complete geometric inventories are threshold-independent. We do
NOT reuse a profile discard at the old5/7 threshold: every discard
depending on an edge length is reconsidered in Section3 below.

Put N=21, use unit A20 normals n_i=(e_i-e_(i+1))/sqrt2, and write

    C_ij=2min(i,j)(21-max(i,j))/21,
    u_i=omega_i/sqrt(C_ii), h_i=x_A.u_i.

Write the strict interior point as E(p-x), and suppose
e=||x_A||^2<=11/14. Actual edges give h_i>=tanh(L_i/2)>0.
For DISTINCT actual opposite walls,

    F(h_i)F(h_j)>=u_i.u_j, F(t)=2t/(1+t^2).

Since 0<h_i<1 and log F is strictly concave, the outer pairs1/2,
19/20 have sum>53/50, and the inner pairs3/4,17/18 have sum>13/10.
These are the exact pair seals in the pinned5/7 source, whose smaller
N18 correlations already suffice. Explicitly,

    F(53/100)^2=112360000/164070481<137/200<sqrt(8/17),
    F(13/20)^2=270400/323761<209/250<sqrt(7/10).

## 2. Flipped weights force an actual common pair

Use the heavier INNER weights

    L=u1+u2+2u3+2u4,
    R=2u17+2u18+u19+u20.

The exact six off-diagonal terms of L^2, and positive rational
upper bounds for them, are

    2sqrt(19/40)<69/50,   4sqrt(3/10)<11/5,
    4sqrt(17/80)<37/20,  4sqrt(12/19)<16/5,
    4sqrt(17/38)<67/25,  8sqrt(17/24)<337/50.

Each follows by squaring positive rationals. Including the diagonal
10 gives L^2,R^2<561/20. Also

    sqrt(L.R)=1/sqrt20+sqrt(2/19)+2/sqrt6+4/sqrt17
      <9/40+13/40+41/50+49/50=47/20.

Consequently

    (L+R)^2<2*(561/20)+2*(47/20)^2=13429/200<68.

If all four selected pairs had distinct actual opposites, then
x_A.(L+R)>2*(53/50)+4*(13/10)=183/25. This contradicts Cauchy:

    (183/25)^2-(11/14)*68=673/4375>0.

After reversing the normal numbering, a common actual pair is1/2
or3/4. Relabeling assumes no ambient chamber symmetry.

## 3. Complete retained-foot list at the stronger target

The pinned uniform inventory includes all zero-own-contact profiles,
all foreign Schur-cost branches, and the finite/affine exclusions;
none of those arguments uses the bound e<=1/2 or e<=5/7. At the new
target an edge cosine>25/3 suffices, since its half-tangent square
then exceeds11/14.

For common3/4, the only adjacent double contacts are native p=3/4/5
at3 and q=3/4 at4. Both3 are affine. If q>=4, actual direction3
has cosine at least

    (18+17sqrt2)/sqrt14>836/75>25/3.

Here sqrt2>7/5 and sqrt14<15/4. For p5,q3, actual direction4 has
cosine

    (3phi+4)phi^2 sqrt(17/21)>4576/225>25/3,

using phi>8/5, phi^2>13/5 and sqrt(17/21)>8/9. The remaining p4,q3
profile has Q=0 and supplies actual foot C; it is retained.
For the single5@3 profile, the cosine is

    3phi^3 sqrt(17/21)>32/3>25/3,

using phi^3>4. The single5@1 profile supplies foot F. Thus every
surviving common3/4 supplies C or F. Common1/2 supplies L/H/F:
the formerly excluded H/F profiles are retained here, while all
other common1/2 cases are exactly the threshold-independent L rows
or the finite/affine/foreign-cost impossibilities of the pinned
inventory. No finite full-star assumption is imposed on spectators.

We recall the FULL point data needed below. Each retained foot is

    xi=A p-U, xi^2=-k<0, g=(A-sqrtk)/sqrt(U^2),
    d_U=U^2/B^2, B>0.

It is a future timelike positive combination of actual normals,
feasible against EVERY actual wall. After spectator clearing, U is
purely in A20. Actual farthest comparison therefore gives

    x_A.(U/sqrt(U^2))>=g.

The pinned exact formulas prove the following bounds for ALL N>=18,
in particular N21; those scalar bounds do not assume e<=5/7:

| Foot | Raw support | g^2 | d_U |
| --- | --- | --- | --- |
| L: positive cleared A2 dual foot | first3 | >49/100 | >5/9 |
| H: 4@1+3@3,Q0 | first4 | >29/50 | >4/9 |
| F: isolated end5,Q0 | first5 | >3/5 | >1/3 |
| C: 4@3+3@4,Q0 | first5 | >13/20 | >1/3 |

Every raw vector has sum sqrt2 B and centered form
U=U_raw-(sqrt2 B/21)1_21. L is the cleared A2 foot with
8/81<z0<1/3 and U^2/B^2=2/3-2/21+3z0^2. H/F/C are the positive
actual combinations displayed in the pinned source. Thus every
non-L foot satisfies g^2>29/50 and d_U>1/3. We retain all four
possibilities rather than treating their edge-only floors as enough.

## 4. Either common far pair gives a second actual foot

If19/20 is common, the reflected full end inventory supplies L/H/F.
If17/18 is common, the reflected central inventory, including the
stronger threshold checks above, supplies C/F. The raw supports of
the left and right feet are disjoint (first and last at most five
coordinates). Correct NORMAL reversal has raw right sum -sqrt2 B',
so their CENTERED product is POSITIVE:

    U.U'=2BB'/21,
    rho=(2/21)/sqrt(d_U d_U').

All three possible type groupings now contradict e<=11/14 by
sum-direction Cauchy:

* L/L: rho<6/35 and both demands>7/10, hence

      e>(7/5)^2/(2+12/35)=343/410>11/14.

* L/non-L: rho^2<12/245<1/16, so rho<1/4. The demands exceed
  7/10 and19/25 respectively, since (19/25)^2<29/50. Hence

      e>(73/50)^2/(5/2)=5329/6250>11/14.

* Non-L/non-L: rho<2/7 and both squared demands>29/50, hence

      e>(4*29/50)/(2+4/7)=203/225>11/14.

This covers BOTH far pairs and every retained L/H/F/C combination.
The positive reversal mean has not been replaced by a negative
correlation estimate.

## 5. Both far pairs distinct contradict the first actual foot

It remains that17/18 and19/20 are BOTH distinct. Their pair laws give

    x_A.R>2*(13/10)+53/50=183/50.

Because every left raw support lies before17 and has sum sqrt2 B,
the correlation with the SAME weighted vector R is

    U.R/B=K=2sqrt(8/357)+2/sqrt63+2/sqrt399+1/sqrt210
      <3/10+13/50+11/100+7/100=37/50<3/4.

All four comparisons are positive rational square tests. Put
rho_R=(U.R)/sqrt(U^2)=K/sqrt(d_U); this is not a normalized unit-
unit correlation, since R is not unit. We already have R^2<561/20.

For L, rho_R^2<(37/50)^2*(9/5)=12321/12500<1 and the actual
unit-foot projection v exceeds7/10. Perpendicular Cauchy gives

    x_A.R<=rho_R*v+sqrt((11/14-v^2)(561/20-rho_R^2)).

On v>=7/10, 0<=rho_R<=1, this expression decreases with v and
increases with rho_R: both derivative signs reduce to
v^2*(561/20)>(11/14)*rho_R^2, immediate from these bounds.
Consequently

    x_A.R<7/10+sqrt((207/700)*(541/20))
          =7/10+sqrt(111987/14000)
          <7/10+29/10=18/5<183/50.

For the radical comparison, 111987/14000<8<(29/10)^2.

For a non-L foot, d_U>1/3 and K<3/4 give
rho_R<(3/4)sqrt3<4/3, while v^2>29/50. Separate perpendicular
Cauchy therefore gives

    x_A.R<(4/3)sqrt(11/14)
             +sqrt((11/14-29/50)*(561/20))
          <6/5+121/50=181/50<183/50.

Both rational seals are explicit:

    (6/5)^2-(4/3)^2*(11/14)=68/1575>0,
    (121/50)^2-(36/175)*(561/20)=1507/17500>0.

Thus both far pairs cannot be distinct either. Section4 exhausted
the common alternative, proving ||x_A||^2>11/14.

## Scope and contribution

The flipped four-pair weights and coupling to BOTH far pairs are
necessary new ingredients; this does not substitute a stronger
target into an old weaker forcing inequality. It uses actual full
feet, keeps arbitrary spectators, and explicitly retains H/F/C when
their old edge discard is unavailable. The theorem is geometric
and componentwise. Any global dimension composition remains a
separate reviewed obligation.

Finite57_research derived the exact rational window/correlation
seals and the complete four-foot closure, following root's proposed
outer1:inner2 weights, both-far-pairs strategy, and independent
general-N foot identities. The predecessor proofs retain their
documented contributions. No accepted source was changed and no
machine or formal proof was run for this note.
