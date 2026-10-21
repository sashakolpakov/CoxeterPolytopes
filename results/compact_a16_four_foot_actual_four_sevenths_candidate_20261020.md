# Actual A components of rank at least sixteen have energy above 4/7

20 October 2026. Supporting paper-proof candidate by the source-audit
agent, held for full root and outside-author review. No scientific
execution, symbolic search, optimizer, Lean, TeX or predecessor edit.
No new ambient dimension bound is promoted by this note.

**Proposed theorem.** At a vertex farthest from a strict interior point
of a compact hyperbolic Coxeter polytope, every whole A_r component,
`r>=16`, has actual squared projection norm `>4/7`, with arbitrary
finite spectators.

For `r>=17` the accepted `>5/7` theorem is already stronger. The only
new case below is A16, `N=r+1=17`. The argument retains all four actual
feet and uses relaxed constants appropriate to this smaller target.

## 1. Exact inherited inputs and normalization

All source paths below are under `results/`.

* Complete common1/2 and common3/4 finite-white inventories for every
  `N>=9`, zero-own-contact exclusion, arbitrary spectator-cost bounds,
  actual wall pair law and the cleared positive L foot:
  `compact_large_an_actual_half_uniform_transfer_20261020.md`,
  SHA-256 `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`.
* Full L/H/F/C actual-point and raw-vector formulas:
  `compact_large_an_actual_projection_five_sevenths_20261020.md`,
  SHA-256 `7a57b66e1b9b68948c6fe14c493527ff3ac0b963144f904ff512ead723e192b0`.
  Its numerical N18 bounds are **not** imported at N17.
* The flipped window and both-far-pairs strategy occur in the held A19
  source `compact_a19_actual_projection_eleven_fourteenths_20261020.md`,
  SHA-256 `e3fe1ad1696788c4b8e5392c0c45d78aa1ab68b30d54728d09b426ecc2146a26`.
  Every new N17 scalar comparison is supplied below.

Use unit A16 normals `n_i=(e_i-e_(i+1))/sqrt(2)` in sum-zero R17,
dual Gram `C_ij=2 min(i,j)(17-max(i,j))/17`, and unit duals
`u_i=omega_i/sqrt(C_ii)`. Write the strict interior point as `E(p-x)`
and suppose, contrary to the proposed conclusion, that

`e=||x_A||^2<=4/7`.

Put `h_i=x_A.u_i`, with `0<h_i<1`. For distinct actual opposite walls,
the inherited pair law gives
`F(h_i)F(h_j)>=u_i.u_j`, where `F(t)=2t/(1+t^2)`.
Concavity of `log F` therefore gives the following convenient calibrations:

* outer pairs `(1,2),(15,16)` have sum `>1`, since
  `F(1/2)^2=16/25<sqrt(15/32)`;
* inner pairs `(3,4),(13,14)` have sum `>6/5`, since
  `F(3/5)^2=225/289<sqrt(39/56)`.

The radical comparisons follow by squaring: respectively
`256*32<625*15` and `50625*56<83521*39`.

## 2. The flipped window forces an actual common pair

Set

`L=u1+u2+2u3+2u4`, `R=2u13+2u14+u15+u16`.

The six off-diagonal terms of `L^2` are

`2sqrt(15/32), 4sqrt(7/24), sqrt(13)/2,`
`4sqrt(28/45), 4sqrt(13/30), 8sqrt(39/56)`.

Their respective positive rational upper bounds are
`11/8,13/6,181/100,79/25,66/25,167/25`.
Each squares correctly. Including the diagonal10 gives
`L^2=R^2<16699/600<28`.
The reflected cross product has

`sqrt(L.R)=1/4+sqrt(2/15)+2sqrt(3/14)+4/sqrt(13)`
`<1/4+37/100+93/100+111/100=133/50<8/3`.

Consequently `(L+R)^2<56+128/9<72`. If all four selected pairs
had distinct actual opposites, their weighted demand would give
`x_A.(L+R)>34/5`, contradicting Cauchy because

`(34/5)^2-(4/7)*72=892/175>0`.

After reversing the numbering, one of1/2 or3/4 has a common actual
opposite. This relabeling assumes no ambient chamber symmetry.

## 3. Complete retained profiles at the new cosine threshold

The inherited inventories cover N17. Their finite-white, affine,
zero-own-contact and spectator-cost exclusions are independent of the
contrary energy threshold. An edge cosine exceeding `11/3` now
contradicts `e<=4/7`, because `(cosh L-1)/(cosh L+1)>4/7`.

For common3/4, adjacent double contacts have `p=3/4/5` at3 and
`q=3/4` at4. Both3 give an affine principal. If `q>=4`,

`cosh L3>=sqrt(3/34)*(14+13sqrt(2))>161/20>11/3`.

Here `sqrt(3/34)>1/4` and `sqrt(2)>7/5`. If `p5,q3`,

`cosh L4=(3phi+4)phi^2 sqrt(13/17)>2288/125>11/3`,

using `phi>8/5`, `phi^2>13/5`, `sqrt(13/17)>4/5`.
The `p4,q3` profile has Q0 and is retained as C, whether or not
its edge alone could exclude it. Single5@3 has
`cosh L4=3phi^3 sqrt(13/17)>48/5>11/3`; single5@1 is retained as F.
Thus a central common pair supplies C or F.

The common1/2 inventory supplies L/H/F: L consists of its admissible
adjacent contacts, isolated native>=7, or isolated5 with Q>0;
H is4@1+3@3,Q0, and F is isolated5@1,Q0. H/F are deliberately
retained rather than discarded by the older half-energy argument.
Every other end-common alternative is excluded by the same complete
threshold-independent finite/affine/foreign-cost inventory. This uses
arbitrary whole spectators, never an assumption that the full star is finite.

## 4. Honest N17 bounds for all four actual feet

Each foot is a positive combination of actual normals, future timelike,
and feasible against EVERY actual wall. Write

`xi=A p-U`, `xi^2=-k<0`, `g=(A-sqrt(k))/sqrt(U^2)`, `d=U^2/B^2`.

Its spatial U is purely in A16, with raw support in the first at most
five coordinates, raw sum `sqrt(2)B`, and centered form
`U=U_raw-(sqrt(2)B/17)1_17`. Actual farthest comparison gives
`x_A.(U/sqrt(U^2))>=g`. The following deliberately coarse bounds suffice:

| Foot | Raw support | Squared demand | Density |
|---|---|---|---|
| L, positive cleared A2 dual foot | first3 | `g^2>9/20` | `d>1/3` |
| H,4@1+3@3,Q0 | first4 | `g^2>1/2` | `d>1/3` |
| F, isolated end5,Q0 | first5 | `g^2>1/2` | `d>1/3` |
| C,4@3+3@4,Q0 | first5 | `g^2>1/2` | `d>1/3` |

**L.** Its full formulas are
`a=28/51`, `A^2/B^2=a+z`, `k/B^2=z-3z^2`,
`d=a+3z^2`, with `8/81<z<1/3`. Density is immediate.
The inequality `g^2>9/20` is equivalent to
`2523z^2-720z+121a>0`. Since `a>4/9`, it follows from the
quadratic `2523z^2-720z+484/9`, whose discriminant is `-72976/3`.
All spectator clearing and positivity conditions are those of the pinned
L construction; none is changed by this scalar estimate.

**H.** Put `s=sqrt(2)`. The positive foot is
`c+(s-1/2)n1+(s-1)n2+(s/2)n3`. At N17,

`A^2=2(10+7s)/17`, `k=s-5/4`, `B=(3+s)/2`,
`d=(165-12s)/(17(11+6s))`.

Thus `A^2-9k=(845-556s)/68>0` using `s<3/2`, proving
`g^2>1/2`. Density exceeds1/3 because `308>138s`.

**F.** The pinned positive combination is
`c+sum_(j=1)^4 ((j+1)phi-2j)n_j`, with

`A^2=(8phi-9)/17`, `k=34phi-55`, `B=phi/2`,
`d=4(926-570phi)/(17(phi+1))`.

Use `809/500<phi<40451/25000`; their positive sqrt(5) square
tests have margins19 and8401 in the respective integer numerators.
In particular k is positive.
`A^2-9k=(8406-5194phi)/17>0`: substituting the upper phi bound
leaves integer numerator47506. Hence `g^2>1/2`.
The new density comparison is exactly `11095>6857phi`;
the same upper bound leaves numerator2493. This checks the tight
N17 density directly, rather than importing the N18 table.

**C.** Its actual positive combination is
`c+n3+(2-sqrt(2))n4`, with selected products zero against c,n3
and `1-sqrt(2)<0` against n4. At N17,

`A^2=(51+39s)/17`, `k=3s-4`, `B=2+3s/2`,
`d=2(119-12s)/(17(17+12s))`.

Now `A^2-9k=(663-420s)/17>0` using `s<3/2`, so
`g^2>1/2`. Density exceeds1/3 because `425>276s`.
The full-vector construction, including feasibility against unknown
walls and its Q0 conclusion, is unchanged from the pinned C foot.

## 5. Either far-common pair contradicts the first actual foot

If15/16 is common, the reflected inventory supplies L/H/F; if13/14
is common, it supplies C/F. Left and right raw supports are disjoint
at N17. Correct NORMAL reversal gives opposite raw sums, hence

`U.U'=2BB'/17`, `corr(U,U')=(2/17)/sqrt(dd')<6/17`.

Every demand has square `>9/20`. Sum-direction Cauchy therefore gives

`e>(4*9/20)/(2+12/17)=153/230>4/7`.

The last difference is `151/1610>0`. This one estimate covers every
L/H/F/C pairing, retaining the positive centered correlation.

## 6. Both far pairs distinct also contradict the first foot

Their pair laws give `x_A.R>17/5`. For ANY retained left raw support,
the common raw-sum formula gives

`U.R/B=K=2sqrt(8/221)+2sqrt(3/119)+2/sqrt(255)+1/sqrt(136)`
`<2/5+1/3+13/100+9/100=143/150<1`.

These four bounds are positive rational square tests. Since `d>1/3`,
`rho_R=(U.R)/sqrt(U^2)<sqrt(3)`. Let
`v=x_A.(U/sqrt(U^2))`; then `v^2>9/20` and `v<=sqrt(e)`.
Orthogonal decomposition and Cauchy give the simpler bound

`x_A.R <= rho_R*v+sqrt((e-v^2)(R^2-rho_R^2))`
`<sqrt(12/7)+sqrt(17/5)<4/3+15/8=77/24<17/5`.

The two radical seals square to `108<112` and `1088<1125`;
the last comparison is `385<408`. No monotonicity derivative
argument or numerical optimization is needed. Thus neither far-common
nor both-far-distinct alternatives are possible, proving A16's proposed
`>4/7` bound. The accepted rank>=17 theorem completes the uniform statement.

## Conditional research consequence and status

In the conservative dimension23 frontier, A16 has spectator rank7.
After an I debit the remaining six ordinary ranks give
`4/7+6/14=1`; with one rank-two7 defect they give
`4/7+5/14+9/100=713/700>1`. These are strict contradictions if
this candidate and the existing spectator accounting are accepted.
The held A19/A20 and B22/D22/whole-D23 sources require their separate
root acceptance before any global composition or bound promotion.

This note adapts the existing four-foot strategy; it credits the
inherited actual-point and inventory arguments to their pinned sources.
Its new content is the N17 window, relaxed pair calibration, explicit
F/C densities and the uniform four-type scalar closure. No published
input, actual geometric premise, or proof-execution status is enlarged.
