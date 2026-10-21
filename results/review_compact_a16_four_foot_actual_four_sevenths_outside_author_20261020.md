# Outside-author review: actual A_r projection above 4/7 for r >= 16

20 October 2026. Independent mathematical source review only. No scientific
worker, symbolic computation, optimizer, Lean or TeX was executed; no proof
source was edited. **PASS at the written-proof level.** This accepts the
component theorem below, not a new ambient dimension ceiling.

## Frozen scope and source bindings

Reviewed `results/compact_a16_four_foot_actual_four_sevenths_candidate_20261020.md`,
SHA-256 `a44e1d7f154ff1fab351737e229293778fc4d793b34b247448c2c2ac23a76f6c`
(229 lines). The theorem is actual squared projection `||x_A||^2 > 4/7`
for every A_r component, r >= 16, at an actual farthest vertex of a compact
hyperbolic Coxeter polytope, with arbitrary finite spectators. Its only new
rank is A16: the accepted r >= 17 theorem is stronger.

I read the full candidate and independently checked the mathematical
composition against these unchanged inherited sources, including their
actual-point rather than abstract-Gram hypotheses:

* `compact_large_an_actual_half_uniform_transfer_20261020.md`,
  `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7`;
* `compact_large_an_actual_projection_five_sevenths_20261020.md`,
  `7a57b66e1b9b68948c6fe14c493527ff3ac0b963144f904ff512ead723e192b0`;
* corrected `compact_universal_a7_two_end_actual_projection_20261020.md`,
  `5ae61f571761f517222bef38f35468462ec8fb9ad9d879850c419a6c316828a5`;
* the positive acute projection, reciprocal shadow and orthogonal-component
  opposite obstruction sources pinned in that A7 source, lines25–34.

All paths in this list are under `results/`. The held A19 source supplies
strategy attribution only: no A19 theorem or scalar conclusion is needed
for this proof. Its separate acceptance, and the held D23/A20 composition,
are outside this review.

## Actual geometry and inventory coverage

The normalization at candidate lines35–44 is consistent: `N=17`, unit old
normals, inverse Gram `C_ij=2 min(i,j)(17-max(i,j))/17`, and
`z=E(p-x)`. For any feasible future point `xi/sqrt(k)`, with
`xi=A p-U` and `xi^2=-k`, farthest comparison gives
`x_A.U >= A-sqrt(k)` whenever U is purely in A. No required energy
conclusion is supplied as an input.

The inherited finite-white classification is applied only to replacements
of already actual opposite walls. Positive acute projection supplies the
actual endpoint when the full replacement is positive definite. A foreign
component has at most one contact; its cost q is at least1/4, with no cost
in(1/4,1/3). Reciprocal shadow plus the orthogonal-component obstruction
gives the necessary strict bound `q<A^2`. Thus spectator-cost exclusions
do not assume spectator rank bounds or a finite full star.

The common1/2 and common3/4 inventories in the half-transfer source,
lines93–178, remain exhaustive at N17. For the central double case, the
two retained old sides have enough rank that every nonadjacent double
contact is either already nonpositive in an actual white replacement or
has both bonds3 and a principal affine cycle. At3/4, the native labels
are exactly p=3/4/5 and q=3/4. The remaining single-contact and end-common
cases use finite/affine/cost exclusions independent of the contrary energy
target. In particular, end H and F are retained rather than silently
discarded using the old half bound.

Candidate lines88–104 correctly change the cosine cutoff to
`(1+4/7)/(1-4/7)=11/3`. All three displayed central edge lower bounds
exceed11/3. The residual p4,q3 profile has Q=0: its retained left A3
cost is3/4, leaving Q<1/4. Its actual foot is C. Single5@3 is excluded
at the new cutoff and single5@1,Q0 is retained as F. Thus every initial
and reflected common pair supplies exactly one of the advertised feet.

## Independent scalar and actual-foot checks

The distinct-opposite law is used only for distinct walls. Concavity of
`log(2t/(1+t^2))` on(0,1) gives the pair sums in lines47–53 from the
stated exact correlations. I checked all six terms in the flipped
window's squared norm, its rational upper bounds, and the reflected cross
product. Hence `L^2=R^2<28`, `(L+R)^2<72`, and the all-distinct demand
`x_A.(L+R)>34/5` contradicts energy at most4/7 with margin892/175.
Reversal only relabels the A diagram and assumes no chamber symmetry.

All four feet at lines114–176 are actual feasible future points. For H,
the selected products are zero at c,n1,n3 and negative at n2. For F,
they are zero at c,n1,n2,n3 and `3phi-5<0` at n4. For C, the products
are zero at c,n3 and `1-sqrt(2)<0` at n4. Their coefficients are positive.
Every other wall has nonpositive product by all-wall acuteness. Their
positive time coefficient and negative norm give actual future points.
Q=0 makes H/F/C purely A spatially. L clears each foreign component with
its nonnegative inverse-Gram coefficients before adding the two positive
own coefficients. Its products at c,n1 and the foreign normals are zero,
at n2 negative, and at every other actual wall nonpositive. This retains
all spectators rather than dropping their cost.

I rederived the N17 formulas rather than importing the N18 table:

* L: `a=28/51`, `d=a+3z^2>1/3`; `g^2>9/20` is equivalent to
  `2523z^2-720z+121a>0`. Replacing a by4/9 yields the stated negative
  discriminant. The inherited `8/81<z<1/3` gives timelikeness and all
  positivity conditions independently of this estimate.
* H: `A^2=2(10+7sqrt(2))/17`,
  `d=(165-12sqrt(2))/(17(11+6sqrt(2)))`. The inequalities
  `A^2>9k` and `d>1/3` reduce exactly to the displayed positive signs.
* F: `A^2=(8phi-9)/17`,
  `d=4(926-570phi)/(17(phi+1))`. The upper rational phi bound leaves
  numerator47506 in `A^2-9k` and2493 in `11095-6857phi`. In particular
  the tighter N17 density is proved, not assumed from N18.
* C: `A^2=(51+39sqrt(2))/17`,
  `d=2(119-12sqrt(2))/(17(17+12sqrt(2)))`. Its pressure and density
  reduce to `663>420sqrt(2)` and `425>276sqrt(2)`.

For H/F/C, `A^2>9k` is exactly the positive-root comparison giving
`g^2>1/2`; `g>0` follows from `U^2=A^2-k>0`. Thus every type has
strict demand squared above9/20 and density above1/3.

## Exhaustion after the first common pair

For either far common pair, the reflected inventory is valid at the
same energy threshold. Its raw support is disjoint from the initial
support. Correct reversal of the normal basis changes the raw-sum sign;
centering therefore gives **positive** product `2BB'/17`. Consequently
the normalized correlation is `<6/17`, not zero or negative. The summed
actual demands imply energy `>153/230>4/7`, covering every L/H/F/C pairing
without assuming their walls or scalar parameters are independently
realizable (candidate lines178–191).

If both far pairs are distinct, their demand is `x_A.R>17/5`. On the
first five coordinates R is constant in the form used by the candidate;
the raw-sum identity therefore gives its exact K at lines198–199 for
every foot, irrespective of individual raw coefficient signs. Each
rational square bound is valid and K<1. Thus `rho_R<sqrt(3)`. Writing
`v=x_A.U/||U||`, we have `v^2>9/20`, `v<=sqrt(e)` and
`e-v^2<17/140`. The perpendicular decomposition then gives
`x_A.R < sqrt(12/7)+sqrt(17/5)<77/24<17/5` since `R^2<28`.
This does not require an unjustified monotonicity argument. The two
far-common alternatives and the both-distinct alternative exhaust the
remaining actual walls (lines193–213).

## Acceptance boundary

No blocker or missing mathematical branch was found. The accepted new
substance is the actual A16 component estimate, and hence its union
with the existing stronger r>=17 theorem. Standard finite spherical
diagram classification, compactness excluding affine principal diagrams,
and ordinary acute compact Coxeter geometry remain the named geometric
inputs; no new script result or arithmetic classification is used.

The conditional spectator calculations at lines217–220 are correct
(`4/7+6/14=1`, `4/7+5/14+9/100=713/700>1`). They do not establish all
other dimension23 branches. A compact ceiling22 must wait for root's
separate held-source acceptance and complete forest composition.
