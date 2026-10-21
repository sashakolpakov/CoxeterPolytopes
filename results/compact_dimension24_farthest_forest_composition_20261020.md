# Compact Coxeter dimension at most24 from actual component projections

20 October 2026. Complete global composition by root, held for
separate mathematical audit. Paper proof only: no scientific,
Lean or TeX execution. All predecessor sources are preserved.

**Theorem.** Every compact hyperbolic Coxeter polytope has dimension
at most24, conditional on the complete source-bound acceptance of
the component proofs explicitly listed below. They are internal
proved lemmas, not external literature assumptions. This source
supplies the final global accounting.

## 1. Input boundary

External inputs are the classical finite spherical diagram
classification and the published general compact ceiling29.
The farthest-vertex normalization gives a single actual spatial
vector x with squared norm<1; orthogonal finite components have
additive squared projection norms. No arithmeticity or Coxeter
property of arbitrary metric faces is assumed.

All paths in the following table are under results/:

| Input | Source and exact SHA256 |
| --- | --- |
| Ordinary small-component rates; universal short-A1/rank-two exception accounting; actual B/D>=9 rate>9/14; component caps A<=20,B<=22,D<=23 in the ambient range<=29 | `compact_forest_ordinary_a_parameter_and_dominant_classical_frontier_20261020.md`, `dd6f8b0de2d7d0bc835be3eb4e320013aa40792407fb2a747f35bcf04385d5ad` |
| Actual A_r rate>1/2 for every r>=7 | `compact_large_an_actual_half_uniform_transfer_20261020.md`, `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7` |
| Actual A_r rate>5/7 for every r>=17 | `compact_large_an_actual_projection_five_sevenths_20261020.md`, `7a57b66e1b9b68948c6fe14c493527ff3ac0b963144f904ff512ead723e192b0` |
| Generic localized B/D prefix theorem,6<=M<=16 | `compact_classical_generic_prefix_energy_20261020.md`, `f637d7294d17e37581330ca508c914c23a08ce259845614b055e45815d2790ef` |

The generic theorem supplies in particular:

    B>=19 or D>=20: component energy>.90715>.90;
    B>=21 or D>=22: component energy>.9401975>.94.

The first input is used for its component conclusions and complete
exception proof, not for an assumed ordinary-A parameter. There
is no dependence on the stronger standalone ordinary A8 theorem,
the A20 four-sevenths lemma, or the predecessor global26 conclusion.

## 2. The dominant component and the exact permitted debit

Suppose25<=n<=29. Call an A component of rank>=7 or a B/D component
of rank>=9 dominant. Without one, all ranks pay their strict ordinary
rank/14 rate except at most one short A1 or instead one isolated
rank-two7 block of energy>9/100. Discarding the possible short A1
gives total>(n-1)/14>1; with the rank-two defect the bound is
>9/100+(n-2)/14>1. Both are impossible, so a dominant component
exists.

Every dominant component pays strictly more than1/2. Two would
already contradict||x||^2<1. Thus it is unique. Only after this
step may its spectators be declared ordinary apart from the same
one universal exception. Write its rank as r and their rank as
S=n-r. Their three possible costs are

    S/14,  (S-1)/14,  9/100+(S-2)/14.             (1)

Empty ordinary sums are zero; strictness in the S=2 sole-defect
case comes from its own>9/100 bound.
The exception alternatives are incompatible with each other,
independently of the dominant type, by the pinned actual-wall
isolation and native product gap. A possible short type J may be
discarded conservatively as one rank; no extra funding is assumed.

## 3. A dominant A component is impossible

Its uniform half bound and(1) first give S<=7. Indeed S>=7
contradicts the norm in the fully ordinary case; S>=8 does so
after either debit. Consequently r=n-S>=25-7=18. The stronger
A>=17 theorem therefore applies and pays>5/7.

Apply the same three cases again. Without a debit, S>=4 gives
5/7+4/14=1. With a short A1, S>=5 leaves at least four ordinary
ranks and gives the same contradiction. With one rank-two defect,
S>=5 gives

    5/7+9/100+3/14=713/700>1.

Therefore S<=4 in every case. This forces r=n-S>=21, contrary
to the accepted arbitrary-spectator A cap20. No dominant A case
survives for25<=n<=29.

## 4. A dominant B or D component is impossible

Its uniform actual9/14 rate and(1) give S<=5. Hence r>=20.
The accepted component caps give r<=22 for B and r<=23 for D.

If r20 or21, the generic prefix rate is>.90 and S=n-r>=4.
Even after discarding a short A1, at least three ordinary ranks
pay>3/14>.10. With a rank-two defect the cost is
>9/100+2/14>.10. The fully ordinary case is stronger. Thus the
total squared norm exceeds1.

If r22 or23, the generic prefix rate is>.94. Also
S=n-r>=25-23=2. After a short-A1 discard at least one ordinary
rank pays>1/14>.06. With a rank-two defect, that block alone
pays>9/100>.06. Again the fully ordinary case is stronger, so
the total norm exceeds1. These exhaust the possible B/D ranks.

Thus no dominant type is possible, a contradiction. Dimensions
25 through29 are excluded. The external compact ceiling29 then
proves dimension at most24.

The theorem concerns arbitrary compact Coxeter polytopes. It
does not alter the already stronger internal arithmetic compact
ceiling20, nor produce a new example. Compact examples remain
known through dimension8.

## 5. The next boundary, as a necessary condition only

The same estimates leave only two possible dominant configurations
for a hypothetical dimension24 polytope:

* A20 with spectator rank4, containing a short A1 or one isolated
  rank-two7 defect. The fully ordinary case costs>5/7+4/14=1.
* D23 with one isolated short A1 of type I. A single ordinary A1
  already costs more than the remaining.06; type J requires two
  A1 components and cannot fit.

For the second assertion, S<=5 first gives r>=19. The generic
M13 floor>.88864 removes r19 with S5, M14/M15 floors remove
r20/21, and the M16 floor removes r22 with S2, after every
allowed debit. B23 and larger D are capped. This is a necessary
farthest-diagram condition, not an exclusion of the two remaining
configurations or a theorem of dimension at most23.

All component inequalities concern orthogonal projections of one
actual vector. The proof assumes no rank/14 rate for a large
component, no census of all hyperbolic polytopes, and no arbitrary
metric-face inheritance. Written proofs, exact source integrity,
script execution and conditional Lean coverage must be reported
separately. Further lowering and a compact dimension9 construction
remain active research objectives.
