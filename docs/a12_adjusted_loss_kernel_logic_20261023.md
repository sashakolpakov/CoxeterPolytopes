# A12 adjusted losses: conditional scalar logic

23 October2026. **Prepared source only: no compilation or actual axiom
audit is credited.** The source is
`formal/CoxeterProofs/A12AdjustedLossKernel.lean`. It imports only Std and
changes no existing module. Ten declarations have requested axiom prints;
those requests are not execution evidence.

## Accepted mathematical input and purpose

The written theorem is
`results/arithmetic16_a12_both_typed_seven_rank_two_adjusted_loss_exclusion_candidate_20261023.md`
(a35d9f56...), with complete different-author review6a44888d and root gate
`results/review_arithmetic16_a12_both_typed_seven_root_acceptance_and_twelve_group_frontier_20261023.md`.
It excludes original exact-dimension16 group8, leaving twelve necessary
groups. The global arithmetic ceiling stays16. The older executed
thirteen-group frontier kernel retains its original scope.

This new module addresses the nontrivial CONTINUOUS scalar argument and
the simultaneous finite-list budget. It does not merely assume the
component-energy floor or encode the deletion of a row.

## Algebraic setting

The metric quantities live in a linearly ordered commutative ring R with
Lean.Grind's `OrderedRing` compatibility, not in Nat. All fractions are
cleared by positive integer constants. Radical b is specified by

    b >= 0, b²=u(1-u), 0<=u<=1.

No square-root implementation, real analysis library, field enumeration,
transcendental function, geometric object or real instance is imported.
These ring/order typeclasses and all theorem parameters are assumptions;
printed foundational axioms alone do not enumerate them.

## Derived statements

1. `circle_upper` derives the continuous Cauchy estimate from the exact
   circle identity and a strict rational square gap. For arbitrary a,c,r,
   the identity is

       [a(2u-1)+2cb]²+[c(2u-1)-2ab]²=a²+c².

   Nonnegativity of the second square gives the desired upper bound when
   r>=0 and r²>a²+c². No maximum value is a premise.
2. `radical_half` derives 2b<=1 from the primitive radical equation and
   the square(2u-1)². It supplies the first interval estimates.
3. `outer_adjusted_polynomial` proves

       264u+80u²+2000b<1201

   on the WHOLE closed interval. This is the adjusted loss with a=.212
   and beta=.08. The ranges are5u<=2 and5u>=2. On the second, the
   circle bound uses a212,c1000,r1023 and a positive square gap1585;
   the correction25(2u-u²)>=16 follows from interval factors.
4. `inner_adjusted_polynomial` proves

       698u+61u²+4000b<2402

   on the WHOLE closed interval. These are exactly the coefficients
   a=.205,beta=.0305. The proof splits at1/2 and7/10. The middle chord
   20u²<=24u-7 and circle data1928,10000,10185 leave the strict margin
   corresponding to .6003125<.6005. The final range derives
   5b<3-u from(2u-1)(13u-9)>0. Endpoints are included.
5. `individual_adjusted_loss` accepts either actual individual role:

       beta>=.08 and loss<.212u+b; OR
       beta=.0305 and loss<.205u+b.

   It DERIVES loss-(beta/2)(2u-u²)<.6005. Actual positive outer
   coefficients .15745 satisfy the first role. Inner coefficients are
   exactly .0305. The source uses each foot's own coefficient.
6. `Foot` carries those primitive interval, root and role data plus
   beta>0. `finite_gain` proves strict positivity of the total gain for
   every nonempty finite list by structural induction, and nonnegativity
   for the empty list. It uses

       credit=1000 beta²(2u-u²),
       debit=beta(2000 loss-1201).

   Their difference is1000 times the actual completed-square gain.
   There is no assumed nonadjacency or upper bound on list length.
7. `simultaneous_foot_completion` takes ONE pressure, ONE squared norm
   and ONE nonnegative completed square. In cleared form the inputs are

       2000 base-sumDebit <= 2000 pressure,
       1000 normZ <= 1000 base-sumCredit,
       0 <= energy-2pressure+normZ.

   For the empty weak list, a strict reference pressure is retained as
   an explicit geometric input. These hypotheses imply energy>base.
8. `reference_energy` calculates20000000 base=14204501 from rational
   reference heights1.8795,1.9405 and the reference norm identity.
9. `typed_seven_budget_contradiction` checks the final rational sum
   using the actual defect floor>.09 and CONDITIONAL rank-two floor>.2.
10. `complete_typed_seven_scalar_contradiction` combines the entire chain:
    individual primitive foot data, the one-vector estimates, calculated
    reference, conditional spectator floors and original norm<1. A final
    A12 energy bound or norm contradiction is not an additional premise.

## External premises and exact limit

The following mathematical obligations remain external, even after a
successful compilation and separate axiom audit:

- Actual compact acute geometry, one selected farthest vertex/vector,
  all-wall acuteness and positive-cleared feasible points.
- The complete fresh N13 wall inventory, every finite/dotted/foreign
  alternative, the actual K7 field lock, cyclic integrality, conjugate
  Schur bounds and helper norm/trace exclusion.
- Strong finite-own caps and the radical pressure/loss inequality for
  every remaining actual isolated weak cap.
- The correct reference coefficients and the fact that its norm and
  actual pressures give exactly the three simultaneous budget inputs.
  The favorable old-Gram sign for adjacent feet is part of this bridge.
- The empty-list strict reference pressure. Distinct actual component
  projections are orthogonal parts of the SAME original vector.
- Either ACTUAL typed defect N7 or O7; its own energy>.09; and its
  separately proved rank-two spectator floor>.2. No universal ordinary
  rank-two floor>.2 is assumed or asserted.
- A real/ring instance and interpretation of the continuous variables.

The finite list is algebraic. It is not an enumeration or a construction
of actual walls. No number field, polytope, dimension instance, full
geometric formalization or new example is supplied. This module does not
verify the six papers, the entire twelve-group frontier, or arithmetic15.
