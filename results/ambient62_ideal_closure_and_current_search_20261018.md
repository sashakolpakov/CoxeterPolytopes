# Dimension 61, positive ideal margins, and the current searches

The completed unrestricted dimension bound is **61** for finite-sided,
finite-volume real hyperbolic Coxeter polytopes. Its
[proof](fresh_rank15_rank16_dimension61_20261018.md) is complete
within the recorded internally checked computer-assisted framework, with no
identified missing mathematical step. The internal reviews are agent reviews;
external refereeing and proof-assistant verification are absent. The
[163 status note](dimension163_proof_status_20261017.md) makes the same
distinction for the earlier bound. The gap to known examples is a reason to
strengthen these results, not an unresolved step in those completed proofs.

This checkpoint resolves every ideal-vertex obstruction at ambient dimension
62 and lowers the complete finite maximum from 441 through 438 to 434,
equal to finite curvature. Every ideal vertex has surplus at least 73,
and the global extremum correction is 4805. The resulting contradiction
is

    0 >= 0*v_f + 73*v_i + 4805 > 0.

The previous bound 62 excludes every higher dimension, completing bound 61.
Compact arithmetic 21 and arbitrary compact 29 are unchanged. No new
finite-volume example is established. The
[preceding checkpoint](dimension62_three_transfers_raw39_and_wall59_20261018.md)
and all its artifacts remain frozen.

## Complete ideal positivity at dimension 62

Two new local theorems remove long-pair payments without adding any finite
coupon. The [C14 theorem](fresh_order14_c_inner_endpoint_projection_transfer_20261018.md)
removes the inner-endpoint pair of a whole affine C13: all six affected
maps require a physical norm-two terminal whose forced projection is too
large. Its cost falls from 84 to 83.

The [B14 theorem](fresh_order14_b_terminal_residual_transfer_20261018.md)
removes the two endpoint/near-terminal pairs of a whole affine B13. Of its
108 affected maps, 92 have a direct projection obstruction. In the remaining
16, the doubled terminal bond contradicts the Cauchy inequality in the
positive-definite orthogonal complement of the Lorentzian core. All other
possible contacts are included in that argument. Its cost falls from 84
to 82. Both theorems preserve the independently established short-path
alternative and all triangle payments.

The [composition and scalar certificate](fresh_order14_ambient62_ideal_margin_closure_20261018.md)
combines them with the earlier D13, B16 and D17 transfers. The four order-14
family costs are now A70/B82/C83/D82. Every one of the 2,641 terminal
rank-61 cusp profiles has positive margin. The normalized minimum is
2741/371712 at (10,10,11,15,15); the absolute minimum is 73 at the single
rank-61 factor. Only the proved concentration of ranks below 9 is used,
so arbitrary repeated larger factors remain in the complete domain.
Consequently every actual ideal vertex in dimension 62 has surplus at
least 73. This excludes the all-ideal case, with the completed unrestricted
bound supplying the higher-dimensional tail.

## The first finite improvement: maximum 438

The [new finite D17 theorem](fresh_finite_d17_five_unit_reduction_20261018.md)
reduces its cost from 125 to 120. Two long pairs become zero and six become
halves. The six affected selected contexts have a radical that forces any
extra supporting root either to give a retained short path or to be
orthogonal to the entire core. In the latter case a distinct existing
whole-small/D12 half pays the other half. No new coupon is introduced.
All 216 inherited literal interfaces are debited simultaneously under the
current ideal rule, and the complete nonselected interfaces are retained.

The [complete rank-62 count](fresh_finite_d17_cost120_rank62_budget_20261018.md)
includes the current D13 coupon increments and E6/D11 total 4. It checks
35,721 simultaneous type comparisons, all 1,227 large-component vectors,
and the complete residual Bellman certificates. The actual formal forest
2D15+2D16 attains 438. Curvature is 434, so a deficit of four remains.
Formal attainment is not a claim that this forest occurs in a polytope.

The five offending representative alternatives are 2D15+2D16, A1+3D15+D16,
3D15+D17, 2D16+D30 and D16+D46. Their actual type ties must also be covered;
in particular A15, B15 and D15 all have current cost 105, and A16 and D16
both have cost 114. Lowering only a representative's cost would not prove
the corresponding bound for the actual types.

The finite and ideal changes compose: the D17 theorem explicitly preserves
the B14/C14 short-path escapes, and their projection exclusions consume no
changed finite D17 long pair. Their combined necessary inequality is

    4*v_f >= 4805 + 73*v_i.

This intermediate inequality is not a contradiction for arbitrary vertex
counts. The next two completed local theorems resolve that obstruction.

## Rank15 and rank16 transfers close the finite comparison

The [complete rank15 inventory](fresh_finite_rank15_long_pair_support_scout_20261018.md)
identifies the remaining support and literal obligations. The resulting
[rank15 theorem](fresh_finite_rank15_existing_half_transfer_20261018.md)
lowers A15, B15 and D15 from 105 to 104. Two A15 radical cases use a
distinct existing whole-A2/A13 endpoint half. Four D15 cases use whole-core
projection sealing and a distinct retained whole-D16 half. B15's removed
endpoint occurs in no selected or literal case; the theorem also covers
every nonselected, universal and later interface. No coupon is added.

The [rank16 exchange](fresh_rank16_existing_payment_exchange_20261018.md)
lowers A16 to 113 and D16 to 225/2, while restoring two A17 halves and
raising A17 to 115. The A16 alternatives use those restored whole-A17
halves; the D16 alternatives use existing whole-small/D13 halves. The
rank15 theorem's D16 pairs (1,14)/(1,15) are unchanged by this exchange.
Both changes are applied simultaneously to all 216 literal records.

The [complete successor count](fresh_finite_rank15_rank16_rank62_budget_20261018.md)
uses a conservative D16 allowance of 113 to dominate A16. All changes lie
in the fixed component counts; every residual item and all 921 certified
residual tables remain unchanged. The doubled fixed-count shift is
-2*c15-2*c16+2*a17. Checking all 1,227 alternatives gives upper maximum
434, attained by the actual formal forest D17+D45 at 120+314. The A17
increase is fully counted.

The [final composition](fresh_rank15_rank16_dimension61_20261018.md)
checks the independent ideal priorities, preserved short-path routes,
disjoint selected cases, mutually retained alternative payments and
simultaneous literal debits. Finite surplus is zero; ideal surplus is
at least 73. The single positive extremum constant 4805 then excludes
dimension 62, including its compact case. There is no remaining proposed
local reduction in this proof.

## A new global arithmetic solve

The [enlarged compact arithmetic H13 model](fresh_cross_vertex_current_full_model_20261018.md)
over Q(sqrt5) incorporates 213 actual exports from twelve proved packages,
including their cross-vertex incidence costs. Positive scaling identifies
209 distinct columns, 14 of which charge two vertex types. The full model
has 439,619 variables, 82,806 face obligations and 6,558 vertex budgets.
All old matrix data are preserved.

One bounded solve, followed by exact primal and dual recovery, gives

    12.35984031684756 <= finite-model optimum <= 12.359855171582.

Both endpoints exceed the dimension-13 target 12. The previous optimized
lower endpoint was 12.3685672578. The old fixed-face-mass repair below 12
does not establish feasibility at the newly optimized face masses. The
dual verifier reads every vertex cost of a column and checks all 439,619
inequalities; it does not reduce a two-vertex column to a single owner.

A separate full-matrix replay verifies both exact certificates without
running the optimizer. A separate stdlib saved checker verifies all new
export origins, primitive scaling, both-owner costs and new-column
certificate arithmetic. Its scope does not independently reconstruct the
old matrix. Complete implicit charge cones, faces, finite forests and
ground fields remain outside this finite model. No arithmetic dimension
bound changes.

## A changed construction seam still fails to close

The [wall66 seam subdivision](fresh_wall66_seam_subdivision_20261018.md)
changes incidence by inserting a bridge J between the dotted companions
X and Y. It samples sixteen representative prefixes, obtains 95 bridge
caps and 2,685 subsequent Z caps, and checks every required finite edge,
corner and whole cap boundary. No disk closes: Y/Z and J/Z are dotted in
every sampled completion. The closest J/Z product is approximately
-2.034872, still beyond the parallel value -2 in this normalization.

Every resulting target cone retains nine open ends. The reduction from
the original ten was caused by deleting old walls, with a debit of seven
finite supports and seven compact-face certificates. These have not been
restored; the lost prism remains an obligation. The saved checker verifies
the exact old debit, selection, disk links and all final failures, but does
not independently reconstruct the numerical normals or cones. These are
bounded numerical experiments, not exact exclusions of continuous families
or an ambient finite-volume chamber.

## Validation

The completed packets have generation, separate fresh replay and internal
source/report reviews at their recorded scopes. The
[aggregate audit](ambient62_ideal_closure_and_current_search_20261018_audit.json)
records their immutable hashes and run identifiers. Its
[program](../scripts/check_ambient62_current_search_checkpoint_20261018.py)
checks file integrity, syntax and links; it does not add a mathematical
or formal proof check. No existing frozen artifact is changed, and the
example-search objective remains active.
