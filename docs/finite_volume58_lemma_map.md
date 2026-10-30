# General finite-volume bound58: proof and formalization map

For current reading, see the [current dimension-bounds manuscript](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

This map concerns **every finite-sided, finite-volume real hyperbolic Coxeter
polytope**. The completed internal theorem is `dimension ≤ 58`. It imposes
neither arithmeticity nor compactness, a simple-cusp condition or a facet
bound. The separate ambient58 research frontier is an attempt to prove57;
its current payment inequalities are not dependencies of the completed58
proof. Construction and compact-arithmetic experiments are also independent.

The authoritative endpoint is
[the dimension58 composition](../results/fresh_a15_d15_dimension58_20261018.md).
The supporting argument develops this proof
spine. [Reproduction instructions](finite_volume58_reproduce.md) distinguish
arithmetic replay, source review and formal verification.

## Exact dependency spine

|Layer|Internal statement|Frozen implementation/report|
|---|---|---|
|Final contradiction|No dimension59 object satisfies all preceding inequalities; completed59 excludes dimensions≥60|[dimension58 binder](../scripts/fresh_a15_d15_dimension58_20261018.py)|
|Finite arithmetic|Every actual finite rank59 Coxeter forest has payment≤783/2|[complete finite closure](../results/fresh_a15_d15_ambient59_forest_closure_20261018.md)|
|Last local changes|A15 `(0,12),(2,14)` and D15 `(0,11)` halves become zero, with four explicitly counted ordinary coupons|[A15 theorem](../results/fresh_finite_a15_three_coupon_zero_transfer_20261018.md), [D15 theorem](../results/fresh_finite_d15_a4_a11_half_zero_transfer_20261018.md)|
|Simultaneous baseline|Actual-corner priority, finite6, all prior finite/ideal transfers,1890 positive cusp profiles|[A14/C-D16 baseline](../results/fresh_a14_cd16_ambient59_frontier_20261018.md)|
|Ideal multiplicity|A marked pair in an affine rank r factor occurs V/(r+1) times at its actual cusp|[factor-fiber theorem](../results/fresh_cusp_factor_payment_framework_20261017.md)|
|Ideal scalar reduction|Low-rank concentration preserves normalized cost, decreases normalized margin and V|[all-rank convex-order proof](../results/fresh_distance14_unrestricted_dimension156_20261017.md), used only on ranks below9 in this endpoint|
|Global geometry|The two-face deficit is bounded below by finite and ideal shape curvatures plus the two-extrema correction|[shape-curvature proof](../results/fresh_cusp_link_shape_curvature_20261017.md), strengthened here using every h coefficient≥1|
|Higher-dimensional tail|All dimensions≥60 are already excluded|[completed59 theorem](../results/fresh_depth6_dimension59_20261018.md)|

Each frozen JSON has a `sources` map. This is the authoritative recursive
dependency graph; a direct link in this table is not a claim that all older
local proofs have disappeared. In particular, the completed59 tail is an
**internal theorem to translate**, not a published assumption. The baseline
contains an extensive but finite local proof inventory. The compact spine
describes its interfaces, not an automatic compression of all geometry into
a hash check.

## Quantified internal lemma interfaces

The following are proof obligations for the internal development. Their
labels are suggested Lean/module names, not declarations of axioms.

1. **`ActualCornerSupport`.** For every actual two-face F of a finite-volume
   Coxeter n-polytope, its n−2 supporting roots W are positive definite;
   adjoining one side gives a positive-definite edge principal. At each
   ideal corner, exactly one irreducible affine factor varies, and the
   other factor pieces are elliptic. This must allow nonsimple cusps.
2. **`AffineFactorFiber`.** For any ideal vertex with affine factor ranks
   `r : List Nat`, all positive and summing to n−1, put
   `V = ∏ (r_i+1)`. For every marked pair of the i-th full affine factor,
   the number of actual two-face corners representing it is exactly
   `V/(r_i+1)`. A connected affine subdiagram used in a local proof must
   first be shown to be the whole actual factor.
3. **`LorentzRadicalSeal`.** Let a norm-two acute core span a subspace
   containing a timelike vector. Every Gram radical is an actual zero
   linear combination of its roots. If all permitted nonzero contacts of
   an additional supporting white must meet two opposite radical signs,
   either the specified retained short path exists or every such contact
   vanishes. The latter conclusion seals all paying components, including
   unmarked witness factors; it is not merely an isolated-node test.
4. **`LorentzProjectionExclusion`.** For every allowed contact vector q
   into a nondegenerate Lorentz core G, a prospective norm-two root has
   positive-definite residual square `2−qᵀG⁻¹q`. A projection greater than2
   is impossible. Equality forces zero residual and exact products with
   every other core-span root. This supplies the saturated-root exclusions.
5. **`CompleteSelectedSupport`.** Every remaining actual quadrilateral,
   after the independent ideal/short-path branches and complete low,
   disjoint, protected and adaptive routes, lies in the enumerated native
   selected domain with its actual marked component embedded injectively
   and inducedly. The order15 domains have536 patterns each,482 A15 maps
   and870 D15 maps; failed searches are part of completeness. Finite saved
   corners alone do not establish this universal statement.
6. **`LocalPayment59`.** There exists one nonnegative intrinsic corner
   weight rule w such that every triangular two-face has total≥2 and
   every quadrilateral total≥1. Its finite cost is exactly the typed rule
   exported below, and its affine factor costs are bounded by c(r).
   The proof uses `ActualCornerSupport`, the local seals/exclusions and
   a well-founded priority order. Every later replacement must retain the
   earlier chosen output or supply a new independent replacement. It must
   not be assumed that all short ideal paths themselves have unit weight.
7. **`TypedForestDomination59`.** For every finite Coxeter forest of rank59,
   replacing each component by its specified representative dominates its
   base, every pairwise interaction and the two guarded context charges.
   This includes exceptional types and arbitrary-label rank-two factors.
   There are180 actual type representatives and32,400 ordered pair checks.
8. **`ForestDomain59`.** Every representative forest decomposes into one
   of961 vectors of fixed ranks9–16 and fixed A17, one separately fixed
   rank5 multiplicity, and residual ranks1–4,6–8,17–59. Rank17 residual
   means D17, not a second copy of fixed A17. There are2574 expanded rows.
   Residual pairwise interactions vanish. Their exact items include all
   interactions with the fixed components.
9. **`LastItemUpper`.** For a finite list of positive-rank items `(r,w)`
   and table V with V(0)=0, if `V(n)≥V(n−r)+w` for every eligible item and
   n, then every item multiset of total rank n has cost≤V(n). Induction on
   the number of items proves this. An attaining predecessor is needed
   only for sharpness, not for an upper bound. Unreachable ranks must be
   represented and verified if present; these completed59 tables contain
   rank1 and hence every rank is reachable.
10. **`ContextCover59`.** Every residual forest is covered by the no-context,
    reserve-one-A1, or reserve-two-A1 branch, with all A2 multiplicities
    fixed explicitly and no-A2 residual tables. The old charge is
    `a2(a2−1)/4` when an A10 and at least two A1 occur; the newer charge is
    `a2*d5` when an A9 and at least one A1 occur. Upper guards allow the
    corresponding rank9/rank10 representative. If both guards hold, the
    two-A1 branch includes both charges. Presence of several witnesses
    does not multiply either charge.
11. **`FiniteBudget59`.** The complete domain, domination, context cover
    and5148 Bellman certificates prove charge≤6264/16=783/2 for every
    finite vertex. This is a representative upper, not an asserted exact
    actual maximum. The representative maximizer has actual cost355.5.
12. **`ProductPolynomial`.** For every positive rank list r summing to58,
    `H=∏(1+…+t^{r_i})`, `V=∏(r_i+1)=Σh_j` and
    `G=Σ min(j,58−j)h_j`. All coefficients are positive, symmetric and
    at least1. The curvature is `(5G−58V)/2` and the envelope is
    `V Σc(r_i)`, not merely `Σc(r_i)`.
13. **`LowNineConcentration`.** If `1≤a≤b<9`, the replacement
    `(a,b)→(a−1,b+1)` at fixed total rank increases the centered sum of
    independent discrete uniforms in convex order. Therefore it does not
    increase normalized curvature. Since c(r)=r/2 for r≤9, it preserves
    normalized cost. It decreases `(a+1)(b+1)` to `a(b+2)` and increases
    the sum of squared ranks, so iteration terminates at a profile with
    at most one positive rank below9. Other ranks are not concentrated.
14. **`TerminalIdeal59`.** Every one of the1890 terminal profiles has
    normalized envelope margin≥1/340>0 and absolute envelope surplus≥91/2.
    This is verified from exact polynomial products and rational c(r),
    with a proof that the profile list is the entire terminal domain.
15. **`PositiveSurplusTransfer`.** If `V≥V'≥0` and `m≥m'>0`, then
    `Vm≥V'm'`. Apply this only after terminal positivity is established.
    Together with the factor-fiber bound `C_actual≤VΣc(r_i)`, every actual
    ideal vertex has `κ−C_actual≥91/2`. Negative terminal margins would
    not permit this argument.
16. **`ShapeCurvature59`.** Let `D5=Σ_F(5−s(F))`. The local index formula
    and global minimum/maximum vertices give
    `D5≥(783/2)v_f+Σκ(v)+4350`. The extra4350 is counted once. The two
    extrema give `5 Σ_{j=30}^{58}(2j−58)=5*870`; every h coefficient is≥1
    for the simplex/product-of-simplices vertex figures here.
17. **`PaymentUpper`.** From nonnegative weights and the triangle/quad
    lower totals, `D5≤Σ corner weights≤(783/2)v_f+ΣC(v)`. Faces with at
    least five sides have nonpositive deficit and need no positive lower
    payment. This is a finite-sum double-counting theorem.
18. **`Exclude59`.** Subtract the last two inequalities and use the cusp
    surplus: `0≥4350+(91/2)v_i>0`, contradiction. This includes v_i=0.
19. **`Dimension58`.** `Exclude59` plus the internally proved completed59
    tail implies n≤58. It would be incorrect to infer the higher tail
    merely from excluding dimension59 or to assume induced faces inherit
    every required Coxeter/volume property without a theorem.

## Published mathematical inputs versus internal proofs

The external boundary must contain exact statements and hypotheses from
published sources, not names standing for the repository's new results.
The frozen reports identify these primary inputs:

* Classification of irreducible finite and affine Coxeter diagrams,
  including the arbitrary-label rank-two families and proper-principal
  definiteness properties.
* Compact Lannér order at most5 and quasi-Lannér order at most10. The
  [high-rank census](../results/fresh_high_rank_quasilanner_census_20261017.md)
  records Chein's primary classification and distinguishes the published
  ceiling from the new146,262-attachment reproduction of26 diagrams.
* Standard finite-volume cusp-link and simple-at-edges structure, and
  published containment of a parabolic subdiagram in its ideal vertex.
  The repository also gives its own positive-kernel/acute-cone argument.
* Khovanskii's exact incidence reduction and simple-polytope index formula
  (theorems2,6,9 as cited by the shape-curvature report).

The stronger4350 correction, factor-fiber payment multiplicity, convex-order
reduction, all current local weights, arbitrary-support seals, complete
native-domain arguments, finite contextual budgets, positive ideal margin
extension and both completed58/59 conclusions are **internal**. They must
be proved or reduced to checked certificates; none belongs in an external
axiom file merely because it is lengthy.

## Plain certificate boundary for Lean

The [exporter](../scripts/export_finite_volume58_arithmetic_20261018.py)
writes plain data.
Rationals are numerator/positive-denominator pairs, finite costs are integer
sixteenths, and type maps use integer indices. The export contains:

* all180 actual types, their representative map and actual/upper base and
  pair-cost matrices;
* all5148 residual item/rank/weight/value/predecessor tables;
* all2574 fixed rows and every contextual option, with fixed costs;
* all1890 factor-rank lists, h-polynomial coefficients and exact moments;
* all58 normalized ideal envelope values, source hashes and final constants.

A hash establishes identity, not mathematical validity. A kernel checker
must reconstruct each domain, verify item identities and inequalities,
check the context-cover theorem and polynomial products, and derive the
upper/margin statements. The export deliberately labels itself
`EXPORTED_FROZEN_DATA_NOT_A_FORMAL_PROOF`. Geometry and the higher tail
are excluded from the arithmetic data boundary and remain separately
named internal proof obligations.

## Formalization status

The parent reports a building Lean4.32.2 project with kernel-checked generic
last-item induction, selected-deficit budgeting, integer-gap and tail
composition lemmas, without `sorry` or custom axioms in those modules.
This is a partial verification milestone. This document does not assert
that all5148 tables, all1890 polynomial certificates, the geometric
interpretation, every local theorem or the higher-dimensional tail have
already been connected in Lean. End-to-end kernel verification of58
remains pending; the current completed theorem is the internally reviewed
computer-assisted mathematical proof described by its frozen endpoint.
