# Independent review of the exact finite-corner budget

The [exact finite-forest budget](fresh_exact_finite8_forest_budget_20261017.md) passes a complete mathematical/source read and a separate certificate check with no producer arithmetic imports. The current unchanged rule has exact abstract-forest maximum **548 at rank69** and **553 at rank70**. These are universal upper bounds for actual finite vertices; attaining abstract diagrams do not assert ambient realizability. No new face coverage or dimension bound follows from this budget alone.

The [review source](../scripts/fresh_exact_finite8_forest_budget_20261017_review.py) independently checks the type domination, complete large-component domain, every Bellman alternative and an attaining predecessor, and every global maximum and attaining forest. Final generation and a separate fresh saved replay both terminated with exit0/PASS; the [saved certificate](fresh_exact_finite8_forest_budget_20261017_review.json) binds the frozen primary source, data and report.

## Complete type domination and coupon count

The finite Coxeter classification is inherited: the only unbounded-rank connected families are A, B and D; the exceptional types have rank at most8, and arbitrary I2(m) has rank2. The whole-component rule pays every internal pair through rank17. At larger rank its inherited exact base costs are8r−34 for A/B and8r−32 for D. These all-rank path/fork identities were read as proofs, not inferred from a finite sample.

The reviewer independently counts twice the cross-payment between any two component types. A small component of rank r≤5 contributes3r against D11–D16, r against A13/A15, and2r against A17. Add2 for each whole-component pair A1/D12, A2/B10 or A4/D13. The first addition upgrades two short-leaf halves, the second adds two A2-mark halves, and the third adds two A4-endpoint halves at the distinct D13 long-arm neighbor. These roles are disjoint from the appropriate previously counted pairs.

This verifies the claimed exact representative-forest formula. It also proves componentwise domination: replace ranks1–9 by A, rank10 by B10, ranks11–16 by D, rank17 by A17 and larger ranks by D. Internal cost cannot decrease, nor can the cross-payment against any other component. Because all interactions are pairwise and nonnegative, replacements one at a time preserve the inequality for an entire forest. Every representative is itself a finite Coxeter type, so domination gives both a valid upper search space and attainable abstract maximizers.

As an additional independent check, the reviewer constructs228 type profiles through rank75 and verifies all51,984 ordered partner comparisons. One symbolic rank-two profile covers every I2(m),m≥4; I2(3)=A2 is handled separately because of the special A2/B10 coupon. The symbolic profile's cost and broad small-component eligibility do not depend on m. This check is not a new proof of finite Coxeter classification; that classification and the inherited same-component rule remain mathematical inputs.

## Bellman certificates prove every residual optimum

The producer enumerates large-rank counts successively. The reviewer instead enumerates the number of large components and all nondecreasing rank tuples from10 through17 with total at most75. Exact set equality verifies all1,608 count vectors, their ranks, base costs and377 distinct coupon keys.

For each coupon key, the remaining ranks are1–9 or18–75. Their costs are additive once the large counts are fixed. The reviewer derives each residual cost by summing its independently counted interactions with a representative large-component multiset, rather than calling the producer's residual formula.

For every saved state s and every allowed component rank r≤s, it checks

    saved_F(s) >= saved_F(s-r) + independently_counted_cost(r).

It also checks that the saved last rank is permitted and attains equality. There are **885,950 upper-bound alternatives** and **28,275 attaining predecessors** across377 tables. Induction from F(0)=0 proves both directions: no residual forest exceeds the saved value, and following positive-rank predecessors constructs a forest attaining it. This is certificate verification, not a heuristic optimizer or a sample of forests.

The reviewer then checks **22,504 global alternatives**: for every total rank and every feasible large count vector, the reported optimum dominates its exact large cost plus the residual-table value. It independently counts the full pair cost of each reported maximizing forest and verifies equality with the selected alternative. All75 positive-rank maxima are thus certified; the zero-rank base case is also checked.

## Consequences and limits

At rank69 the witness is `A1+4 A17`. Four A17 components cost544 internally; their existing half-coupons against A1 add4, giving548. The finite curvature coefficient is544, so this exact abstract-forest budget still exceeds it by4. Improving only the estimate of this same rule cannot eliminate that abstract obstruction. It does not follow that this vertex diagram occurs in an actual finite-volume polytope.

At rank70 the witness is `A2+4 A17`, with cost553. The finite curvature coefficient is560, leaving slack7. This refinement may be combined with a separately proved ideal-vertex and face-payment argument; this package does not supply those arguments.

The unchanged four coupon packages are bound directly: the [base leaf-nine rule](fresh_finite8_leaf9_half_coupon_budget_20261017.md), [broad A17/D11 rule](fresh_finite8_a17_d11_coupon_budget_20261017.md), [A4/D13-neighbor addition](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.md), and [D12/B10 additions](fresh_finite8_d12_b10_coupon_budget_20261017.md). No geometric payment is added or removed here.

## Verification status

The independent complete comparison, final hash-bound generation and separate fresh saved replay all terminated with exit0/PASS. The full producer source and proof, together with the inherited coupon definitions and all-rank base-count arguments, received a mathematical read with no identified gap. A separate root review read the complete independent checker and its mathematical argument, including the alternate component-count enumeration, coupon type domination, Bellman upper bounds and attaining predecessors, and all global bounds and witnesses; that review did not duplicate this certificate arithmetic. The 16 direct source bindings, report links and Python syntax were checked. No full raw-forest enumeration, physical vertex realization, external refereeing or proof-assistant formalization is claimed.
