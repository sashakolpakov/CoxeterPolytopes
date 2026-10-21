# Compact arithmetic Q(√5): whole-component hub weights

The [source](../scripts/fresh_sqrt5_dimension23_coverage_20261015.py) and [complete certificate](fresh_sqrt5_dimension23_coverage_20261015.json) pass generation and saved replay. They bind the completed preceding coverage, the final independently reviewed finite queue and the independent whole-component weight review.

The certified statement is precise: a compact arithmetic hyperbolic Coxeter polytope with ground field exactly Q(√5) has dimension at most 23. The [preceding binder](fresh_sqrt5_dimension24_coverage_20261015.md) excludes dimensions 25–29. This separate composition excludes dimension 24 using the [whole-I2(5) rule](fresh_sqrt5_dimension24_i25_cross_20261015.md) and its [independent exact review](fresh_sqrt5_dimension24_i25_cross_20261015_review.md). No previous certificate is changed.

## Exact whole-vertex counts

Let a be the number of entire A1 components and b the number of entire I2(5) components in a spherical vertex diagram. Connected-pair bonuses and the normalized A1 cross bonuses depend on a and on whether b is zero. A disconnected A1-to-component angle receives u/a, as before.

The additional cross family joins a root of an entire I2(5) to a root in another entire component, whose type is neither A1 nor I2(5). Such an angle receives v_a/(2b). For each destination root there are exactly 2b source roots, so their sum is v_a. A1-to-I2(5) pairs belong only to the A1 family; I2(5)-to-I2(5) pairs receive no cross bonus. Thus the two families do not overlap, and there is no double counting. The root-orbit multiplicities and both denominator cancellations are verified independently.

For fixed exact (a,b), the remaining forest contains neither of these two component types. The independent type-first dynamic program sums its additive component costs, then adds the a A1 and b I2(5) costs. All 157 feasible contexts in rank 24 are included, with all 77 allowed finite component types. The exact maximum is

`267527/5000 = 53.5054 < 54`,

with margin `2473/5000`. This is a feasible upper bound, not an optimality claim. Floating-point optimization only proposed the rule; the displayed bound and all payments use exact rational arithmetic after upward grid rounding.

## Face payments and arbitrary spectators

The short-path and entire-H4 baseline is unchanged. All bonuses are nonnegative. The independently reviewed 76 triangle markings still receive at least two, with the same spherical-extension stability. Any persistent short/H4 quadrilateral payment is also retained.

For a k-node zero-touch core, the external rank in dimension 24 is 26−k. The binder uses one common pair (e,f) of external A1 and I2(5) counts at all four corners. It checks every nonnegative pair with `e+2f≤26−k` and remaining rank zero or at least two. The remaining forest has neither A1 nor I2(5); A2, B2 and A-type components show that every rank at least two is possible for this conservative enumeration. A remainder of one is impossible. Each corner's own counts are added before evaluating its presence flag and denominators.

The [complete-weight independent review](fresh_sqrt5_dimension24_i25_cross_20261015_review.md) directly checks 2,516 prior-catalogue/full-finite-queue rows in 65,300 spectator cases. It includes all 1,356 rank-feasible finite singular stopping states from the complete 2,070-obligation list. Expanded and previously paid states are not omitted. Its minimum face slack is 1/170000.

The current binder additionally replays every older coverage branch under this new rule: all original sealed cores, all 147 singular first attachments, their representative/subcore maps, the 51-state singular frontier, all seventeen continuation-closure states, and every rank-feasible dotted-queue singular state. Every finite-queue singular row is also reevaluated directly. An over-rank row is skipped only for the dimension-24 payment; it remains part of the complete arithmetic queue and the preceding higher-dimensional coverage. No old long-pair payment is assumed to transfer.

The arithmetic and geometric completeness argument is unchanged: conditional sealing makes marked components entire unless a touching root is present; every such root is followed through the reviewed exact domains or creates a persistent baseline payment. Orthogonal whites are precisely the spectators just enumerated. Nonsingular no-touch cores are rank-excluded. The complete finite and dotted queues use the justified 31-node ceiling for dimensions at most 29, which also covers dimension 24. All shared-white minimal unions, exact dotted parameters, literal root/child maps and singular stopping alternatives remain source-bound dependencies of the preceding proof.

## Strict dimension-24 contradiction

With every triangle paid two, every quadrilateral one and all weights nonnegative, double counting gives average vertex weight at least `5 f₂/f₀ − binom(24,2)`. The symmetric nonnegative h-vector of a simple compact polytope gives `f₂/f₀>66`: the minimum paired coefficient is 66, while h₀=h₂₄=1 contributes a strictly larger extreme coefficient. Hence the average vertex weight is strictly greater than `5·66−276=54`.

This contradicts the exact upper bound 267527/5000. Together with the completed preceding coverage and the known compact bound 29, it yields the stated field-specific conclusion. Other fields and noncompact or nonarithmetic classes require separate arguments; this certificate makes no claim about them or about priority.

```sh
python3 scripts/fresh_sqrt5_dimension23_coverage_20261015.py
python3 scripts/fresh_sqrt5_dimension23_coverage_20261015.py --verify
```

The final certificate records 2661 zero-touch rows, of which 1897 fit dimension 24, and checks 60735 common spectator-count cases. It verifies 471 literal maps and 73 recursive manifests containing 722 source-hash references. Generation, saved replay, syntax and local-link checks pass.
