# Complete Q(√5) bound-21 composition

The computer-assisted proof candidate gives **dimension at most 21** for
compact arithmetic hyperbolic Coxeter polytopes with ground field Q(√5).
The [binder](../scripts/fresh_sqrt5_h3_hub_coverage_20261016.py) and
[certificate](fresh_sqrt5_h3_hub_coverage_20261016_bound21.json) combine the
completed [dimension-23 exclusion](fresh_sqrt5_hub_count_coverage_20261016.md)
and all preceding higher-dimensional exclusions with the new dimension-22
rule. The statement remains subject to external mathematical validation;
no literature-priority claim is made.

Generation and fresh saved replay of the final certificate both pass.

## The additional H3 payment

The [independent exact review](fresh_sqrt5_dimension22_h3_cross_20261016_review.md)
checks the literal graph rule without producer or optimizer imports. Let
`a,b,c` be the counts of entire A1, I2(5), H3 components in a finite vertex
diagram. The coefficient contexts use exact `a`, `min(b,2)`, and whether
`c>0`. The angle denominators always use actual source-root counts.

- A1-to-root angles divide their coefficient by `a`; destinations exclude A1.
- I2(5)-to-root angles divide by `2b`; destinations exclude A1 and I2(5).
- H3-to-root angles divide by `3c`; destinations exclude all three source types.

The three families are disjoint. In particular, I2(5)-to-H3 pairs belong to
the second family. At fixed exact counts, the source multiplicities cancel
the denominators, and the vertex total remains additive over its components.
The separate review checks all 71 finite types, all 410 exact count triples,
all marked-root multiplicities and every saved maximizing forest. The
maximum is **21977/500 = 43.954**, with margin **23/500** below 44.

The unchanged short-path/H4 baseline and nonnegative bonuses retain all
triangle and persistent quadrilateral payments used by the arithmetic
attachment queues. The 76 triangle markings are checked again in the binder.

## Complete direct payments

The independent weight review covers all 2,151 optimization-input rows and
89,255 common spectator cases. This includes every one of the 1,113
dimension-feasible singular finite-queue stops, after accounting for all
2,070 stops. Expanded states are retained. All 2,397 original literal-input
maps and their exact dotted parameters are rechecked.

The coverage binder additionally traverses every earlier sealed branch,
singular first extension and continuation, dotted stopping state, and
complete finite-queue stopping state. It recalculates their payments with
the new rule, rather than inheriting previous numerical values. All 2,661
zero-touch rows remain explicitly represented; 1,634 fit dimension 22 and
receive **91,726** common external-count checks. The remaining 1,027 rows
are excluded at this target dimension by corner rank. The output includes
471 literal map checks and 79 recursive manifests with 1,023 hash references.

For a k-node core, residual spectator rank is `24−k`. The binder independently
reconstructs every nonnegative triple `(e,f,g)` with
`e+2f+3g≤24−k` whose leftover rank is zero or at least two. A2 and A3 remain
available to fill every leftover rank at least two. It requires one case
for each triple, no duplicates, four nonnegative corner weights, and an
exact sum at least one. The same external triple is added to all four
corners before their individual counts and denominators are evaluated.

The geometric completeness premises are the frozen minimal-union
classifications, complete one-white domains, all singular zero-touch
alternatives, exact dotted parameters, and justified dimension-29 queue
boundary. Conditional sealing makes marked components entire unless a
touching extension is followed. All orthogonal remaining components occur
in the spectator enumeration. No new arithmetic queue cutoff or numerical
distance sampling is introduced.

## Dimension-22 contradiction

Every triangle receives two, every quadrilateral one, and other weights
are nonnegative. Hence the average vertex weight is at least
`5 f₂/f₀−binom(22,2)`. The symmetric nonnegative h-vector has minimum paired
coefficient 55, with strictness supplied by `h₀=h₂₂=1`. Thus `f₂/f₀>55`,
and the average vertex weight is strictly greater than `275−231=44`.
The exact vertex upper bound 43.954 contradicts this.

The predecessor already excludes dimensions 23–29. Combining the new
exclusion with it and the published compact ceiling 29 gives the stated
Q(√5) bound 21. Other fields and noncompact cases are handled by the
[separate all-field composition](fresh_arithmetic_dimension21_composition_20261016.md).

```sh
python3 scripts/fresh_sqrt5_h3_hub_coverage_20261016.py
python3 scripts/fresh_sqrt5_h3_hub_coverage_20261016.py --verify
```
