# Independent review of the uniform-chain cap obstruction

The [uniform-chain obstruction](fresh_h20_uniform_chain_cap_gate_20261016.md) passes this independent arithmetic and proof review. The [review checker](../scripts/fresh_h20_uniform_chain_cap_gate_20261016_review.py) reproduces the complete 54,744-profile cost census and independently certifies the two continuous physical blockers. No gap was found in the all-label reduction or its stated construction-family scope.

The [review certificate](fresh_h20_uniform_chain_cap_gate_20261016_review.json) binds the primary PY/JSON, the literal physical inputs, and the frozen [independent 133-wall review](fresh_h20_s61_multicap_20261014_exact_review.md). It does not re-enumerate the whole 133-wall chamber.

## Independent arithmetic and coverage

The support Gram matrix is reconstructed from the 21 physical coordinates, using the separate radical-tower arithmetic of the frozen [bridge-edge review](fresh_h20_bridge_edge_closure_20261014_review.md). This is a different implementation from the uniform-gate producer's field helper. The new review does not import the uniform-gate producer or its interval class.

A reversed-coordinate rational Gauss–Jordan inverse and reversed LDL elimination verify positivity and the full inverse. Literal graph traversal gives component ranks 10,4,2,1,1, with only simple bonds and path components. The inverse is entrywise nonnegative. Only walls 17 and 31 have diagonal 1/2, and both are isolated A1 components; the next smallest diagonal is 2/3, and the four smallest sum to 7/3.

The cost census uses a squarefree-radicand dictionary in Q(√2,√3,√5), with multiplication by gcd reduction. This differs from the producer's tensor-basis multiplication. It traverses supports and labels in reverse order and evaluates diagonal terms plus unordered cross terms. Every sign uses outward rational square-root bounds. It reproduces:

| Cost | Profiles |
|---|---:|
| below 2 | 97 |
| exactly 2 | 8 |
| at least 137/66 | 54,639 |
| total | 54,744 |

After canonical serialization, the independently calculated full-profile digest also agrees with the primary. The eight minimizing profiles all have label 3 and choose one wall from each of `{1,10}`, `{17,31}`, and `{61,87}`.

## All finite labels and nonfinite support contacts

The analytic reduction is exhaustive, not a finite-label cutoff:

- Every nonzero finite contact has magnitude at least 1. Four contacts already cost at least 7/3, above 141/70.
- For labels at least 7, the cubic for `2cos(2π/7)` gives squared contact greater than 81/25. The review checks its negative value at 31/25 and the required strict lower bounds. A high label on a nonisolated node costs more than 54/25. On an isolated A1 together with another contact it costs more than 53/25. A lone finite A1 contact costs less than 2.
- Thus one, two or three finite contacts with labels 3–6 are the entire remaining finite case. The zero-contact profile has cost zero and cannot satisfy the norm equation.
- A nonfinite support contact has magnitude at least 2. It costs at least 8/3 away from the isolated A1 nodes, or at least 5/2 if another nonzero contact accompanies an isolated A1 contact. Therefore only a single contact at 17 or 31 remains. The norm equation makes that contact strictly dotted.

The verified gap `(2,137/66)` contains the whole required interval `(2,141/70)`. All strict inequalities and the analytic high-label bounds are checked with rational arithmetic.

## Geometry and the whole-interval blockers

The frozen chain consists of 26 norm-two roots in the three-dimensional orthogonal complement of the positive-definite support. Three consecutive spatial vectors are noncollinear. The review separately checks a positive rational interval for `2 sin(δ)(1−cos(δ))` using the inherited isolated parameter interval. Equal finite angles therefore force the residual direction to be `t x`, where `x²=−1`.

For a common nonright finite angle, `1/a ≤ t < 2/a`, with `a²>280`. The resulting support cost is `2+t²<141/70`. For a common right angle, the proposed normal is a nonzero negative combination of the support normals with nonnegative coefficients. It cannot have negative pairing with any common strict interior point. This handles the right-angle case without a hidden choice of the old center.

The two remaining axes are

`r_j(t)=t x−√(4+2t²) N_j/2`, for `j=17,31`.

The review reconstructs X and its norm directly from the earlier physical edge recipe, checks orthogonality to every support root, and normalizes `x=X/√(−X²)`. Each blocking root is positively normalized to square 2; this matters because not all stored old roots have norm 2. Its own Fraction interval arithmetic, with 320-bit outward square-root bounds, proves positivity throughout the closed interval `[1/a_hi,2/a_lo]`:

| Axis | Blocking old wall | Outward decimal enclosure |
|---|---:|---:|
| 17 | 83 | [0.88227, 1.15264] |
| 31 | 11 | [0.95007, 1.48126] |

The certificate stores exact rational endpoints. Both intervals are strictly positive, so each axis contradicts acuteness. This is a continuous calculation, not sampling of integer angles. The closed interval also covers the uniform parallel endpoint `t=2/a`; the same proof excludes that limiting uniform case.

## Scope and independence limits

The conclusion concerns a coherent acute extension meeting all 26 chain walls at one common finite Coxeter angle, with arbitrary real coefficients and arbitrary Coxeter contacts to the 18 support walls. The parallel endpoint observation uses the same closed-interval proof. Unequal chain angles, general nonfinite chain contacts, changed retained walls, and arbitrary H20 completions are not excluded by this review.

The new census arithmetic, rational inverse, graph reconstruction, interval combination and proof audit are independent. The physical radical engine is reused from an earlier independent review, explicitly source-bound; it is not newly reimplemented here. The previously reviewed chain parameter isolation and full 133-wall construction are inherited inputs. No independent full chamber census or new finite-volume claim is made.

```sh
python3 scripts/fresh_h20_uniform_chain_cap_gate_20261016_review.py
python3 scripts/fresh_h20_uniform_chain_cap_gate_20261016_review.py --verify
```

Generation and fresh saved replay both pass. The review source, certificate and report are final and frozen; the primary files were not changed.
