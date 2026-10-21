# Independent review of the dense high-label gate

The [separate checker](../scripts/fresh_h9_dotted_central_20261014_gate_review.py) and [certificate](fresh_h9_dotted_central_20261014_gate_review.json) review the [primary local gate](fresh_h9_dotted_central_20261014_gate.json) without producer or project imports. They reconstruct all 444 case-5 vertices from the missing-face clutter and check 44 injective placements, 308 actual-face premises and all 24 deduplicated condition sets. This is a local obstruction, not a complete exclusion of case 5.

Normalize all roots to square one. Let A=I2(a), B=I2(b), with a≥6 and b≥7, and let x have nonright contact with both roots of A. The leaf y has a nonempty contact with B: otherwise its missing own triple would have positive-definite Gram matrix and could not be a missing face of a compact Coxeter polytope.

The actual positive-definite faces A+B, B+x and A+y force all cross-block and foreign-leaf contacts to be right. Indeed, any nonright contact to a block of label at least six yields a three-node principal path with determinant at most 1−3/4−1/4=0; an additional edge only worsens positive definiteness. These deductions use literal actual faces, not an assumed symmetry.

Choose a nonright own contact of y, and either own contact of x. The corresponding actual four-root face is the path Ai–x–y–Bj. Its determinant is

```text
(1−alpha²)(1−beta²)−delta² ≤ 9/16−delta²,
```

where alpha,beta≥1/2 and delta is the x–y cosine. Since cos²(π/5)>9/16, this positive face forces the x–y label to be at most four, including the right case. Hence delta²≤1/2. All four possible path premises are checked in each literal placement.

Project x and y off the positive rank-four span A+B. The foreign products vanish, so the residual off-diagonal is still −delta. For u=cos(π/a)≥√3/2, the dense own contacts of x give projection cost at least

```text
(1/2)/(1−u) ≥ 2+√3.
```

For y, a single own contact at least 1/2 gives cost at least 1/[4(1−v²)], where v=cos(π/b). The strict bound cos²(π/7)>4/5 makes this cost strictly greater than 5/4. The checker independently expands

```text
T7(z)−1 = (z−1)(8z³+4z²−4z−1)².
```

The cubic changes sign from −29/125 at 3/5 to 7 at 1. Its largest root is cos(2π/7), so cos²(π/7)>4/5. Monotonicity covers every integer b≥7, with no label cutoff.

Both residual diagonals are negative, with magnitudes at least 1+√3 and strictly greater than 1/4. Their determinant is therefore strictly greater than

```text
(1+√3)/4 − 1/2 = (√3−1)/4 > 0.
```

The residual two-plane is negative definite, contradicting Lorentz index one. All other labels and dotted distances remain free subject to the checked hypotheses. Fraction arithmetic and an independent quadratic-pair calculation verify the displayed thresholds and identities; the preceding monotonicity and Schur-complement argument explains their all-parameter scope.

Generation and saved replay pass. No traversal result follows from this review alone.

```sh
python scripts/fresh_h9_dotted_central_20261014_gate_review.py --verify
```
