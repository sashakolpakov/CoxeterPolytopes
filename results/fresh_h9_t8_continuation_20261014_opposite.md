# Opposite-root metric family in T-shaped case 8

The saved opposite-root mandatory assignment has no compact hyperbolic Coxeter realization in dimension nine. This remains true when both high block labels range over **every finite integer at least seven** and the two own nonright anchor labels range over **every finite integer at least three**. All other saved mandatory contacts stay fixed. This does not exclude the complete incidence type.

The [producer](../scripts/fresh_h9_t8_continuation_20261014_opposite.py) and [exact certificate](fresh_h9_t8_continuation_20261014_opposite.json) bind the [realized incidence](fresh_h9_tree_cube_wedges_20261014.json) and [933-node continuation state](fresh_h9_t8_continuation_20261014_case_08.json). The argument uses no choices for dotted distances and allows every omitted meeting contact that survives necessary actual-vertex constraints.

## Literal contact domains

Use unit-diagonal Gram matrices. The positive support is

```
S = [1,2, 4,5, 7,8, 12,13].
```

Its four blocks are respectively `I2(m), B2, B2, I2(n)`, with finite `m,n >= 7`. The anchors are 11 and 17. Anchor 11 is right to 12 and has own contact magnitude `alpha` at 13; anchor 17 is right to 13 and has magnitude `beta` at 12. Here `alpha=cos(pi/a)`, `beta=cos(pi/b)`, for arbitrary finite integers `a,b>=3`. The third normal, 14, has magnitude `1/2` at roots 2,4,5 and zero at the other support roots.

The certificate verifies all containing actual vertices, rather than assuming that the local normals meet:

- All 24 cross-block contacts are right. A nonright edge attached to the first high block gives a nonpositive three-path minor. Between the other blocks, the relaxed `4–3–4` chain has an explicit positive null vector. Restoring larger labels cannot make it positive definite.
- For anchor 11 retain root 13 and remove root 12; for anchor 17 retain root 12 and remove root 13. These actual supports give own cost at least `1/4`. Any first-high-block nonright contact alone costs more than one. In the two B2 blocks, at most one foreign contact is possible, and it must have label 3.
- If there is a foreign contact, the own contact must also have label 3: own label at least 4 costs at least `1/2`, and the foreign contact costs at least `1/2`. Consequently a variable `alpha>1/2` or `beta>1/2` occurs only for a foreign-free anchor.
- The two seven-support tests for normal 14 have fixed cost `1/4+1/4`. Every additional nonright contact costs at least `1/2`, so its stated projection is forced.
- The actual path `[11,12,17]` excludes anchor mutual label at least 6; labels 2,3,4,5 remain. Two newly computed actual-vertex LDL failures exclude labels at least 4 between normal 14 and either anchor. These latter contacts are therefore labels 2 or 3.

Failure under a contact increase follows from monotonicity of a symmetric nonpositive-off-diagonal Gram matrix: on a connected block, a least-eigenvalue eigenvector can be chosen nonnegative. These are compact-vertex tests; no noncompact realization claim is made.

## Exact continuous determinant gate

Put `t=cos(pi/m)^2`, `c=cos(pi/n)` and `B=1-c^2`. The positive support has determinant `(1-t)B/4>0`. Its rank-eight Schur complement would have rank at most two in an H9 realization.

Let `u,v` be the foreign projection costs (zero or `1/2`), and let `p` be the mutual foreign projection (zero, `1/2`, or `sqrt(2)/4`). Write `d=delta+p`, where `delta` is the mutual anchor cosine. Let `ac,bc` denote the pairings of the two residual anchors with the third residual normal; they are nonpositive. Specifically each is the negative of a label-2/3 cosine, minus `1/2+sqrt(2)/4` when its foreign contact is at root 4 or 5. The third residual square is `-Q`, where

```
Q = sqrt(2)/2 + 1/[4(1-t)].
```

For the two-anchor determinant `A/B`, its numerator is

```
A = [(1-u)(1-v)-d^2]B
    - alpha^2(1-v) - beta^2(1-u)
    - 2 d c alpha beta + alpha^2 beta^2.
```

For the three-residual determinant `N/B`, its numerator is

```
N = -Q A - B[(1-u)bc^2 + (1-v)ac^2 + 2d ac bc]
    + alpha^2 bc^2 + beta^2 ac^2 - 2c alpha beta ac bc.
```

The script verifies these identities symbolically and, independently within the producer, checks the fixed-contact determinant by polynomial convolution.

The second-kind Chebyshev identity gives

```
64x^3 - 80x^2 + 24x - 1 = 0,  x=cos(pi/7)^2.
```

The polynomial is increasing above `3/4`, and its value at `81/100` is `-559/15625`. Thus `t>81/100`, `c>9/10`. At `t=81/100`, write `Q0=sqrt(2)/2+25/19>2`.

For `alpha=beta=1/2`, the five foreign choices for each anchor, four mutual labels, and two choices for each third contact give 400 cases. The certificate reconstructs the two quadratic polynomials `A(c),N(c)` over `Q(sqrt(2),sqrt(5))` and checks **all three Bernstein coefficients** on the closed interval `[9/10,1]`: those of `A` are strictly negative and those of `N` strictly positive. This is an exact continuous interval proof, not a finite-label cutoff. The endpoint `c=1` is used only after multiplying the positive denominator; actual finite labels have `c<1`.

## All larger own-contact labels

A free own magnitude can vary only on a foreign-free anchor. For such an `alpha`, `u=p=0` and `d=delta`. The checked derivative identities are

```
A_alpha = -2 alpha(1-v-beta^2) - 2 delta c beta,
N_alpha/2 = alpha[Q(1-v-beta^2)+bc^2]
           + c beta[Q delta-ac bc].
```

The first expression is strictly negative: the other anchor either has no foreign contact and `beta<1`, or has foreign contact with `v=1/2,beta=1/2`.

For the residual determinant there are three exhaustive cases.

1. One anchor is free and the other has a foreign contact. Here `beta=1/2,v=1/2` and `|ac|<=1/2`. Setting `z=|bc|`, dropping the nonnegative `delta` contribution gives

   `N_alpha/2 >= Q/8+z^2/2-z/4 = Q/8-1/32+(z-1/4)^2/2 > 7/32`.

2. Both anchors have no foreign contact and the mutual label is nonright. Then `delta>=1/2`, `ac bc<=1/4` and `Q>2`; the second bracket is positive, as is the first. Both derivatives are positive.

3. Both anchors have no foreign contact and the mutual label is 2. In this case

   `N = Q[c^2-(1-alpha^2)(1-beta^2)] + (alpha bc-beta ac)^2 + 2alpha beta(1-c)ac bc - B(ac^2+bc^2)`.

   Since `alpha,beta>=1/2`, `c^2>81/100`, and `ac^2+bc^2<=1/2`, it follows directly that

   `N > 2(81/100-9/16)-19/200 = 2/5`.

Thus either decreasing the free own contacts reduces to the exact 400-case table, or the direct positive bound applies. The two-anchor determinant remains negative as the free contacts increase. Finally the full three-residual determinant increases with `t`, because its derivative is `-A/[4B(1-t)^2]>0`. The positive determinant at the lower `t` bound therefore excludes the entire two-high-label family. The full eleven-normal Gram has rank eleven, contradicting the ambient rank bound ten.

## Replay and scope

```
python scripts/fresh_h9_t8_continuation_20261014_opposite.py --verify
```

Generation and saved replay pass. The producer uses the frozen exact quadratic arithmetic for finite-vertex LDL guards and a separate SymPy algebraic field for the 400 polynomial tables. The [independent review](fresh_h9_t8_continuation_20261014_opposite_review.md) uses separate Fraction arithmetic for all 400 determinant/Bernstein cases and literal vertex guards, plus a separate rational-polynomial check of the general determinant and derivative identities. Its proof review checks the unrestricted own-label bounds. This result excludes one explicitly specified continuous mandatory-contact family, not the whole T-shaped incidence or all compact H9 polytopes.
