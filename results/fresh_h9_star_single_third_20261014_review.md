# Independent review of the single-contact third-leaf theorem

**PASS.** The [independent checker](../scripts/fresh_h9_star_single_third_20261014_review.py) and [certificate](fresh_h9_star_single_third_20261014_review.json) verify every new literal domain guard and independently reconstruct all eleven determinant polynomials and their 78 strictly positive coefficients. No producer or SymPy arithmetic is imported. The field operations use the previously independent `Fraction` basis \(1,\sqrt2,\sqrt5,\sqrt{10}\); the kernel checks separately use rational quadratic arithmetic.

The [primary theorem](fresh_h9_star_single_third_20261014.json) excludes a specified family on star case 0: the owner block has a finite label at least 4, the other three blocks have finite labels at least 7, the third leaf has own contacts \((2,k)\) with \(k\ge5\), and the two anchors share one own label 3 and one right contact. The third's owner roots and the shared last roots may be reversed. All noncore contacts are free. This is **not** a whole-incidence exclusion.

## Literal guards and complete domains

The review checks the 378 saved vertex incidences as input, the 23 minimal nonfaces, the nine maximal induced core faces and all 55 meeting core pairs. It then independently reconstructs the following 48 positive-kernel witnesses:

| Role | Count |
|---|---:|
| Cross-block right-angle guards | 24 |
| High middle blocks force anchors right | 8 |
| An owner contact forces owner label at most 5 | 4 |
| Third-anchor labels are at most 3 | 2 |
| A contact at the third's nonright root forces that third-anchor pair right | 2 |
| Shared anchors are mutually right | 1 |
| At most one owner contact across both anchors | 4 |
| The third cannot contact both anchors nonorthogonally | 1 |
| A contact at the third's nonright root forces the other third-anchor pair right | 2 |

For every witness, its support is contained in its displayed actual vertex. More than support containment is checked: every retained relaxed edge must follow from an explicit fixed label, lower bound, or the particular contact being tested. Thus the conditional cycles do not silently assume unknown contacts. The exact relaxed Gram, positive kernel and zero matrix product are independently reconstructed.

All three eight-root projection budgets are checked literally. An anchor retains own cost \(1/4\), allowing at most one foreign contact, necessarily label 3. The third retains cost at least \(\cos^2(\pi/5)>1/2\); any foreign contact would add at least \(1/2\), contradicting positive definiteness. The high-block paths and conditional cycles then give the complete seven-pattern residual domain.

With no foreign owner contact, the third-anchor cosine pairs are \((0,0),(1/2,0),(0,1/2)\). If an anchor contacts the third's right owner root, the same three pairs remain. If it contacts the nonright owner root, both third-anchor contacts must be right. Both choices of the anchor carrying the foreign contact have literal guards. Renaming these two residual vectors in the formula does not assume a global incidence automorphism.

## Unbounded tail and independent polynomial reconstruction

The review reconstructs \(T_7\) by the rational Chebyshev recurrence and checks
\[
T_7(u)-1=(u-1)(8u^3+4u^2-4u-1)^2.
\]
For \(u=\cos(2\pi/7)>1/2\), the cubic vanishes and is strictly increasing on \([1/2,1]\). Its value at \(31/50\) is strictly negative. Hence every finite last label \(m\ge7\) satisfies
\[
D=\frac1{1-\cos^2(\pi/m)}>\frac{100}{19}.
\]
No numerical cosine estimate or finite upper cutoff is used.

For each residual pattern, the review multiplies all three rows of the physical residual Gram by \(d=1-c^2\), expands its determinant over all six permutations, and divides the resulting polynomial by \(d^2\). This independently produces \(d\det G\); it does not use the producer's symbolic determinant simplification. Every saved unshifted expression agrees exactly.

The three no-foreign cases permit the full continuous owner range \(c\ge\sqrt2/2\), with \(c<1\). A foreign contact forces the integer owner label to be 4 or 5, giving eight further cases. The independent binomial expansions after the specified shifts reproduce ten positive coefficients in each no-foreign case and six in each foreign case: **eleven polynomials, 78 coefficients**, all strictly positive, including their constants.

The positive support has rank eight in a ten-dimensional Lorentz space, so its orthogonal complement has dimension two. The three projected residual vectors must therefore have determinant zero. The certified polynomial positivity, together with \(d>0\), makes that determinant strictly positive, yielding the contradiction. There is no arithmeticity or coefficient-field assumption on a putative realization.

## Verification and transfer scope

Generation and saved replay pass:

```sh
python scripts/fresh_h9_star_single_third_20261014_review.py --verify
```

The four source orientations are checked on the full vertex and minimal-nonface lists. The saved 26,241-visit assignment satisfies one of them. The review inherits the incidence realization; it does not reconstruct its primal polytope or rerun any predecessor determinant bank.

The prepared transfer binder was also source-reviewed. It requires the frozen independent theorem and twelve-core-map reviews, rechecks the complete induced face complexes, and verifies all 612 new path, cycle and budget support occurrences. Transporting the four orientations gives 48 additional predicates. This is valid local induced-core transport, not a claim of full-nerve symmetry. The binder's generation and any subsequent search are separate artifacts.
