# Independent review of the dense-third core obstruction

**PASS.** The [independent checker](../scripts/fresh_h9_star_third34_20261014_review.py) and [certificate](fresh_h9_star_third34_20261014_review.json) reconstruct the new domains and polynomial arithmetic of the [primary theorem](fresh_h9_star_third34_20261014.md). They use standalone `Fraction` arithmetic in \(\mathbb Q(\sqrt2)\) and sparse multivariate polynomials, with no producer, SymPy or predecessor-bank arithmetic imports.

The theorem applies to the literal star case 0. The four two-root blocks have finite labels at least 4, with the anchor block at least 7. The third leaf has own labels at least 3 and 4, in either order. Each anchor has exactly one own label 3 and one right contact; their nonright contacts can be shared or opposite. All noncore contacts remain free. This is an exclusion of that family, **not** of the entire incidence.

## Literal domain checks

The review verifies the source manifests, all 378 vertex incidences used as the geometry input, the nine maximal induced core faces, and the fact that every core pair meets. It independently reconstructs all 30 relaxed Gram matrices and their strictly positive null vectors:

- 24 affine 4–3–4 paths force the four blocks mutually right.
- Two affine 4–4 paths bound each third-anchor label by 3.
- Four affine all-3 triangles are used conditionally: if an anchor contacts an owner root nonorthogonally, its third contact must be right.

Every displayed support is contained in its saved actual vertex. The four projection budgets are checked literally, including both retained last-root choices. An anchor's retained cost is \(1/4\), and a foreign block has inverse diagonal at least 2 and nonnegative inverse entries. Thus an anchor has at most one foreign contact, necessarily label 3. The third leaf's retained cost is at least \(1/2\); any foreign contact would make the cost at least 1 and contradict strict positive definiteness.

No unused leaf contact is needed by these arguments. Both owner reversals and both last-root reversals are checked against the complete vertex and minimal-nonface lists, producing eight correctly stated families.

## Independent continuous arithmetic

Write \(d=1-c^2>0\), and let \(N_X,N_F,N_T\) be the numerators of the primary quantities \(X,F,T\). Instead of repeating the producer's rational-expression simplification, the review uses the two-dimensional determinant identity
\[
N_XN_F-N_T^2
=\frac{1-c^2}{4}(\alpha n_1-\beta n_0)^2.
\]
It independently derives the cleared polynomial
\[
\frac32(N_X-d)+\frac{(\alpha n_1-\beta n_0)^2}{4}
-N_F-2\delta_{\max}N_T-\delta_{\max}^2d
\]
for each of the six possible owner-contact count pairs. All saved unshifted expressions agree. Expanding with the binomial theorem after
\(\alpha=1/2+x\), \(\beta=\sqrt2/2+y\), \(c=\sqrt2/2+z\)
reproduces all 66 coefficients exactly: each polynomial has 11 strictly positive coefficients, including a positive constant.

A separate formal eight-variable expansion checks the determinant difference
\[
(X-1)(Y-7/2+F_M+2\mu)
+(\delta_{\max}-\delta)(\delta_{\max}+\delta+2T).
\]
Its first term is strictly positive because \(X>1\), and the last-block contribution satisfies \(Y>7/2\) for both shared and opposite single-contact anchors at every finite label at least 7. The second term is nonnegative by the conditional triangle bounds. Both residual norms are negative, and their determinant is positive; hence they would span a negative-definite plane in a Lorentz space of index one. This is a contradiction.

The finite-label parameter has no upper cutoff, and no arithmeticity or coefficient-field assumption is imposed on a putative realization.

## Verification and limits

Generation and the saved replay pass:

```sh
python scripts/fresh_h9_star_third34_20261014_review.py --verify
```

The certificate independently verifies the six new polynomials and every new literal domain witness. It inherits the saved incidence realization and does not repeat its primal construction. It neither runs a new label search nor claims a whole-incidence exclusion. The saved 20,837-visit assignment satisfies a verified family.

The prepared composition source was also read: combining the old exact \((3,3)\) third profile with the new alternatives “at least \((3,4)\)” and “at least \((4,3)\)” covers every pair of finite integer third labels at least 3. Transport through the twelve previously reviewed induced-core maps is legitimate when all new literal supports are checked, as that binder does. This last observation is a source/proof review; the composition's own generated certificate and search continuation are separate artifacts.
