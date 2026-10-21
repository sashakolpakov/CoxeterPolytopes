# Supporting-white extensions of the three-node L25 overlap

This is an exact local lemma for one remaining distance-11 configuration. It does not establish a complete distance-11 payment rule or a new dimension bound. All roots below have square 2 and all diagrams are principal diagrams of an actual acute Coxeter face star.

Let two copies of the frozen census diagram L25 share their three terminal nodes in reverse order. Its edges are

```
0-1, 0-2, 2-3, 0-4, 4-5, 5-6, 6-7, 7-8, 8-9,
```

all labelled 3. The second copy has the node map

```
(0,1,2,3,4,5,6,7,8,9) -> (10,11,12,13,14,15,16,9,8,7).
```

Assume the two pairs of opposite black roots have no connecting path of length at most 11 whose internal nodes are white. Then every additional supporting white root is orthogonal to this 17-root core, with at most one exception. The exceptional root has contact -1 with node 3 alone, or with node 13 alone. In either case it lies in the core span. A preferred adjacent corner therefore has an isolated marked component of order at most **16**, even after all supporting whites are included.

## Exact core and the nine marked contexts

The [checker](../scripts/fresh_distance11_extension_review_20261017.py) reads only the literal L25 diagram from the [frozen exact census](fresh_high_rank_quasilanner_census_20261017.json); it imports no producer arithmetic. Deleting nodes 0 and 10 gives a positive definite 15-root block. Its Schur complement is

\[
 \begin{pmatrix}-1/15&-1/10\\-1/10&-1/15\end{pmatrix}.
\]

Thus the full Gram matrix \(G\) is nonsingular of signature \((16,1)\), with determinant \(-2\).

The checker examines all \(\binom72^2=441\) choices of marked pairs outside the common path. Exactly nine have no marked path of length at most 11: each pair is one of \(\{1,2\},\{1,3\},\{2,3\}\), with the second pair translated by 10. All four adjacent corner diagrams in each case are positive semidefinite. An unspecified extra edge between the two exclusive core parts always creates a path of length at most 11, so the displayed core has no unexamined free cross contacts under this hypothesis.

Choose the preferred marked node 2 for the pair \(\{1,2\}\), node 1 for \(\{1,3\}\), and node 2 for \(\{2,3\}\); use the translated choice on the right. Delete the other black node on each side. The resulting preferred corner has a connected marked component of order 15 and marked distance 12.

For each finite component of this corner, the checker computes its inverse diagonal. A new square-two root cannot have a nonzero Coxeter contact at a node whose inverse diagonal exceeds 2: the squared projection already exceeds its norm, since the inverse of a spherical Coxeter Gram matrix is entrywise nonnegative. For an affine component, its positive kernel forces every contact of a further root to vanish whenever the extended corner is positive semidefinite. After also allowing contacts at the two omitted black roots, the possible contact domain is contained in

\[
 E_L=\{1,2,3\},\qquad E_R=\{11,12,13\}.
\]

Every allowed node is within one edge of its own marked pair. A root contacting both ends would therefore produce a marked white path of length at most 4. Every nonorthogonal supporting white is consequently supported at just one end.

## The only two one-end directions

On either end the corresponding block of \(G^{-1}\) is

\[
 J_E=\begin{pmatrix}7/2&4&2\\4&6&3\\2&3&2\end{pmatrix}.
\]

Let \(b\) be the nonpositive vector of native contacts of a proposed root with the core. Every nonzero magnitude is at least 1. The core span contains a timelike vector, so its orthogonal complement is positive definite. Projection therefore gives the necessary inequality

\[
 b^T G^{-1}b\leq2.
\]

All entries of \(J_E\) are strictly positive. Its first two diagonal entries exceed 2, and its last is 2. Hence the only nonzero one-end possibilities are exactly \(b=-e_3\) and \(b=-e_{13}\). They attain equality, so their orthogonal residual is zero: the actual roots are uniquely

\[
 r=-G^{-1}e_3,\qquad s=-G^{-1}e_{13}.
\]

Their coordinates in the ordered core basis are

```
r = (-4,-2,-3,-2,-3,-2,-1,0,1,2,6,3,4,2,5,4,3),
s = (6,3,4,2,5,4,3,2,1,0,-4,-2,-3,-2,-3,-2,-1).
```

Direct exact multiplication gives \(r^2=s^2=2\), the prescribed sole contacts, and

\[
 (r,s)=-2.
\]

Both cannot be supporting roots of the quadrilateral: the supporting diagram \(W\) is positive definite, whereas their two-root Gram matrix is singular. A second root in the same direction is also impossible for distinct normalized facet normals. Every remaining supporting white that is orthogonal to the core is automatically orthogonal to the exceptional root, because that root lies in the core span. Thus there is at most one extra node in the marked component, and no indirect extension through other whites.

## Corner and orientation checks

For all nine contexts and either one-root choice, the checker verifies positive definiteness of the supporting white diagram and positive semidefiniteness of all four adjacent corners. Among the 18 preferred-corner variants the marked component is finite of order 16 in six cases, finite of order 15 in four, and affine of order 15 in eight. An extra node disconnected from that component does not enlarge it.

These are feasible local acute contact sets. For example \(x=-G^{-1}\mathbf1\) has square \(-1001/2\), contact -1 with every core root and contact -13 with each of \(r,s\). This establishes coherent orientation for the local calculation; it is not a finite-volume realization claim.

For an actual ideal corner, the [factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md) identifies its varying affine factor inside \(W\) plus the two side roots. Additional incident walls outside that corner star do not enlarge the marked component relevant to the face payment. This is the same geometric scope used in the [frozen distance-12 proof](fresh_selective_distance12_dimension136_20261017.md) and its [independent review](fresh_selective_distance12_dimension136_20261017_review.md).

## Verification scope

Exact generation and a fresh saved replay both passed. The [saved certificate](fresh_distance11_extension_review_20261017.json) records all nine marked contexts, exact preferred-corner contact restrictions, both unique roots, and all 18 extension variants. The arithmetic uses only Python integers and `Fraction`. The six direct source hashes, Python syntax and all local report links passed. The parent independently read the mathematical argument and found no gap; that review did not duplicate the exact corner census. No claim is made here about the other distance-11 core overlaps or an ambient dimension bound.
