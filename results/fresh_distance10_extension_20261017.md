# Bounded supporting extensions for two distance-10 overlap cores

For the two overlap families below, every actual supporting-white extension leaves a preferred marked corner component of order at most **17**. Both families permit nonzero Euclidean residuals; their supporting diagrams need not be isolated unchanged. The proof controls arbitrarily many additional whites through a residual norm budget, not a finite search cutoff.

This is a local result under the assumption that no path of length at most 10 joins the two marked opposite pairs with white internal nodes. It is not a complete distance-10 classification, an ambient dimension bound or a finite-volume realization claim. The [exact checker](../scripts/fresh_distance10_extension_20261017.py) and [certificate](fresh_distance10_extension_20261017.json) use the literal L22 and L25 diagrams from the [frozen census](fresh_high_rank_quasilanner_census_20261017.json).

## Common geometric setup and residual budget

All normals have square 2. Let G be the Gram matrix of one of the stated nonsingular Lorentzian cores. The actual supporting white diagram W is positive definite; W plus any one black root is positive definite; W plus an adjacent black pair is positive semidefinite. These are the actual face-star hypotheses used in the [completed distance-11 proof](fresh_selective_distance11_dimension126_20261017.md).

For an additional supporting white e with core contacts b, orthogonal projection onto the core span has square `b^T G^-1 b`. The core contains a timelike vector, so its orthogonal complement is positive definite. Write

\[
 e=e_{\rm core}+w,\qquad w^2=2-b^TG^{-1}b\geq0.
\]

At a preferred corner, the nonnegative inverse of each finite component restricts a new nonzero contact to inverse diagonal at most 2. The magnitude of such a norm-two Coxeter contact is at least 1. The positive kernel of an affine component forces every extra contact there to vanish. The omitted black nodes are allowed separately. These tests give the explicit end domains below.

We will use the following elementary fact. Suppose at most one supporting white e contacts the core, its residual has square \(\rho\), and all remaining supporting whites are core-orthogonal. Their Gram matrix E is a positive-definite principal submatrix of W. If c is the contact vector of e with those roots, Euclidean projection gives

\[
 c^T E^{-1}c\leq\rho.
\]

Each contacted connected component costs at least 1/2: projecting onto any one contacted square-two root already costs at least 1/2. Hence \(\rho\leq2/3\) permits at most one contacted component.

If that connected component has at least three nodes, choose a connected three-node principal subdiagram through a contacted node. When the node is an end of a path with bond magnitudes a,b≥1 and possible extra contact c≥0, its inverse diagonal is

\[
 \frac{4-b^2}{8-2(a^2+b^2+c^2)-2abc}\geq\frac34.
\]

The denominator is positive. After clearing it, the inequality becomes
`6a²+2b²+6c²+6abc−8>=0`, immediate from a,b≥1. If the contacted node is the middle of the path, the corresponding difference is `6a²+6b²+2c²+6abc−8>=0`. Projection onto a larger positive-definite component cannot cost less than projection onto this three-node subspace. Thus a contacted component of rank at least three costs at least 3/4.

For a two-node connected component with bond magnitude t≥1, the endpoint inverse diagonal is `2/(4−t²)>=2/3`. Equality requires t=1. If both nodes are contacted, positivity of the inverse gives a larger cost; at the minimal simple diagram the cost is already 2. Therefore:

* if \(\rho=2/3\), the only nonempty attached component is A1 or A2, all contacts simple, and e meets only one end of A2;
* if \(\rho=1/2\), only one A1 with a simple contact is possible.

This bounds the **whole** attached component, including roots connected indirectly to e. Every other component is orthogonal to both e and the core. No rank cutoff on the ambient supporting diagram was used.

## Four-node reversed L25 overlap

Use the norm-two L25 edges

```
0-1, 0-2, 2-3, 0-4, 4-5, 5-6, 6-7, 7-8, 8-9,
```

all labelled 3. The second copy has map

```
(0,1,2,3,4,5,6,7,8,9) -> (10,11,12,13,14,15,9,8,7,6).
```

The resulting 16-root Gram has determinant -3 and signature (15,1). Deleting centers 0 and 10 leaves a positive-definite block and Schur matrix

\[
 \begin{pmatrix}-1/18&-1/9\\-1/9&-1/18\end{pmatrix}.
\]

Each black pair is one of {1,2}, {1,3}, {2,3} in its own copy. Choose preferred root 2, 1, 2 respectively. In all nine contexts the preferred marked component has order 14. Its finite inverse or affine kernel restricts contacts to the end sets {1,2,3} and {11,12,13}. A root contacting both ends would give a marked path of length at most 4. The exact checker also confirms that any unspecified direct cross edge between the exclusive core parts would give a path of length at most 10. The core Gram is consequently fixed under the hypothesis.

The full inverse block on either end is

\[
 \begin{pmatrix}
 2&2&1\\
 2&10/3&5/3\\
 1&5/3&4/3
 \end{pmatrix}.
\]

All entries are positive. The norm budget and the Coxeter contact alphabet therefore leave only two one-end types:

* **A:** a sole contact -1 at node 1 (or 11), projection square 2 and zero residual;
* **B:** a sole contact -1 at node 3 (or 13), projection square 4/3 and residual square 2/3.

Node 2 is impossible because its inverse diagonal exceeds 2. Two nonzero contacts at one end already cost more than 2. The next possible single-contact magnitude after 1 is sqrt(2), which also exceeds the budget at either remaining node. Some contexts eliminate type A already at the preferred-corner test; retaining it in the following analysis is conservative.

The cross and same-end inverse entries give complete control of multiple core-contacting whites. The A root is uniquely in the core span. Its product with a same-end B root is 1, contradicting acuteness; with an opposite-end B root it is -2, contradicting positive definiteness of W. Opposite A roots have product -3, also impossible in W. A second copy of the same A root cannot occur in a positive-definite supporting diagram. Thus an A root excludes every other core-contacting white. It is automatically orthogonal to all core-orthogonal whites and adds at most one node.

For two B roots at the same end, Cauchy in the positive residual space gives product at least `4/3−2/3=2/3`, violating acuteness. For B roots at opposite ends, the product is at most `−4/3+2/3=−2/3`, so there is a nonright edge between them. That edge gives a marked white path of length at most 5, forbidden by hypothesis. Hence there is at most one B root in total.

Apply the residual lemma with \(\rho=2/3\). Only an A1 or A2 tail can attach to that B root, with all bonds simple. At most three new connected roots occur: B plus at most two core-orthogonal whites. The marked component therefore has order at most `14+3=17`. Other supporting whites cannot connect to it.

The checker verifies all allowed local variants: six context-permitted A choices and 54 B choices with tail lengths 0, 1 or 2. Every supporting edge star remains positive definite and all four adjacent corners remain positive semidefinite. Six maximal-tail variants attain marked component order 17 at this local Gram level; no finite-volume chamber realization is inferred.

## Reversed L22/L25 terminal edge

Number L22 by 0 through 8 and map the L25 nodes by

```
(0,1,2,3,4,5,6,7,8,9) -> (9,10,11,12,13,14,15,16,8,7).
```

The 17-root Gram has determinant -2 and signature (16,1). Its two-center Schur matrix is

\[
 \begin{pmatrix}-5/36&-1/9\\-1/9&-1/18\end{pmatrix}.
\]

Use the same three black-pair choices and preferred roots. The nine contexts cover both orientations after exchanging the cores. Preferred marked components have order at most 15; a detached node 4 in some contexts is retained in the contact analysis. The allowable contact domains are contained in {1,3,4} on the L22 end and {10,11,12} on the L25 end. Contacts at both ends produce a marked path of length at most 5. Every extra direct cross edge gives a path of length at most 10.

The inverse block on the L22 end is

\[
 \begin{pmatrix}7/2&3&3/2\\3&4&2\\3/2&2&3/2\end{pmatrix},
\]

and the L25 end has strictly positive entries and diagonal minimum 4. Thus no L25-end contact is possible. On the L22 end the only possibility is a sole contact -1 at node 4: projection square 3/2, residual square 1/2. Every other node, multiple contact or larger single Coxeter contact exceeds the norm budget.

Two roots of this type would have product at least `3/2−1/2=1`, violating acuteness. There is therefore at most one. The residual lemma with \(\rho=1/2\) permits only one attached A1. At most two new connected whites occur, and the marked component has order at most `15+2=17`.

All 18 local variants, with tail length 0 or 1, pass the exact supporting-edge and adjacent-corner checks. Three attain marked component order 17. In contexts where node 4 is detached from the marked component, its extension does not enlarge that marked component.

## Actual ideal corners and validation scope

At an actual ideal corner, the [factor-fiber incidence lemma](fresh_cusp_factor_payment_framework_20261017.md) places the varying affine factor inside W plus the two side roots. The foregoing bounded connected component is therefore the component relevant to the actual face payment; additional cusp walls outside that star cannot enlarge it. The supporting-diagram argument uses W positive definite throughout and does not assume the entire ambient cusp is simple.

Exact generation passed: 18 preferred contexts and 78 permitted extension variants, with exact inverse blocks, signatures, residual endpoint Grams, all supporting-edge positive-definiteness checks and all adjacent-corner semidefinite checks. The general residual bound is a written proof using every indirect supporting root, not a claim extrapolated from those 78 variants. Generation and a fresh saved replay both passed. Seven direct source hashes, Python syntax and all local report links passed. The parent independently read the full mathematical argument and found no gap; that review did not duplicate the finite matrix calculations. Source, data and report are final/frozen.
