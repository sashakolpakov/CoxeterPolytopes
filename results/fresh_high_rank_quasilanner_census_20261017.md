# Exact reproduction of the order-six-to-ten quasi-Lannér census

The [exact source](../scripts/fresh_high_rank_quasilanner_census_20261017.py)
examines **146,262 attachments to 24 affine source diagrams**, retains 84
labelled attachments, and identifies **26 diagrams up to labelled graph
isomorphism**. The [saved census](fresh_high_rank_quasilanner_census_20261017.json)
records every retained attachment, its explicit isomorphism, every vertex
deletion and the relevant diagonal-decrease thresholds.

| Number of nodes | Nonisomorphic diagrams |
|---|---:|
| 6 | 12 |
| 7 | 3 |
| 8 | 4 |
| 9 | 4 |
| 10 | 3 |
| Total | 26 |

These are **classical diagrams, not new examples or a new classification**.
The result reproduces the higher-order part of
[M. Chein's 1969 classification, printed pp.15–16](https://www.numdam.org/item/M2AN_1969__3_3_3_0.pdf).
Chein's printed p.4 states the order-ten ceiling, attributing it to Vinberg,
and describes the finite/affine principal-subdiagram criterion. Those
classification facts are published mathematical inputs here; a search
restricted to ten nodes would not itself establish the ceiling. A
[repository copy of the primary PDF](../data/references/fresh_finite_edge_cusp_tradeoff_20261017_chein.pdf)
is also available through the separately frozen finite-edge package.

## Why the finite search covers these orders

Use the usual norm-two Gram convention: diagonal entries are 2, a missing
edge denotes label 2, label 3 has product −1, and label 4 has product
`-sqrt(2)`. A minimal indefinite diagram has indefinite full Gram matrix and
positive-semidefinite every proper principal submatrix. It is connected,
since otherwise a smaller connected component would already be indefinite.

For orders six through ten, the compact Lannér branch is absent: the
published compact classification has at most five nodes. Consequently at
least one proper principal submatrix is singular positive semidefinite and
contains a connected affine component A. Its kernel has strictly positive
coordinates. If another node is joined to A, the enlarged block cannot be
positive semidefinite: its affine kernel has a nonzero pairing with the
new node. Indeed a vector of zero quadratic form in a positive-semidefinite
matrix must lie in its full kernel, while the positive coefficients and
nonpositive Coxeter contacts make this pairing strictly negative.

Connectivity supplies such an adjacent outside node. Minimal indefiniteness
therefore forces A together with that node to be the whole diagram. Thus
every diagram under consideration is an irreducible affine diagram of
order `n-1` with one additional, nonorthogonal node. This also explains why
an arbitrary disconnected affine forest is not needed as a source.

There are no edge labels 5 or larger. Any edge in a connected diagram with
at least six nodes can be extended to a connected induced subdiagram on
five nodes, by adjoining adjacent vertices. This is a proper principal
subdiagram and must be finite or affine. The standard finite/affine
classification has no connected five-node diagram containing an edge
label at least 5. The smaller exceptions, such as the noncrystallographic
finite diagrams or the rank-two families, cannot survive this extension.

Dotted pairs are already indefinite two-node subdiagrams, so are excluded
immediately. A parallel pair is a connected affine two-node subdiagram;
adjoining any adjacent node makes a proper indefinite subdiagram by the
same positive-kernel argument. Parallel bonds are therefore also excluded
at these orders. The attachment alphabet is exactly `2,3,4`, including
orthogonality; there is no numerical cutoff on a hidden label or dotted
distance.

The complete affine source list has ranks four through eight:

- rank 4: affine A4, B4, C4, D4 and F4;
- rank 5: affine A5, B5, C5 and D5;
- ranks 6, 7 and 8: affine A, B, C, D and the corresponding E type.

These are the 24 source models. For a source on `n-1` nodes the program
tests all `3^(n-1)-1` nonzero attachment vectors. Their total is 146,262.

## Exact acceptance, signature and isomorphism checks

Every one-node deletion is decomposed into connected components. Each
component must be a recognized finite or affine Coxeter diagram. Finite
and affine trees are recognized by their complete labelled tree codes;
the affine A cycles are recognized separately. This supplies exact
positive semidefiniteness of every one-node deletion, and hence of every
proper principal submatrix. Rejection depends on this classification,
not on approximate eigenvalues.

For every newly encountered isomorphism type the source reconstructs the
full Gram matrix over exact rational numbers and `sqrt(2)` and proves its
determinant negative. Interlacing with the positive-semidefinite
one-node deletions then gives precisely one negative eigenvalue. The
nonzero determinant excludes a radical, so the full signature is
`(n-1,1,0)`.

Refinement is used only to select possible isomorphism matches; the actual
merges use a complete labelled graph-isomorphism test. All 84 retained
attachments have an explicit saved map to one of the 26 final diagrams.
The classical order counts are asserted only after enumeration, as a
comparison check; they are not used to prune or accept candidates.

The export distinguishes finite from affine node deletions. For a diagram
with a unique affine deletion j, it also records, for every other deleted
node b, the exact threshold

\[
 \tau_b=\frac{1}{(G_{\widehat b}^{-1})_{jj}}.
\]

All these thresholds are checked positive and rational. They describe
positive-semidefinite rank-one diagonal decreases at j. Applying them to
an actual simplex face with additional white supporting roots still
requires the appropriate geometric and Schur-complement argument.

## Scope and frozen validation

The finite/affine classifications, compact Lannér order bound and general
order-ten ceiling are explicit published inputs. The computation
exhausts the stated attachment domain and reproduces the classical
higher-order diagrams. It does not identify an induced geometric simplex
face with a Coxeter simplex, settle arbitrary extra supporting roots or
black/white partitions, or establish a new ambient dimension bound.

**FINAL/FROZEN:** generation 76418 and fresh saved replay 93878 passed;
the root agent confirmed the latter terminated with exit code zero.
The saved replay log ends with the complete `146262 / 84 / 26 / 24`
PASS result. No additional full replay was started to prepare this report.
Python syntax, all four direct source-manifest hashes and local report
links were checked. The frozen script and JSON were not edited.

| Artifact | SHA-256 |
|---|---|
| Source | `ab08dbff4912545541bba7117ce66848269e527e4f53270dcde138ef27812399` |
| Saved JSON | `6fe2aa3a6b5cd9ef70e25c2efea1d484f2f2c4a17219a9e1ea1966c8bc0e2729` |

The saved complete rejection-trace digest is
`04bbd355d85dc1f1ba411e7e179b12af09fb89fb5f6d5e66dcf4ab565a6c8ce3`.
This is a reproducible classical census and exact data package, not a
novelty claim or a formal proof-assistant verification.
