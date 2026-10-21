# Exact subdivision families in the unpaid frontier

**Nine unbounded necessary principal-star families, plus one bounded
family.** The [source](../scripts/fresh_sqrt5_subdivision_families_20261015.py)
and [certificate](fresh_sqrt5_subdivision_families_20261015.json) start from
ten literal 12-node graphs in the
[immutable batch610 export](fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.md).
Insert k white nodes on their common white edge 7–8, replacing that edge
by a path of k+1 edges, all labeled 3. Nine graphs retain four spherical
corners, physical inertia `(N-2,1,1)` and conjugate inertia `(N-1,0,1)` for
every integer k>=0. The remaining graph permits exactly k=0,...,10 under
these tests. Generation and saved replay pass.

This is an algebraic family theorem, not a finite-rank extrapolation. It
explains recurring classical components in the BFS and supplies precise
rank-dependent diagrams for prospective weights. It does **not** prove
that every subdivided graph is sealed, has no affine principal subdiagram,
is an actual future-oriented quadrilateral, or occurs in a compact
polytope. The active BFS remains partial and unchanged.

## Uniform Schur identity

Let H be the scaled Gram matrix of one base graph, and let the two edge
roots have common squared norm d. They have mutual product `-d/2`. Write
`z=e_7-e_8`. The inserted path block is `(d/2)Cartan(A_k)`; its first and
last inverse diagonal entries are `2k/[d(k+1)]`, and its first-to-last
entry is `2/[d(k+1)]`. Eliminating the inserted nodes changes the original
core to

`H - [d k/(2(k+1))] z z^T`.

Every base kernel has equal coordinates at 7 and 8, so this update preserves
its relation. Delete any coordinate on which that kernel is nonzero to
obtain an invertible grounded matrix J. Define

`q = (d/2) z^T J^-1 z`,

with z restricted to the grounded coordinates. The determinant multiplier
is `1 - [k/(k+1)] q`. If q<=1, that multiplier is strictly positive for
all finite k>=0. Along the continuous interval `0<=k/(k+1)<1`, no eigenvalue
can cross zero, so the grounded inertia stays unchanged. Adding back the
kernel and the positive-definite inserted path gives the asserted full
inertia. The same argument applies under the conjugate field embedding.

For each of the four corner matrices, the corresponding inverse cost is
also at most one in both embeddings. They are positive definite at k=0
and remain so for every k. Thus the corner assertion uses exact matrix
identities, not an assumption that a classical diagram pattern continues.
The checker independently reconstructs all base matrices and costs, and
directly verifies the full Schur identities at k=1,2,3 as normalization
checks; the Cartan inverse formulas prove the arbitrary-k statement.

| Base index | Grounded cost q | Allowed k from both embeddings |
|---|---|---|
| 343, 344, 367, 368, 395, 396 | 1 | Every k>=0 |
| 353, 377, 430 | `(67-sqrt(5))/76` | Every k>=0 |
| 333 | `(37-2sqrt(5))/38` | Exactly 0 through 10 |

For graph 333, the conjugate inequality is
`k < 1/(q'-1) = 2+4sqrt(5)`, which lies strictly between 10 and 11.
At k>=11 its conjugate has a negative direction. Its physical signature
and all four spherical corners alone would not detect this obstruction.

## Dimension-specific pattern table

The certificate includes 191 literal rows, with their complete corner
forests, marked-pair orbit keys and both external-A1 payments: k=0,...,19
for each of the nine unbounded families, and k=0,...,10 for graph 333.
Their corner ranks are `10+k`, covering the admissible ranks through 29.
All current payments are zero. This finite table is complete for those
parameter ranges; the displayed upper rank is a target-dimension table,
not a restriction in the all-k lemma.

These patterns can be used conservatively when designing weights. To turn
them into a new arbitrary-extension coverage theorem one must still prove
sealing, control extensions of their marked components, or run a separately
scoped exhaustive search at the target dimensions. No new dimension bound
is asserted here. The unchanged general BFS should not be called complete
merely because many consecutive ranks show the same pattern.

```sh
python3 scripts/fresh_sqrt5_subdivision_families_20261015.py --verify
```
