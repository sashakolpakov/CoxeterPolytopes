# All-parameter gate for 280 specified transverse-edge cap families

**2026-10-14. No new finite-volume example is claimed.**

The [79-wall chamber](fresh_h20_e8_cusp_completion_20261014.md) has a
completed H12 face supported by E8 and 35 finite vertices on that face.
Its saved strict spacelike escape survives every compatible cap in the
280 families below. This explains why these simplest edge extensions
cannot finish the chamber. It does not exclude caps with several E8
contacts or other residual directions.

The [producer](../scripts/fresh_h20_e8_transverse_edges_20261014.py),
[exact certificate](fresh_h20_e8_transverse_edges_20261014.json), and
[independent interval replay](../scripts/fresh_h20_e8_transverse_edges_20261014_review.py)
with its [result](fresh_h20_e8_transverse_edges_20261014_review.json) are saved.
Both replays pass. The second reconstructs the weights and certificate
coefficients using separate sparse field multiplication and rational
radical intervals; it imports no producer code. The inherited chamber
and the completeness of its H12 vertex list remain certified dependencies.

Let V be one of the 35 future timelike vertices, let omega_j be an E8
fundamental weight, and put d_j=omega_j squared and q=-V squared.
The norm-two family

    r(c) = -c omega_j + sqrt((c² d_j - 2)/q) V

is orthogonal to the 19 old walls supporting the transverse edge obtained
by removing E8 node j at V. We take the positive square root. Its pairing
with the removed norm-two E8 root is -c. A Coxeter relation therefore
requires either c=1 or c>=sqrt(2): the latter includes all higher finite
labels, parallel walls, and every disjoint-wall distance.

Write u=sqrt((d_j-2/c²)/q), which increases with c. For the saved escape E,

    r(c).E / c = A - beta u,
    A = -omega_j.E > 0,   beta = -V.E > 0.

For an old wall N_i, write a=-omega_j.N_i and b=-V.N_i>=0. Acuteness
requires a-bu<=0. Exact inequalities divide all 35×8 families as follows:

| Families | Certificate |
| ---: | --- |
| 231 | Already A-beta u<0 at c=1; monotonicity covers every larger c. |
| 44 | One retained wall forces u>=a/b>=A/beta, or is nonacute for every parameter. |
| 5 | Every c>=sqrt(2) retains the escape. At c=1, four are nonacute; the fifth gives a cosine strictly between 0 and 1/2 against wall 11, hence no Coxeter angle. |

The exceptional fifth family is vertex 1, E8 position 4. Its forbidden
cosine is approximately 0.406879, but the certificate uses exact radical
inequalities. No finite label cutoff or numerical tolerance enters the
proof. The other four are vertices 17, 21, 25 and 29, E8 position 0.
All indices are zero-based in the linked source records.

Thus no finite collection of compatible caps from these families removes
the spacelike point E. If the resulting chamber has nonempty interior,
mixing E with an interior point gives a strict spacelike escape whenever
needed, proving infinite volume. General E8-coupled caps are outside this
result and remain the constructive search direction.

Reproduce with:

```sh
python3 scripts/fresh_h20_e8_transverse_edges_20261014.py --verify
python3 scripts/fresh_h20_e8_transverse_edges_20261014_review.py --verify
```
