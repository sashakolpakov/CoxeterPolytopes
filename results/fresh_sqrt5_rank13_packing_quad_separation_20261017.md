# Exact quadrilateral separation of the simplex-packing rule

The current generic simplex-packing rule fails quadrilateral coverage: 77,076 of the 81,183 determinant-flexible canonical full-star aggregates in this saved dotted domain are unpaid. Among these, 10,396 failures use only vertex forests already represented in the optimized model. The successful simplex certificate remains valid on its stated partial domain.

The [source](../scripts/fresh_sqrt5_rank13_packing_quad_separation_20261017.py) and [result](fresh_sqrt5_rank13_packing_quad_separation_20261017.json) bind the frozen [weight rule](fresh_sqrt5_rank13_simplex_packing_repair_20261017.md), [full-star aggregation theorem](fresh_quad_product_aggregation_20261017.md), [local dotted cores](fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.json), and [finite spectator contexts](fresh_sqrt5_rank13_dotted_b1_spectators_20261017.json). No optimizer or new face search is used.

| Exact census | Count |
| --- | ---: |
| Marked source cores / connected core isomorphism classes | 1,399 / 1,335 |
| Raw spectator contexts / canonical full stars | 119,698 / 108,805 |
| Raw flexible / fixed contexts | 89,329 / 30,369 |
| Canonical flexible / fixed full stars | 81,183 / 27,622 |
| Flexible aggregates paid / unpaid | 4,107 / 77,076 |
| Unpaid aggregates using only represented forests | 10,396 |
| Finite corner forests reached by the flexible domain | 6,521 |
| Zero-extended forest types encountered | 4,688 |

Fixed-determinant aggregates are counted but are not evaluated in this pass.

## Exact full-star and corner map

All 1,399 source cores are connected. Core isomorphisms preserve every finite bond label and both exact rational coefficients of every dotted squared product. Color refinement only filters candidates; each retained map checks all pairs, including nonedges. The complete positive-definite two-root deletion graph and every quadrilateral `K2,2` partition are checked under each map.

The connected core is the unique Lorentz component of the full star. Its orthogonal finite spectators are identified by their complete finite-component multiset. This gives a full-star isomorphism key and removes duplicate markings before testing payments.

If the core has `n` normals, its inertia is `(n-2,1,1)`. No principal block obtained by deleting only one root can be positive definite: eigenvalue interlacing leaves a nonpositive eigenvalue. Therefore both deleted roots of a full positive-definite vertex block must lie in the core. Adding orthogonal finite spectators creates no new deletion edges or quadrilaterals.

For a `K2,2` with sign sides `{i,a}` and `{j,b}`, the vertex omitting `i,j` has released corner ports `a,b`. The source reconstructs that entire finite vertex forest and its exact marked-root-pair orbit. Every quadrilateral contributes four occurrences exactly once. Its full-star demand is the number of distinct `K2,2` blocks, even when larger product faces overlap.

The new simplex-block bonuses pay only their actual simplex flags. They contribute **zero** to every quadrilateral aggregate here. The evaluated payments consist solely of the frozen rule's pair-orbit weights, with zero extension on absent forests.

## Determinant classification for a corank-one star

A connected `K`-valued root gauge with the same normalized Gram differs by one common totally positive scale and a diagonal congruence over `K`. For a singular core, choose any principal basis omitting a root with nonzero kernel coordinate. The nonsingular quotient has rank `n-1`, so its determinant changes by that common scale to the power `n-1`, times a square. The basis choice itself changes the determinant only by a square.

Consequently the core determinant squareclass is invariant when `n` is odd. Each finite spectator contributes its actual component rank. The conservative full ambient rank-14 test is therefore fixed when `n` is odd and every spectator component has even rank; otherwise an odd-rank connected factor leaves determinant compatibility flexible. This is a determinant test only, not a claim of matching Hasse invariants, a common lattice, or simultaneous global realization.

## A represented-forest zero-payment witness

Context 1527 of the frozen spectator catalogue uses source core 0 and spectator forest `A2 + A9`. The core is a compact quadrilateral with consecutive angles

`pi/2, pi/2, pi/3, pi/2`.

Its opposite squared normalized products are `t=(7+sqrt(5))/2` and `u=6+2sqrt(5)`. They satisfy `(t-4)(u-4)=4`; the physical products are strictly dotted and the conjugate products lie below 4. Its certified physical and conjugate inertias are `(2,1,1)` and `(3,0,1)`. The quotient rank is 3, so this context is determinant-flexible.

The four local corner maps, written as omitted roots followed by released ports, are

- `(1,0) -> (3,2)`;
- `(1,2) -> (3,0)`;
- `(3,0) -> (1,2)`;
- `(3,2) -> (1,0)`.

Three corners have forest `2 A1 + A2 + A9`, with the two isolated A1 roots marked. The other has forest `2 A2 + A9`, with both roots of one A2 component marked. Both forests already belong to the optimized model. The exact aggregate inequality is

`3 w(2A1+A2+A9; isolated pair) + w(2A2+A9; A2 edge) >= 1`.

Both weights are zero in the frozen rule. No simplex packing charge applies to this quadrilateral. This supplies a concrete represented-forest failure without any absent-type convention. The JSON also stores the absolutely smallest recorded failure with spectator `B11`, and the smallest represented-forest failure with spectator `A1+B10`; all three use the same four-root core.

## Saved data and scope

The result stores canonical core/corner maps, representative indices into the existing spectator catalogue, all unpaid context indices with exact payments and demands, the complete set of encountered zero-extended forests, and an exact trace digest for every evaluated flexible aggregate. It does not duplicate the original full spectator export.

Generation and a fresh saved replay both passed. These are individually certified local arithmetic quadrilateral configurations, not compact ambient polytopes. Their unpaid aggregates refute coverage by this particular rule. They do not refute the simplex incidence theorem, establish a dimension bound, or make the saved dotted domain a complete face catalogue.
