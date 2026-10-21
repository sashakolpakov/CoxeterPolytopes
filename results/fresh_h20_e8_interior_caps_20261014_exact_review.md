# Independent review of the 87-wall coupled E8 extension

**2026-10-14 — PASS.** The [independent checker](../scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py) reconstructs the extension from its frozen inputs using standard-library rational arithmetic. Its [certificate](fresh_h20_e8_interior_caps_20261014_exact_review.json) binds the [primary source](../scripts/fresh_h20_e8_interior_caps_20261014_exact.py), [87-wall data](fresh_h20_e8_interior_caps_20261014_exact.json) and inherited inputs. The primary saved replay also passes.

This checker imports no producer or project arithmetic. It implements recursive quadratic multiplication in `Q(√2)(√79)(√11)` and determines signs by exact relative-norm comparisons. This differs from the producer's sparse radical multiplication and outward rational interval signs. The independent square classes of 2, 79 and 11 ensure the eight-coordinate field encoding is injective.

The replay reconstructs the E8 inverse from the physical roots and both cap recipes. Their E8 profiles are `{3,6}` and `{4,5}`, with weight norms 66 and 68. The H12 residuals are positive combinations of the two specified source vertices, hence future timelike; the exact coefficients give norm two, the prescribed E8 contacts, the old-wall-19 contacts `−1` and zero, and mutual orthogonality. All four height sign choices are then reconstructed and compared with the eight saved normals. No assumption of global old-bank symmetry replaces the explicit pair checks.

The complete 87-by-87 Gram is recomputed. All 660 new unordered pair relations pass: 144 right, 20 at π/3, and 496 strictly dotted. The old 79 Coxeter relations remain an inherited dependency. All 87 norms and strict-center inequalities pass. Every saved facet point is independently reconstructed as the center's orthogonal projection onto its wall. Its norm and future orientation are correct, and its pairing with every other wall is strict. Thus all 87 listed inequalities are actual facets.

For each new cap, its complete finite-neighbor list has 21 walls. The independent checker repeats all 210 choices of 19 companions, for 1,680 tests in total, and checks the rational physical Gram by its own elimination. Exactly the saved eight distinct positive-definite rank-20 supports survive. Their physical vertex points are independently reconstructed and satisfy all 87 inequalities, with exactly the specified 20 equalities. Their types are four `A11+A7+A2` and four `A11+E7+A2`.

This is a complete finite-vertex census **incident to the eight new caps**: a finite Coxeter vertex has a spherical, hence simplicial, local chamber, so its 20 facets must appear in those finite-neighbor subset tests. It is not a complete ambient vertex census. The proof does not infer actuality merely from a positive-definite Gram; the physical projected points and all remaining strict inequalities are checked separately.

Four pairs of these vertices share a positive-definite rank-19 support. Their future timelike segments satisfy the shared equalities and all other inequalities strictly in the interior. They therefore give four actual finite edges. This verifies those four edges, not a complete ambient edge census. In particular, the cap with profile `{4,5}` already supplies a finite corner individually; the other cap completes its adjacent edge.

All eight new walls are strict on the 35 inherited H12 face vertices and on both previously completed cusp rays. Therefore those finite vertices and the two bounded cusp links survive. The full inherited cusp-link and H12 source completeness proofs are dependencies, not repeated here.

The original unflipped cap pair cuts the saved 79-wall escape. A different exact future spacelike point satisfies all 87 inequalities strictly, with squared norm

`1336028658157091827 / 250000000000000000`.

This proves the new chamber still has infinite volume. The result is a coherent finite-edge extension with eight additional finite corners, rather than a finite-volume completion. No unrestricted completion obstruction or complete ambient corner/cusp classification is claimed.

Reproduce the independent saved check with:

```sh
python3 scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py --verify
```
