# Disjoint-cap extensions, exact endpoint obstructions and star symmetry

No new finite-volume Coxeter polytope or new dimension is established. The constructive advance is an exact nineteen-facet H9 chamber with 53 finite vertices and 175 finite edges. Four compatible walls close six old open edges, but a strict spacelike escape still proves infinite volume. This checkpoint preserves the [previous frozen archive](unbounded_star_labels_and_uniform_h8_exclusion_2026_10_14.md).

## Exact edge construction and its limits

The [edge-walking report](fresh_h9_h8_disjoint_walk_20261014.md) starts from the earlier fifteen-facet chamber with 44 finite vertices and 167 finite edges. Its exact frontier has 62 unfinished finite-corner germs, 1,076 spherical extension profiles and 64 symbolic dihedral tails. A separate complete local census finds fourteen profiles giving a bounded ideal endpoint after one additional wall.

Four coherent caps can be written `(p_k,−√K)`, where `K=6+16√5/5` and every `p_k` lies in `Q(√5)`. Their norm is two. All 66 new wall pairs are Coxeter-compatible. They produce nine new finite vertices and close six old germs. The [independent review](fresh_h9_h8_disjoint_walk_20261014_review.md) reconstructs all 171 pairs, nineteen actual facets, all 53 finite vertices, all 175 finite edges and the six positive forward parameters in a separate number-field implementation. Its finite-vertex census filters all 92,378 nine-wall subsets independently of the producer.

The same exact escape survives the quartet. The resulting 127 open finite-corner germs consist of 56 surviving old germs and 71 at new vertices. The [new exact frontier](fresh_h9_h8_disjoint_walk_20261014_next_frontier.json) has 1,065 spherical profiles and fifty symbolic dihedral tails on those 71 directions. These are finite-corner frontiers, not censuses of edges with no finite endpoint.

The first bounded numerical scout found 23,926 actual-endpoint proposals and no escape cutter. A fourteen-profile affine scout likewise found no cutter among its 218 accepted events. The [second scout](fresh_h9_h8_disjoint_walk_20261014_next.md) found 21,024 accepted proposals on the 71 new germs and no cutter. Those bounded tests alone prove no universal exclusion.

The [all-parameter finite-endpoint theorem](fresh_h9_h8_disjoint_walk_20261014_universal_gate.md) strengthens the first result: all 1,076 spherical profiles and all 64 continuous dihedral tails on the original 62 germs miss the literal escape, for every endpoint parameter and every integer Coxeter label. A closed, outward-rounded interval proof terminates in 4,964 boxes. The complete profile classification is inherited; the frozen interval engine is reused rather than independently reimplemented. Intermediate walls, later supports and ideal endpoints are outside this first theorem.

The [next-stage theorem](fresh_h9_h8_disjoint_walk_20261014_next_gate.md) proves the same result on all 71 new germs, with 1,065 finite profiles, fifty continuous dihedral tails and 3,269 boxes. Combined with the first theorem on the 56 old germs that remain open, this covers all 127 current finite-corner directions. No single compatible finite-endpoint extension there cuts the escape. The [affine theorem](fresh_h9_h8_disjoint_walk_20261014_affine_gate.md) separately covers all fourteen original proper-affine profiles for every real parameter and every old-wall label, in fifty boxes. Ideal profiles on the newer supports remain open.

An additional [endpoint scope review](fresh_h9_h8_disjoint_walk_20261014_endpoint_scope_review.md) independently checks the original 62 support matrices, inverses, finite Schur costs, isolated tails and positive affine kernels. It audits the orientation and parameter completeness arguments and the angle relaxation. It explicitly inherits the physical coefficients and interval-tree arithmetic, and does not review the later 71 supports.

The narrower [source-corner theorem](fresh_h9_h8_corner_endpoint_gap_20261014.md) has separate Fraction and SymPy checks for all 41 right-contact corner families, at every real scale. It remains useful as a small independently computed obstruction.

## Changing the three cap heights

The [maximal-height construction](fresh_h9_h8_cap_height_20261014.md) uses the form `Q8⊕<(7+√5)/4>`. Its fifteen actual facets have exactly 62 finite vertices and 220 finite edges, leaving 118 finite-corner germs. The [independent review](fresh_h9_h8_cap_height_20261014_review.md) obtains the same complete finite census by testing all 5,005 nine-wall subsets and solving physical nullspaces.

Changing these three heights cannot finish the chamber. Normalize each exceptional projected direction to minimum nonzero source contact one. Its cap norm must satisfy `0<q≤1/2` or `q=1/(2cos²(π/m))` for an integer `m≥3`; hence `q≤2`. After a common coordinate rescaling, each last-coordinate magnitude is at most one. A literal strict escape for the maximal case remains strict when any height decreases. This proves infinite volume for all independent height choices on these three fixed projected directions, without a label cutoff. New projected directions remain open.

The [maximal-height corner theorem](fresh_h9_h8_maxheight_corner_gap_20261014.md) additionally excludes an escape cutter in every right-contact source-corner family, at every scale. Each of its 41 inequalities has separate Fraction and SymPy checks. The maximal-height form has positive-definite nontrivial conjugate and vertical-coefficient norm `11/4`; neither fact alone establishes reflectivity.

There is nevertheless a concrete arithmetic route in this changed form. The [integral-host certificate](fresh_h9_h8_cap_height_20261014_integral_host.md) rescales the bottom normal by two and checks that all 225 Gram entries and 225 reflection coefficients are in `Z[φ]`. Their fifteen generators span a rank-ten integral reflection module. An [explicit free basis](fresh_h9_h8_cap_height_20261014_host_basis.md), with both module inclusions checked, identifies it as `E8 ⊥ [[−2φ,−1], [−1,2]]`. Its determinant is `−φ(3+φ)`, with field norm `−11`.

The [primitive-root norm theorem](fresh_h9_h8_cap_height_20261014_norm_pool.md) reduces its complete arithmetic reflection search to norms 2 and `2(3+φ)`, up to unit rescaling. Norm-two vectors automatically give integral reflections; the second class additionally requires the stated divisibility filter. Evenness, the prime determinant, two primitive Bezout arguments and the explicit unit calculation prove completeness. This prepares a new Vinberg search; it does not enumerate new roots or establish reflectivity. In particular, no theorem for unimodular hosts is applied to this lattice.

The [independent arithmetic review](fresh_h9_h8_cap_height_20261014_norm_pool_review.md) reconstructs both module inclusions, the form, determinant and positivity in standalone Fraction arithmetic. It also simplifies the second norm class's divisibility filter: modulo `β=3+φ`, all eight E8 coordinates vanish and `t≡6x`. This is a concrete sieve for the future enumeration.

## Compact star search

The [dependency extension](fresh_h9_star_single_right_20261014.md) frees the own contacts of both cut leaves 15 and 16, while retaining its explicit constraints on the positive leaves, core and internal labels. The [independent review](fresh_h9_star_single_right_20261014_review.md) checks all 4,956 actual-vertex witnesses and all 23 bound reductions. It also checks a literal reversal of all 378 vertices and all 23 minimal nonfaces. Transferring the certified predicates under that reversal increases their number from 161 to 173.

A further [finite-core theorem](fresh_h9_star_single_right_20261014_finite_core.md) fixes only the six own contacts at core leaves 11, 14 and 17, all at pi/3, and allows every noncore own contact to vary freely. It excludes all sixteen choices of internal labels 4 or 5. The new dependency guards inherit the frozen 21,952 nonzero rank determinants rather than repeat that bank. The full star incidence remains unresolved.

The [independent finite-core review](fresh_h9_star_single_right_20261014_finite_core_review.md) checks the 27 guards, three budgets and full inherited key coverage. The composed binder has 205 predicates: 173 previous predicates, sixteen finite-core cases and sixteen literal reversals. The [next saved batch](fresh_h9_star_single_right_20261014_finite_core_search_batch_009.json) reaches 1,269 cumulative visits with nine pending frames. Its next state has internals `(4,4,5,≥7)`, opposite-single own profiles `(2,4)/(4,2)` at leaves 8/16, and shared-single profiles `(2,3)/(2,3)` at leaves 11/17. Other own profiles remain 3/3. These changed core contacts and the unbounded last internal are outside the latest finite-core theorem; the positive-leaf right contact also lies outside the preceding unrestricted-cut theorem.

## Reproducibility and continuation

Each exact construction or theorem has a source-bound certificate and a saved replay; numerical scouts retain their bounded scope. The [checkpoint audit](disjoint_cap_edge_extensions_and_star_symmetry_2026_10_14_audit.json) verifies recursive source hashes, preservation of the previous frozen archive and local report links. It is an integrity audit, not an additional mathematical proof.

Continue on genuinely new wall directions and supports, or a changed retained chamber. The current escape and the endpoint obstructions prevent treating a larger finite-vertex count as evidence of finite volume. Preserve the distinction between arithmetic host questions, exact geometric constructions and bounded numerical proposals.
