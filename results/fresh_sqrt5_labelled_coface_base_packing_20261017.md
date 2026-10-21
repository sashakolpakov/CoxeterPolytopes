# Simultaneous coface/subset prism-base packing

All **164,219** exact labelled subset-packing columns pass on the unchanged restricted rank-13 witness. There are 12,079 positive loads and **zero violations**; maximum capacity utilization is exactly one. These constraints therefore leave the saved obstruction

\[
\frac{66908102699388132000000000000000}
{5419878674792788039126804769429}
\approx12.3449447329>12
\]

unchanged. This is a restricted-model diagnostic, not a dimension bound or a globally realizable incidence construction. No optimizer, simplex-mass repair or vertex-mass repair was used.

The [standalone source](../scripts/fresh_sqrt5_labelled_coface_base_packing_20261017.py) uses standard-library exact integer/Fraction arithmetic. The [saved data](fresh_sqrt5_labelled_coface_base_packing_20261017.json) include every orbit and check, the actual base mappings, and complete coface slice demands for a subsequent full packing-cone test. Generation and a fresh `--verify` reconstruction pass.

## Geometric inequality

The [shared-base theorem](fresh_prism_simplex_shared_capacity_20261017.md) states that, for a fixed maximal simplex \(S\) and immediate coface \(U\), the vertex subsets \(B\subseteq V(S)\) that are bases of distinct counted maximal dotted prisms are pairwise disjoint. The [weighted version](fresh_vertex_weighted_prism_base_capacity_20261017.md) retains pointwise vertex labels. Here we keep a coface together with an entire vertex subset.

Fix \(T\subseteq V(S)\), \(|T|\ge2\), and an integer \(r\in\{2,\ldots,|T|\}\). Disjointness gives

\[
\sum_B\left\lfloor\frac{|B\cap T|}{r}\right\rfloor
\le\left\lfloor\frac{|T|}{r}\right\rfloor.
\]

This floor weight dominates the indicator \(1[|B\cap T|\ge r]\) with the same capacity. It counts each maximal-prism base once; it does not sum separately over the vertices of that base. Squares use the invariant convention of four edge bases, and higher-dimensional prisms have two simplex bases.

An immediate coface is specified by removing one supporting white root \(q\) of the complete \(S\)-star. Let \(O\) be an automorphism orbit of simultaneous pairs \((q,T)\). Summing the displayed inequality over that orbit gives capacity

\[
|O|\left\lfloor\frac{|T|}{r}\right\rfloor\alpha_S.
\]

A literal base with coface \(q_B\) contributes

\[
\sum_{(q,T)\in O\, :\, q=q_B}
\left\lfloor\frac{|B\cap T|}{r}\right\rfloor\alpha_Q.
\]

The same actual provenance mass \(\alpha_Q\) is used in every column. There is no coordinatewise minimum over unrelated representatives and no independent marginal symmetrization of core \(q\) and \(T\).

## Complete symmetry and incidence reconstruction

The input [flag maps](fresh_sqrt5_flag_refined_base_capacity_20261017.json) retain all 6,016 bases from 1,858 actual quad provenances. Exactly 355 of these provenances have positive mass in the [current five-port witness](fresh_sqrt5_joint_five_port_diagnostic_20261017.md), whose 2,219 positive face masses are unchanged. Each base is grouped by `(provenance, deleted base normal)`; its full distinct canonical simplex vertices are checked against the frozen [unlabelled scout](fresh_sqrt5_coface_integer_base_scout_20261017.md). Its coface root is the saved opposite base facet.

Every relevant complete star has one connected Lorentz component and orthogonal finite components. All 287 Lorentz cores are treated. The source explicitly distinguishes black spherical-deletion roots, which encode simplex vertices, from supporting white roots. It refines these colors using the complete intrinsic weighted graph, then enumerates every permutation within the final color classes. A permutation is retained only if it preserves every exact edge label, including dotted parameters. This is exhaustive because graph automorphisms preserve each refinement color. The retained set is additionally checked for identity and closure.

There are 603 automorphisms in total, from 1,057 candidate permutations. The stored action generates the complete simultaneous `(q,T)` orbits for core coface roots. For cofaces in a finite spectator component, the core subset action is independent of the rooted finite-component action. Capacity includes every root in its rooted orbit and every identical component copy. The positive coefficient still counts the one literal coface of each actual base. This product rule is used only between orthogonal components; core cofaces always retain the simultaneous action.

For each simplex type, the number of core whites plus all spectator roots is checked to equal \(13-d\). All 4,969 represented simplex types and all subset sizes \(2\le|T|\le d+1\) are included.

| Simplex dimension | Columns | Positive-mass coface slices |
|---:|---:|---:|
| 1 | 25,457 | 382 |
| 2 | 23,528 | 101 |
| 3 | 40,757 | 172 |
| 4 | 74,477 | 215 |

## Sum identities and unchanged witness

For every simplex type, subset size \(u\), and divisor \(r\), the capacities sum exactly to

\[
(13-d)\binom{d+1}{u}\left\lfloor\frac{u}{r}\right\rfloor.
\]

For each actual base of size \(k\), the summed positive coefficient is independently checked against

\[
\sum_h\left\lfloor\frac{h}{r}\right\rfloor
\binom{k}{h}\binom{d+1-k}{u-h}.
\]

These identities check simultaneous orbit multiplicities, all spectator copies and all source provenances. Every saved load is computed with exact rational masses. No negative slack occurs. Neither \(\alpha\) nor \(\mu\) changes, so all inherited signed constraints, full local-simplex packing covers and prior port inequalities remain valid without combining different normalizers or reusing covers at changed masses.

## Interface for the full coface packing cone

The output `coface_slice_demands` contains 34,520 slices, of which 870 have positive block mass. A slice records its simplex face, dimension, complete vertex-root order, coface class and multiplicity, exact capacity, source base indices, provenance coefficients for every subset mask, and current rational `block_marginals`.

For a core coface class, choose canonical \(q_0\). A source base at \(q_B\) is transported by **all** core automorphisms sending \(q_B\) to \(q_0\); its subset is averaged over these maps. The resulting measure is invariant under the stabilizer of \(q_0\). It retains the total incidence mass across the coface orbit, so its packing capacity is

\[
|\operatorname{Orb}(q_0)|\alpha_S.
\]

There is no extra division by this coface multiplicity. For a spectator coface, average the base subset under the whole core group; capacity is its rooted-orbit size times its identical-component multiplicity times \(\alpha_S\). These normalizations make both sides transform together. All source bases occur in exactly one slice, and the sum of slice block masses equals the total actual base mass.

The `vertex_roots` list specifies mask-bit order. `block_coefficients` stores `[[mask, [[provenance, rational_coefficient], ...]], ...]`; `block_marginals` stores positive `[[mask, rational_mass], ...]`, with omitted masks having zero mass. A subsequent exact mixture of disjoint nontrivial subsets can be compared directly with the saved `capacity`. At most five vertices occur, so this stronger cone has a small finite packing domain.

Passing the present subset-floor inequalities does not prove membership in that full cone. Intersecting subset families and other packing inequalities remain possible. No such stronger test or new global optimization is claimed in this package.

**FINAL/FROZEN.** Generation and a fresh saved replay pass. Syntax, all 15 direct input hashes and report links pass. A separate internal conceptual/source read found no gap in the simultaneous symmetry, floor inequality or slice normalization; it did not duplicate the arithmetic scan. In particular, a common random automorphism transports every base of a packing together, preserving disjointness, while its linear marginals are the per-base averages exported here.
