# Exact finite-vertex incompatibilities for the positive raw23 prism family

The 57 positive actual raw23 prism classes of the current arithmetic witness give **294 necessary incidence inequalities**. Exactly 120 of them violate the preceding vertex allowances. Keeping every face mass fixed and increasing only those 120 vertex allowances gives the exact repaired ratio **8.804856730635848…**, down from 11.986190572394824… after the [first six-orientation certificate](fresh_dominant_cap_compatibility_20261017.md). This invalidates that fixed witness for the enlarged constraint family. It is neither an optimized feasible rule nor a dimension theorem.

The [source](../scripts/fresh_dominant_cap_compatibility_family_20261017.py) and [saved exact data](fresh_dominant_cap_compatibility_family_20261017.json) cover every finite corner of those 57 full-star classes. No optimizer, historical 94,722-column matrix reconstruction, or 2.7-million joint-family scan is run here. The argument uses actual normalized wall geometry at a common finite vertex, beyond the earlier separate local simplex packing and prism-base capacities.

## Complete selected domain and literal maps

The selected classes are precisely the positive actual-provenance rows whose canonical prism core has raw source 23 in the frozen input. There are 57 such rows. Each is a 5-prism with ten finite vertices and ten quadrilateral subfaces. Every one of its ten distinct two-root deletions is reconstructed, giving 570 actual corner occurrences in 299 rank13 finite-forest types. The old model's injective actual-provenance-to-face map is retained throughout; no merged row is represented by a newly chosen provenance.

Raw23 has seven roots. Its nonright products are the dotted pair 0–2 with squared normalized magnitude

\[
\tau=(7+\sqrt5)/2,
\]

the simple edges 2–3, 3–5, 5–4, 4–6, and the 5-edge 6–1. At a finite corner, the five retained raw roots form a positive finite forest; the other eight incident roots are the saved orthogonal spectator forest. Every weighted connected component of the retained raw forest is mapped onto a whole matching component of the canonical rank13 vertex forest, using every choice among identical components and every internal weighted isomorphism. The remaining spectator permutations change neither external covectors nor occupied quadrilateral port pairs, so they may be omitted. The source checks the full active incidence map and orthogonality to every unused component.

Each literal star is identified by the **unordered pair of complete external contact covectors** in this fixed thirteen-root frame. Identical keys are deduplicated before any exclusion test. The corresponding occupied quadrilateral port pairs must agree. Complete corner-orientation orbits are checked to be identical or disjoint, and no orbit is allowed to combine distinct selected Q classes. Distinct corners of one prism remain distinct global vertex occurrences; they are not counted as different simultaneous presentations of one local star.

The result is 570 corner orbits and 6,010 distinct literal star orientations, counted separately at their respective vertex forests. All map witnesses and actual omitted-root pairs are saved. The package reproduces the preceding inequality at vertex3137 with its original coefficients1 and2.

## The common finite-vertex equation

Let \(v\) be the future unit timelike vector representing a finite vertex and let \(V\) be the positive definite Gram matrix of its thirteen incident norm2 wall normals. For a genuinely external actual wall with full covector \(p\), write

\[
s(p)=\sum_i(V^{-1}p)_i n_i,
\qquad h(p)=p^TV^{-1}p-2.
\]

The supporting halfspace is strict at this vertex, so \(r\cdot v<0\). Consequently \(h(p)>0\), and the outward normal is uniquely

\[
r=s(p)+\sqrt{h(p)}\,v.
\]

For two such walls the actual product is therefore

\[
r\cdot r'=p^TV^{-1}p'-\sqrt{h(p)h(p')}.
\tag{1}
\]

This fixes the relative timelike sign; it is not an independently adjustable parameter for each local prism. Equal complete covectors reconstruct the same outward wall and are explicitly skipped by the distinct-wall incompatibility tests.

All active raw23 components have labels3 and5, so their inverse Gram entries lie in \(K=\mathbb Q(\sqrt5)\). Passive B/F or other spectator components can have normalized products outside K, but every external covector is zero on them. Only the active blocks are inverted, and the support assertion is checked explicitly. Each covector entry lies in \(K+K\sqrt\tau\), while every squared height is checked to lie in K and to be strictly positive. For each literal star, equation(1) is also checked against its prescribed product between the two omitted roots. Thus every exported orientation uses the same future side.

The comparison is exact. If \(E=a+b\sqrt\tau\) and \(H>0\), the sign of \(E-\sqrt H\) is negative when \(E\le0\); when \(E>0\), it is the sign of

\[
a^2+b^2\tau-H+2ab\sqrt\tau.
\]

The same test with a replaced by a+1 determines whether the product is greater than−1. Rational arithmetic and exact signs in K suffice; there is no floating-point angle tolerance.

## Three exact pair exclusions and the clique inequality

Two distinct literal corner stars cannot coexist at the fixed vertex if any of the following holds:

1. They occupy a common actual quadrilateral port pair. A pair of incident edge ports determines one actual two-face. Its full wall star is unique, whereas the two external covector-pair keys are different.
2. Some pair of distinct reconstructed external walls has positive product, violating acuteness.
3. Some such pair has normalized norm2 product strictly between−1 and0. This lies in the universal Coxeter gap between a right angle and a 3-edge.

The source leaves every pair not excluded by these gates compatible. It does not assert that such a pair is realizable. Short-circuiting stops a failed capacity-one test; every pair needed to certify an exported clique nevertheless has a saved exclusion witness.

If all distinct literal orientations in one complete corner orbit are mutually excluded, at most one actual occurrence of that orbit can be present at the vertex. Identical orientation data are one physical occurrence. Two such orbits are joined only if **every** cross-orientation pair is excluded. Any clique of this orbit graph therefore has local capacity1. Summing over all vertices of forest type Y gives

\[
\sum_Q c_{Q,Y}\,N_Q\le N_Y,
\tag{2}
\]

where \(c_{Q,Y}\) counts actual corner occurrences of Q in the chosen orbits. It is not an orbit size or the number of embeddings. All maximal cliques of the retained orbit graph are enumerated, and identical coefficient columns are deduplicated.

The scan proves capacity1 for560 of570 corner orbits, with378 universal cross-orbit edges. The ten omitted orbits occur in five forests with two H4 components; their unresolved pair tests are not converted into inequalities. In total, 134,445 literal pairs were tested. Of these, 134,435 have exact exclusions: 22,397 shared-quad witnesses, 94,464 open-angle-gap witnesses, and17,574 positive-product witnesses. The ten remaining tested pairs remain allowed.

## Exact repaired witness

The exported 294 columns have `vertex_cost=1` and `positive_face_coefficients=[[face,multiplicity],...]`. Their loads are computed against the unchanged 2,524 positive face masses of the [scaled-core edge repair](fresh_sqrt5_scaled_core_cap_edge_repair_20261017.md). The baseline vertex masses include the preceding six-orientation update at3137. For each forest, the new allowance is the maximum of its baseline allowance and all new column loads. Exactly120 allowances increase.

For example:

| Vertex type | Necessary incidence inequality |
|---|---|
| 2A1+A4+H3+H4 (1844) | \(2N_{26410}+2N_{26448}\le N_{1844}\) |
| 3A1+A3+H3+H4 (982) | \(2N_{26328}+N_{26402}+N_{26448}\le N_{982}\) |
| 3A1+A3+A4+H3 (944) | \(2N_{26328}+N_{26394}+N_{26410}\le N_{944}\) |
| 3A1+2H3+H4 (1325) | \(2N_{26362}+N_{26448}\le N_{1325}\) |

The largest additional allowance increase is approximately0.0236855810 at1844. Their total increase is approximately0.3743885423. The exact objective remains

\[
\frac{3842394855727}{309375000000},
\]

and the new normalizer is approximately1.4105694777192301. Their exact quotient, stored as a rational in the JSON, is8.804856730635848… . Every new column is checked against the repaired allowances. Since all face masses are fixed and vertex allowances only increase, every inherited signed transfer, complete local-packing cover, joint-port inequality, cap constraint, and forced-face incidence remains valid by monotonicity. No changed-alpha cover is being reused.

This is a **fixed-witness repair**, not the optimum of an enlarged LP. The selected positive raw23 cohort is complete for that witness; it is not a complete catalogue of all possible geometric face classes. In particular, the calculation supplies neither a target12 primal certificate nor a compact arithmetic dimension ceiling.

## Validation and reproducibility

The final source and data are frozen at:

- PY SHA256 `330809a628b4fa85f6d4d81f14d652a721874000bbfbd11c5e7267e6d3f1c37b`.
- JSON SHA256 `8deaef7fa9ff7138c317f61509d953f3cf5a3c25a2229fc37e1e880aec23b839`.

Final generation64150 completed PASS in17.37 seconds. Fresh saved replay91910 completed PASS in14.86 seconds: `--verify` reran the complete selected embedding/projection/pair/clique calculation and compared every saved field. All18 direct source hashes, four local report links, and source syntax checks passed. The compact JSON retains all relevant literal maps and pair certificates without a duplicate full model export.

Construction-audit's independent source/proof read found no gap in active component completeness, passive spectator elimination, the radical-sign comparison, future orientation, identical-wall skips, orbit deduplication, or corner multiplicity accounting. Its separately written checker generation passed all570 actual occurrence records, 6,010 distinct keys and orbit partition, 560 internal capacities, 378 universal orbit edges, 294 clique coefficients and maximality, and120 allowance raises. It independently re-evaluates all82 distinct radical expressions used by the saved geometric sign certificates, without producer imports. This is an independent saved-certificate check, not a second raw projection or embedding census.

Root's independent Fraction-only check reproduces all294 loads, all120 maximum allowance raises, and the exact normalizer/ratio without importing producer geometry. Root also read the full producer source and geometric argument and found no gap; that read does not duplicate the pair arithmetic. These are internal checks, not external refereeing or formal proof-assistant verification.
