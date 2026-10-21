# Star contact families and H8 attachment obstructions

**No new finite-volume Coxeter polytope or dimension is established.** This phase advances the unresolved compact-H9 star search, excludes an entire continuous attachment construction, and classifies the rank-four couplings of the fixed mixed H5 source. The [preceding checkpoint](general_endpoint_caps_and_compact_case5_2026_10_14.md) and its frozen artifacts are preserved.

## Compact H9: exact contact-family exclusions

The search targets literal star case 0 in the saved eighteen-facet tree catalogue: 378 actual vertices, eight missing triples and fifteen dotted pairs. It starts with 132 transferred structural gates plus the previously proved star conditions. These placements are checked against the literal incidence; no incidence is silently identified with another catalogue case.

Three successively resumed batches preserve every pending DFS frame and sibling. They reach necessary assignments at cumulative visit counts **21, 40 and 51**. None is a metric realization. Their exact family obstructions are broader than exclusions of the individual tuples:

* [All own contacts π/3](fresh_h9_star_zero_20261014_own_three.md): all four internal labels are unrestricted finite integers. Missing triples and reviewed high-label gates reduce them to sixteen choices in `{4,5}`. The complete Schur census has **21,952 nonzero determinants** and **784 negative anchor determinants**, forcing rank eleven in the ten-dimensional ambient space of H9. A [separate Fraction review](fresh_h9_star_zero_20261014_own_three_review.md) reconstructs every determinant, all thirty vertex guards and all 128 support-budget checks.
* [Opposite single contacts, last internal label at least four](fresh_h9_star_zero_20261014_opposite_single.md): first three blocks retain own contacts `(3,3)`, while the last anchors have opposite single contacts with arbitrary labels at least four. Literal gates reduce internal labels to `{4,5}`. Of **896 continuous polynomial cases**, 832 have all **7,488 Bernstein coefficients** strictly positive; the other 64 fail an actual vertex at the minimum own labels. The [independent interpolation review](fresh_h9_star_zero_20261014_opposite_single_review.md) reproduces all polynomial and Bernstein tables and the vertex failures.
* [The separate last-internal-three branch](fresh_h9_star_zero_20261014_opposite_three.md): missing triples force both opposite single contacts to have labels at least seven. Nineteen actual supports force a residual Lorentz plane and an orthogonal negative third vector, contradicting negative index one. Its [independent literal and symbolic review](fresh_h9_star_zero_20261014_opposite_three_review.md) passes. This case was checked separately before the combined search gate was used.

The [51-visit saved state](fresh_h9_star_zero_20261014_opposite_search_batch_003.json) has internal labels `(4,4,4,4)`, with last own profiles `(2,5)` and `(4,3)`. The first three blocks still have own contacts `(3,3)`. The whole incidence remains open; omitted meeting labels, dotted distances, Lorentzian signature and compactness must all be checked before any survivor is an example.

## H8 attachment: an all-real contact gap and whole-line exclusion

The [uniform source theorem](fresh_h9_h8_facet_cover_20261014.md) concerns the fixed compact Bugaenko H8 chamber with eleven facets and 41 finite vertices. Set `a=(√5−1)/4`. Any future timelike projected wall with contacts in `{0}∪[1,∞)` either has norm magnitude squared at least `4a`, or lies on one of exactly three source-vertex rays. The proof enumerates all **43 independent nine-root subsets**, with 512 binary contact choices per subset, obtaining **3,657 distinct feasible vertices in 1,411 zero patterns**. Convex decomposition and the reverse Lorentz triangle inequality extend the finite calculation to all real projected normals. The three exceptional patterns are one-dimensional rays, with minimal norm magnitude squared a. No finite-label cutoff or numerical sampling is used.

This already excludes the literal cap in two ways: four distinct old centers would occupy only three source vertices, and wall 7 would require a contact strictly between the consecutive permitted values one and √2.

The [exact parameter theorem](fresh_h9_h8_facet_cover_20261014_parameter_gate.md) then restricts any possible finite-volume attachment along the frozen H–E segment to

`0.99674949 < t < 0.996752404`.

It reconstructs all 41 source vertices and encloses 820 distances, using outward dyadic arithmetic. A proved monotonicity formula makes every ball radius and pair of center distances increase along the all-dotted segment. At the lower endpoint the entire ball union is too small to contain the source; at the upper endpoint two mandatory centers are too far apart to fit in it.

The [whole-line theorem](fresh_h9_h8_facet_cover_20261014_line_exclusion.md) splits the remaining interval at `T=0.9967505`. Before T, wall 7 has norm in the forbidden interval `(a,2a)`. After T, four distinct centers all have norm below `4a` and would occupy only three exceptional rays. Both inequalities include T. Consequently **every isometric placement of this source fails for every all-dotted cap on the specified segment**. The [independent scalar and proof review](fresh_h9_h8_facet_cover_20261014_line_review.md) verifies the derivative signs and split margins in separate 192-bit arithmetic; it explicitly depends on the source-gap census and does not independently repeat the full source-distance census.

This stops one concrete route. Other cap directions, other source chambers, additional source facets, nonperpendicular attachments and changed retained walls remain outside the theorem. Preliminary assignment and local-placement diagnostics are archived separately and are not premises of the exclusion.

## Fixed-source couplings and the B4 fold

The [complete coupling census](fresh_h9_alternative_mixed_couplings_20261014.md) examines 4,096 labeled spherical rank-four candidates, reducing 182 positive diagrams to eighteen types. It checks 216 positive-projection profiles and 5,168 pairs at the fixed source's π/5 edge. Uniform inequalities handle both high support labels and high contact labels, so the finite census is not an assumed angle cutoff. The [independent Fraction review](fresh_h9_alternative_mixed_couplings_20261014_review.md) reproduces the support, profile and pair tables.

The π/5 edge forces both source endpoints ordinary, and the π/3 chain propagates this to facets 0–5. Facet 6 can have zero or unit support-weight square. Only D4 and B4 give irreducible nontrivial couplings; every allowed reducible choice leaves an entire support component untouched. This classifies only the fixed seven-facet H5 source with four spherical support facets and one lifted wall per source facet.

The [exact B4 fold](fresh_h9_alternative_mixed_couplings_20261014_b4_fold.md) has an eleven-wall seed with ten finite and twenty spacelike extreme rays. Only one of its twenty polars passes all retained-wall checks. It is exactly the inherited D4 polar up to positive scale. The resulting twelve-facet chamber has twelve finite and 38 spacelike extreme rays; none of its 38 polars passes. A strict exact escape proves infinite volume. Arbitrary caps remain outside this exhausted polar family.

## Verification and continuation

Exact producers and their saved replays pass, with independent reviews as scoped above. The [archive audit](star_contact_families_and_h8_attachment_obstructions_2026_10_14_audit.json) binds the source manifests and local report links and preserves the preceding frozen archive. Archive integrity alone is not a mathematical proof.

Continue the star search from its saved cursor with broader proved last-block contact families. For the constructive route, change the cap direction and use the exceptional source-vertex rays as exact anchors before attempting placements. Repeating arbitrary-height endpoint banks, the fixed H–E attachment line, or the exhausted B4 extreme-ray-polar procedure cannot improve their established outcomes.

```sh
python scripts/audit_star_contact_h8_attachment_checkpoint_20261014.py --verify
```
