# Right-angled dimension eleven: proof map and reproduction

19 October2026. [Current theorem and verification record](../results/right_angled_mixed_dimension11_20261019.md).
The completed internal mixed ceiling is11. The focused
[11-page PDF](../manuscripts/build/right_angled_dimension11.pdf) contains
the complete geometric, scalar and counting argument with three vector
figures. Two full manuscript reads, the zero-warning bounded build and
all three actual figure-page visual checks pass. The earlier
corner/transport sources remain supporting predecessors in this same
mixed-paper family.

|Step|Complete written source|Boundary|
|---|---|---|
|Actual canonical P3 faces and opposition graph|[Weighted-triangle geometry](../results/right_angled_mixed_weighted_triangles_research_20261018.md), [opposition invariants](../results/right_angled_mixed_opposition_moments_research_20261018.md)|All-wall feasibility, actual two-finite-vertex bipyramids, fixed componentwise h and null vector; written geometry.|
|Strict weighted face count|[Weighted-triangle proof](../results/right_angled_mixed_weighted_triangles_research_20261018.md)|DS and actual links, published82 facets in dimension8, elementary face minima.|
|Strict lower2B/f>184460/1619|[Full polygon/three-face coupling](../results/right_angled_mixed_polygon_three_face_coupling_research_20261018.md)|Q8 attachment and exceptional-incidence geometry reviewed independently. Lean checks its final incidence/E/R algebra.|
|All220 formal-triple moments and ratio domain|[New full opposition proof,Sections1–2](../results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md)|Pointwise at every finite vertex; exact monomial expansion is independently replayed.|
|Full-domain quartic positive-part majorant|[New proof,Section3](../results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md)|Rational Bernstein and convexity argument; all finite coefficients replayed. No sampled inequality.|
|Deletion, reciprocal pairing, ordinary upper2B/f<=1420561/12500|[New proof,Section4](../results/right_angled_mixed_positive_part_opposition_cap_research_20261019.md)|Actual subsets and one null vector per component; written real-valued argument. This is an explicit analytic premise in Lean.|
|Contradiction and dimension11|[Current theorem](../results/right_angled_mixed_dimension11_20261019.md)|Conditional count kernel passes; published Dufour12 supplies the higher tail.|

Root and a separate reviewer read the complete new proof and inherited
lower-count geometry. Their reviews are linked in the theorem record.
The earlier candidate source is pinned and preserved; its pending-status
header records the original handoff, superseded by the current theorem.

## Exact certificate

Source checker:
[check_right_angled_mixed_positive_part_opposition_cap_20261019.py](../scripts/check_right_angled_mixed_positive_part_opposition_cap_20261019.py),
SHA256 `50627366713d0a8f2417fc0ed7251a62a8572a702b79a5389d3a2130eb663892`.
The3333-byte saved certificate binds the new source and four earlier
reports. It independently expands the primitive polynomial and220
subset sums, rather than trusting a copied coefficient table.
Generation and separate replay have passed. Reproduce with:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_right_angled_mixed_positive_part_opposition_cap_20261019.py --verify
```

No optimizer, matrix census, enumeration of polytopes or scientific
third-party package is needed. The complete proof that the recorded
coefficients imply full-interval positivity is written in the source.

## Conditional Lean

Source: [RightAngledMixed11.lean](../formal/CoxeterProofs/RightAngledMixed11.lean).
Compile and audit only this module:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledMixed11 --audit-axioms
```

The checked theorems are polygon_coupling,canonical_sum_lower,
strict_opposition_lower,exact_comparison,no_data12 and
dimension_at_most11. The compilation and audit ledgers are under
`formal/.proof-batches/CoxeterProofs_RightAngledMixed11*.json`.
Only propext and Quot.sound are printed; exact_comparison is axiom-free.
All project dependency provenance is capped, with no legacy-cache gap.

`Counts12` exposes the actual Euler/coface identities, strict weighted
row, nonnegative excess and both exceptional-incidence comparisons.
`Data12.positivePartTransport` explicitly assumes the cleared aggregate
real-analytic comparison obtained by the reviewed majorant/moment/deletion
proof. Lean does not declare actual ratios integer-valued, construct
hyperbolic objects or prove the full-interval polynomial inequality.
It checks the final implications after these named interfaces are supplied.

Both commands use the shared serial lock and the40GiB disk/6GiB RAM
launch gate. Research jobs retain59s/1GiB caps; Lean uses one module,
one thread and2048MiB. The cleanup trigger remains20GiB. The new
module was checked after complete mathematical review and exact replay.

## Current complete package

The [source-bound execution record](../results/mixed_dimension11_package_verification_20261019.json)
collects three exact-check sources, the new certificate's separate
generation/replay, the conditional kernel/audit and final PDF/visual
records. The final source is
`ccb9615aa7ea82ed503b179a443c1ffd056d6df94657dcef749ff3f61dc3df1b`;
the PDF is
`d00d8d8f4a927281f7a0a187351b440c7de95ad2680f8d9fae984e8338a6c779`.
The final11-page build took1.649s with zero warnings and peak RSS below87MB.
All TeX intermediates and all nine raster previews have been removed.

To reproduce the declared checks and a temporary PDF in one serial plan:

```sh
python3 -B scripts/verify_manuscript.py mixed11 --mode core --pdf
```

The driver plan has been inspected and includes exactly three bounded
script replays, one conditional Lean module with axiom audit and one
temporary PDF build. Each component already has its source-bound PASS
record above; the plan inspection itself was not labelled a new run.
