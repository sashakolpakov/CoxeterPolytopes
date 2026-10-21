# Compact facet capacity: proof and reproduction map

The [seven-page manuscript](../manuscripts/compact_ra4_facet_capacity.tex)
and [PDF](../manuscripts/build/compact_ra4_facet_capacity.pdf) give a direct
geometric proof of f>=2V+5, where f is the number of facets of any compact
right-angled hyperbolic four-polytope and V is the number of vertices of
any chosen facet. There are two explanatory combinatorial figures.

The proof uses no polytope census or small-facet classification. Its
published background consists of acute-chamber geometry, elementary
hyperbolic polygons and Gauss--Bonnet; it proves the actual-incidence and
planar closing statements internally.

|Internal step|What it establishes|
|---|---|
|Compact normal gap|Distinct walls have product0 or below-1. Orthogonal actual walls have an actual common face by projection from a strict interior point.|
|Perpendicular-edge opposites|Each finite vertex x of a facet has a distinct opposite facet R_x. Its incidences with the facet's neighbors are exactly those through x.|
|Clique footprints|A disjoint facet can meet only a clique of the chosen facet's neighbors, by positive definiteness of a Lorentz-plane complement.|
|Three-dimensional closing|A k-gonal face needs at least2k+2 faces in its compact right-angled three-polyhedron. If there is only one extra face, it meets none of that polygon's neighbors. A dual-disc ear proves both assertions.|
|Four-dimensional capacity|There are1+M+V known facets, where M=V/2+2. Remaining facets meet at most two neighbors. Singleton demands consume distinct singleton extras, so at least M extra facets are needed.|
|Volume consequences|The ridge excess W=v-5f satisfies W>=(V-20)/2, hence f+W>=(5/2)V-5. The orbifold Euler/Gauss--Bonnet formula gives Vol=pi^2(f+W+16)/12.|

A hypothetical polytope of volume strictly below34pi^2/3 has at most48
vertices in each facet; at or below that value the limit is50. Strictly
below68pi^2/3, the limit is104. These are necessary local-size reductions.
They do not prove that the regular120-cell minimizes volume or that its
double has the second volume, and do not imply a120-cell decomposition.

The [direct research proof](../results/compact_ra4_perpendicular_facet_capacity_research_20261018.md)
and [root mathematical/source review](../results/review_compact_ra4_facet_capacity_root_20261018.md)
have separate scopes from PDF compilation. No numerical theorem checker,
complete Lean proof or external human review is claimed.

To reproduce typesetting from the repository root:

```sh
python3 -B scripts/build_illustrated_ra_manuscript_20261018.py manuscripts/compact_ra4_facet_capacity.tex --check-only
```

The shared lock,40GiB/6GiB start gates,59-second/1GiB process-tree cap and
temporary-directory cleanup are built into this serial builder. Check-only
mode does not overwrite preserved artifacts. The
[build ledger](../results/compact_ra4_facet_capacity_build_20261018.json)
records seven pages, two passes, no warnings,84.3MB peak process-tree RSS
and removal of temporary files. The first build encountered a TikZ
foreach parser error; the unfrozen source was corrected before this final
successful build. A PDF build certifies typesetting, not mathematics.

The [freeze manifest](../results/compact_ra4_facet_capacity_20261018_freeze.json)
binds the direct proof, manuscript, PDF, root review, reproducibility map
and build record. The normal, edge and topological interfaces remain
internal proof obligations for future Lean formalization.
