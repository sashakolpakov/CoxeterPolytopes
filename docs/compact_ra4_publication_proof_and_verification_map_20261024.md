# Compact RA4 publication revision: proofs and verification

24 October 2026. This is the same compact right-angled four-volume paper
family. Its new [publication source](../manuscripts/compact_ra4_volume_constraints_publication.tex)
has a complete [independent mathematical review and coordinating acceptance](../results/review_compact_ra4_consolidated_translation_root_20261024.md).
The [editorial-delta record](../results/compact_ra4_publication_source_editorial_delta_20261024.json)
binds it to the reviewed draft. **The replacement PDF is not built.** The
existing canonical TeX/PDF and the fixed workers' original source pins remain
unchanged; the reader-facing edition has not been replaced.

The paper proves facet capacity, neighboring excess, the zero-excess
recognition, the full gap W=0 or W>=17, existence of a regular dodecahedral
facet, the small regular-pentagon exclusions, and the universal constraints
f>=65 and Vol>=49*pi^2/6. It does not prove the 120-cell minimum or a compact
second-volume theorem. Exploratory center-cut, inball and homology arguments
remain outside the paper.

## Common method and proof correspondence

| Statement | TeX proof home | Inputs and reason it applies |
|---|---|---|
| Actual compact wall gap and orthogonal intersections | `lem:gap` | Acute-wall property; positive-principal projection checks every original wall; compactness excludes a feasible null |
| Exact outgoing opposites and clique footprints | `lem:opposite`, `lem:footprint` | Actual finite edges, Lorentz signatures and the proved wall gap |
| Polygon capacity and full Loebell equality type | `lem:polygon` and following paragraph | Simple dual sphere, exact opposite incidences and a complementary-disc ear argument |
| f>=2V+5 for every facet | `thm:capacity` | Singleton extras cannot serve a second neighbor; every remaining extra serves at most two |
| W=v-5f, d(F)=(V-20)/2, sum d=2W | `eq:volume`, `eq:facet-excess` | Simplicity, trivalence and Euler; exact face incidences |
| Volume and neighboring excess | `eq:volume`, `thm:neighbor-excess` | Reflection stabilizer orders and Gauss--Bonnet; distinct selected facets, not disjoint ridge contributions |
| W=0 recognizes the regular 120-cell | `prop:zero` | RA3 twelve-face recognition; matching whole facets and connected propagation |
| W=0 or W>=17 | `thm:gap` | Full G/H graph and marked-polygon proofs, including the degree-two and ten-vertex endpoints |
| Every farthest facet is dodecahedral | `thm:farthest-dodecahedron` | One original interior point, all-wall projection, an explicit quartic code bound and polygon capacity |
| Sharp maximum-wall-height inequality and equality | Same theorem | Antipodal icosahedral normals and propagation at the same point |
| No regular pentagon at 26 vertices | `prop:c26` | Count-one fullerene input, explicit sphere, actual symmetry realization and two exact edge witnesses |
| No regular pentagon at 28 vertices | `prop:c28` | Loebell L7 plus count-two fullerene input; tetrahedral metric and full second-type orbit analysis |
| f>=65; Vol>=49*pi^2/6 | `thm:global-lower` | Dodecahedron-to-positive-facet adjacency, every small positive neighbor excluded, capacity and the earlier W gap |

The shared argument branches only after actual geometry supplies a regular
dodecahedral facet. The graph classification is confined to its small positive
neighbors. The G/H excess argument is proved once, before those exclusions.
No arithmeticity or four-dimensional polytope census is assumed.

## External mathematical premises

1. The acute-chamber wall property: distinct outward unit normals have
   nonpositive products. The stronger compact gap and actual intersections
   are proved in the paper.
2. Foundational Lorentz geometry, convex boundary spheres, surface/Euler
   facts, connected facet adjacency and reflection-orbifold Gauss--Bonnet.
3. A compact right-angled RA3 polyhedron with twelve faces is the regular
   right-angled dodecahedron.
4. Marked compact nonobtuse RA3 metric uniqueness, at the scope of
   [Roeder--Hubbard--Dunbar, Theorem 1.4](https://doi.org/10.5802/aif.2279).
   It realizes symmetries of an ASSUMED ACTUAL polyhedron. It does not
   construct a realization from a formal graph or partial Gram matrix.
5. The exact [Brinkmann--Goedgebeur--McKay Table 2](https://arxiv.org/pdf/1207.7010v2)
   counts: one fullerene at 26 vertices and two at 28 vertices, with the
   stated cubic plane pentagon/hexagon definition. The repository does not
   re-enumerate this classification.

Their primary-source retrieval and scope checks are recorded in the
[preceding source audit](../results/audit_compact_ra4_primary_and_c28_fresh_20261024.md).
These mathematical inputs remain assumptions of conditional formal kernels.

## The finite algorithms as mathematical statements

**Fixed graph verification.** Input is one displayed finite triangle set and
specified label permutations. Construct every edge incidence and vertex link.
Check two triangles per edge, connected cyclic links, overall connectedness
and Euler characteristic two. This proves a triangulated sphere. Check the
entire triangle image for each claimed group element. The C26 generator
closure is finite; the C28 actions have 24 and 4 explicit elements. Degree
counts identify the fullerene duals; distinct hexagon-induced graphs prove
the two C28 types nonisomorphic. Exhaustiveness uses the external counts.

**Regular-pentagon propagation.** Input is a subset of the twelve pentagon
labels and the explicit symmetry action. Add every pentagon with two regular
neighbors and close under that action. Every strict step enlarges a subset of
a twelve-element set, so there are at most twelve. Every inference follows
from the proved two-edge pentagon lemma and actual metric symmetry. The two
outer second-C28 orbits close to all twelve labels; the middle orbit remains
DEJK and is excluded by the written normal-coordinate argument.

**Exact algebra.** Input is the displayed finite set of rational-pair
expressions a+b*sqrt(d), d=17 or 5, in its positive real embedding. Rational
field operations and norm inversion are exact. Opposite-sign comparisons
reduce to a rational square comparison. Exponent halving and determinant
order reduction prove termination. These checks verify the stated equations,
factorizations and edge values; universal geometry remains written.

The retained outputs are C26 dual counts (15,39,26), action order 12 and the
two edge values (7+sqrt(17))/4 for squared cosh and (1+sqrt(17))/4 for cosh;
and C28 dual counts (16,42,28), actions 24/4, the stated propagation and the
two roots in k^2. The older C28 certificate retains its longer normal-fiber
ending. The publication proof uses the shorter actual-nonadjacency ending.

The fixed checkers and saved full certificates are bound in the
[draft source ledger](../results/compact_ra4_consolidation_source_ledger_20261024.md).
Each has an actual historical generation and a separate full fresh-process
replay. No new execution occurred to prepare or review this edition. A hash
comparison or `--verify-code` does not substitute for those mathematical
algorithm statements or a full replay.

## Conditional Lean scope and necessary execution

This new edition currently has **no matching compact-RA4-specific Lean
kernel or actual axiom audit**. Existing right-angled and ideal-RA4 kernels
retain their own scopes. Useful selected logic to formalize next is the
capacity/excess arithmetic, graph-count implications, and final positive
neighbor/volume composition. Geometry, classification, exact census outputs
and any unformalized graph lemmas must be explicit theorem parameters. A
printed foundational axiom list does not enumerate those parameters.

After resource admission, the coordinating agent can build this new source
through the existing entry:

```sh
python3 -B scripts/run_proof_batch_20261024.py tex \
  manuscripts/compact_ra4_volume_constraints_publication.tex
```

This is a pending command, not a build receipt. Run one scientific/Lean/TeX
job under the shared lock; require 40 GiB disk and 6 GiB available RAM before
admission, retain the original abort gates and process-tree/RSS monitoring,
and remove only owned unused intermediates after the batch. Review the
rendered PDF and bind the final source/build hashes before changing the
reader-facing edition. Compile any necessary new Lean module separately
from its actual axiom audit. Existing successful fixed checkers need no
rerun merely because a TeX presentation changed.
