# Independent review of the consolidated compact RA4 translation

24 October 2026. Reviewer and coordinating acceptance: `/root`.

**Accepted as a complete written mathematical translation at its named
published-input and fixed-model boundaries.** The revised source has not been
compiled. This acceptance supplies neither a new execution nor end-to-end
geometric Lean verification, external refereeing, or a minimum-volume theorem.

## Exact source and independence

The reviewed translation is
`manuscripts/revisions/compact_ra4_volume_constraints_consolidated.tex`,
1,785 lines, SHA256
`f51b8c9f5752a7c218b5191589e032333be31fe5017300a8f24d2872b04495ab`.
Its author was `/root/audit_finite57_algorithms`; root did not supply edits to
this translation before review. Root read the entire source, including all six
TikZ figures and the algorithm/external-input discussion, and independently
reconstructed the mathematical steps below. Figures have only a source review:
there is no rendered-page or visual PDF review for this revision yet.

Root also read both complete fixed-model workers: C26, 303 lines,
SHA256 `47108a5ff3f732bf4b54d7893cc925f8e2a9423106c0ac812f690d355a9976e4`;
and C28, 427 lines,
SHA256 `7356fa22ba8095e480c81cd2800f62832764925ffc6c0706be904ec8c0a6bf4e`.
The full source ledger has SHA256
`b56effe49cdf5dfad45d72f0559093b4f93a82c8cf7c80fc3895eec20e6ae80a`.
The preceding 329-line source audit has SHA256
`29880f5978cdb91b95cb95b8fb7fc58357078f41ebd11557981fac7b168b7902`.
These bytes were freshly matched. The historical 27-member freeze and four
runtime authentications retain the preceding audit's scope; root did not
rerun the workers or claim a fresh transitive authentication in this review.

## Shared geometric and counting arguments

The positive-definite two-wall projection has the correct sign and checks
EVERY original wall. It gives an actual non-right ridge for products in
(-1,0); a product -1 gives a nonzero feasible future null, excluded by
compactness. Thus distinct walls have product zero or below -1. Orthogonal
collections meet in actual faces by the same projection. Incident finite
normals then prove simplicity, and right-angled face inheritance is justified
in this particular class, rather than assumed for arbitrary Coxeter faces.

For an outgoing perpendicular edge from a facet vertex x, the normal
`sinh(L)x-cosh(L)a` is unit, has the correct inward side, and is disjoint
from the original facet. Orthogonality to neighboring normals is equivalent
to incidence at x. Equal such normals force equal x. A proposed second actual
wall on the same perpendicular line would have positive product with the
first. These facts justify all distinctness and exact-incidence claims in the
capacity count.

The polygon proof treats both zero extra faces and a unique extra face that
meets an old neighbor. In each case the complementary dual disc has no
interior vertices and an ear with total degree below five. The equality
continuation rules out chords and proves the full Loebell lattice. This is a
proof of the required equality type, not a reference to an unexplained list.

In dimension four an extra facet's footprint is a clique. Three members would
force the already-counted outgoing opposite, so a genuine extra serves at
most two neighbors. An index needing only one extra consumes a distinct
singleton, by the polygon equality argument. The resulting incidence sum is
at least M extras and yields `f>=1+V+2M=2V+5`. It does not assume the relevant
ridge contributions are disjoint.

Trivalence, simplicity and Euler give `W=v-5f`, `d(F)=(V-20)/2=M-12`, and
`sum d(F)=2W`. Neighboring excesses are summed over distinct facets, with
shared ridges correctly retained in the latter identity. Thus
`W>=ceil(3(V-20)/4)`. The reflection-orbifold stabilizer orders give the stated
Euler characteristic and, using the named Gauss--Bonnet input,
`Vol=pi^2(f+W+16)/12`. The 46/96 local cutoffs, including their strict and
nonstrict endpoints, follow from the displayed integer calculations.

## Recognition and the full excess gap

Zero excess makes every facet a regular dodecahedron at the named RA3
recognition input. Matching whole facets, not only their hyperplanes,
propagates across actual pentagons; the inward side fixes the choice. Closure
under all neighbors in the finite model and both connected facet graphs give
equality of the complete halfspace systems. This proves recognition, without
propagation from a single dodecahedron through arbitrary larger facets.

The graphs G and H have distinct roles: G records nonpentagonal ridges with
their integral excess weights; H includes all positive-facet adjacencies.
The no-induced-square argument applies to H. A degree-two G vertex has equal
weights and the full Loebell type. The lateral metric calculation for L6
identifies its unique longest edge. If another lateral facet also has excess
two, compact simplicity identifies its hexagonal-base parent correctly.
Counting the two disjoint rows and bases then raises the bound to W>=20.

Three consecutive regular-length edges close a polygon because the two
orthogonal closing lines are actual walls; the all-wall intersection lemma
prevents truncation. The eligibility lemma uses this fact, and the four-vertex
graph proof is exhaustive. Root checked both branches of the nine-vertex
complement argument, including its two isolated-vertex case.

The five-neighbor proof at excess three exhausts graphs with a universal
vertex, with a triangle and with neither. Its marked-hexagon restrictions use
the actual boundary orders. The P5/C5 endings require the pentagon metric,
so the argument does not assert a false purely abstract ten-vertex
obstruction. Root also checked the t=2,3,4 cases in the ten-vertex complement
lemma. With `3N<=2W<=32`, these yield N=10 and the contradiction proving
`W=0 or W>=17`. All prerequisites occur earlier in this manuscript.

## Farthest facet and the classification-to-metric branch

Projection of the SAME original interior point to a maximizing wall checks
all other walls and lies in the actual facet. All intrinsic heights are
strictly below one. The Legendre polynomial expands exactly as printed; root
recomputed its coefficients and the trace-free tensor contractions through
degree four. The bound N<=14 and its strict equality contradiction give
N<=13. Polygon capacity then forces twelve pentagonal faces and the named
RA3 recognition supplies the dodecahedron.

The paired icosahedral normals give the maximum-height estimate, not a
minimum-height estimate. Root checked `h^2=phi/2`, `T^2>=phi^3`, and the
equality propagation at the SAME point through all maximizing neighboring
facets. No inscribed 120-cell or lower bound on all wall heights is inferred.

The C26 count-one input is applied only after capacity excludes seven-sided
faces. Its explicit dual sphere and actual finite symmetry group justify the
normal ansatz. Root rederived `U^2=2/3`, `beta^2=1/2`,
`(2+3X)(3X-1)=2`, and both distinct edge witnesses. Their physical signs and
strict comparison with the regular pentagon edge exclude both pentagon
orbits.

At C28 the heptagonal equality case is treated separately as L7. The two
explicit fullerene dual spheres have different hexagon-induced edge counts;
the named published count two proves coverage. In the tetrahedral case the
actual fixed-point representation, positive spatial Gram and stabilizers
justify the coordinates; root rederived `B^2=3/4`, `a^2=7/4`, `h^2=3/4`
and the edge cosh 4/3. In the other type, either outer pentagon orbit forces
all twelve regular and contradicts the hexagon closing lemma.

For the remaining orbit the actual Lorentz basis is nondegenerate. Root
recomputed H1 from its four products, its norm equation and its product with
H3. The quadratic in k^2 has exactly the printed roots; the compact wall gap
rejects the smaller one. The surviving `k=phi^2,e=-phi` gives `C.H1=0`,
contradicting checked nonadjacency. This shorter ending is valid independently
of the longer historical fiber proof. It uses actual wall intersection and
never constructs a polyhedron from a partial Gram matrix.

A dodecahedral facet exists. Unless all facets are dodecahedra, connected
adjacency yields a positive neighbor with an ENTIRE regular pentagonal ridge.
Excess one, L6, C26 and C28 are all excluded. Thus its V>=30 implies f>=65;
the already-proved excess gap gives `Vol>=49*pi^2/6`. The recognized 120-cell
satisfies these inequalities separately. This completes the new revision's
claimed universal constraints. It establishes neither its minimum-volume
target nor the proposed compact second-volume double.

## Algorithms, explicit inputs and publication gate

The workers inspect the complete stated fixed triangle sets, cyclic links,
symmetry actions, edge witnesses and exact quadratic-field identities.
Finite permutation closure and monotone subset enlargement halt. Connected
cyclic links, two triangles per edge and Euler characteristic two establish
the stated spheres. Rational-pair arithmetic and its sign algorithm are
exact. The manuscript states their complete finite domains and what their
outputs mean. These are fixed-model algorithms, not fullerene enumerations.

The external inputs remain acute-wall geometry, foundational Lorentz/convex
and surface facts, orbifold Gauss--Bonnet, twelve-face RA3 recognition,
marked compact nonobtuse RA3 metric uniqueness, and the precise two published
fullerene counts. The preceding audit's primary-source checks retain their
scope; no new literature proof verification is claimed here. General compact
RA4 geometry assumes no arithmeticity.

No additional mathematical repair was identified. A publication source may
retain these exact mathematical bytes while changing only the draft-status
banner/date and the final review-status sentence. Bind that editorial delta
and use a NEW filename directly under manuscripts; do not overwrite the
pinned canonical path. A matching bounded PDF build and rendered review
remain mandatory before replacing the reader-facing edition. Any conditional
Lean addition must separately list its geometric, classification and script
parameters and receive compilation plus an actual axiom audit.

No scientific worker, enumeration, test, Lean or TeX process ran for this
review. Existing successful model runs were not repeated. This is another
accepted revision inside the existing six families, not a seventh paper.
