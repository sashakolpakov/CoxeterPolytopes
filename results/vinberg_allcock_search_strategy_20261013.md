# Vinberg, Allcock, and the search strategy

Checked against the primary literature on 2026-10-13. The useful distinction is **distance-ordered wall enumeration versus traversal of chamber edges**. It is not an arithmetic/nonarithmetic division: Allcock's published edgewalking algorithm is itself a lattice algorithm.

The [author-hosted paper](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
is also available when the arXiv HTML endpoint fails. Its introduction explicitly
restricts the method to lattices over Z and describes the remaining work along
each edge as a two-dimensional Lorentzian lattice calculation.

The practical consequence for this repository is to use **both algorithms on
the arithmetic route**, choosing distance batches or complete first-edge
calculations according to the host. The broader construction route is an
**Allcock-inspired diagram-and-metric search**. Its diagram branching is useful
without a chosen lattice, but the arithmetic theorem does not supply its
candidate enumeration or completeness.

Dimension restrictions apply before selecting either algorithm. Esselmann
rules out noncompact arithmetic examples in H20 and every Hn with n>21.
Vinberg rules out all arithmetic examples and all compact Coxeter polytopes
for n>=30. Thus an integral-lattice edge walk cannot discover a finite-volume
H22 example. The active H22 route instead proposes nonarithmetic normals
and verifies their geometry; compact arithmetic H22–H29 hosts over
nonrational fields are not excluded by the dimension-21 theorem. See the
[dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md).

## What the published algorithms establish

Vinberg orders mirrors by distance from a control point and accepts outward roots with nonpositive inner products against the previously accepted roots. Its geometric formulation applies to a specified discrete reflection group; an effective enumeration of that group's mirrors is required. Arithmetic lattices make that enumeration concrete. See [Allcock, Section 3](https://arxiv.org/html/2110.03784v1#S3).

Allcock instead starts at a chamber corner and finds the next corner along an incident edge. It combines spherical diagram extensions with arithmetic in the two-dimensional Lorentzian orthogonal complement. The algorithm must find the first admissible root residue in each relevant lattice coset, compare the endpoint positions, and recover the complete corner stabilizer. The paper explicitly treats lattices over the integers, with rational bilinear form; it does not establish the same algorithm over arbitrary number fields. Its Theorem 9.5 identifies finite volume with exhaustion of the complete unexplored-edge queue for the specified lattice chamber. See [the introduction and Section 9](https://arxiv.org/html/2110.03784v1#S9).

These are chamber algorithms for given input data, not searches through every possible Coxeter polytope or quadratic form.

## What transfers to a combinatorial construction

The following is our proposed adaptation, not a completeness claim from Allcock's paper. In a compact H9 search, begin with an actual finite vertex, delete one of its nine incident walls, and retain the resulting positive definite eight-wall edge diagram. Branch over spherical nine-wall extensions that could describe the other endpoint. Each successful endpoint creates further edges to resolve.

There is a metric variable left after this discrete choice. To see it directly, normalize the edge normals and proposed new normal to square 2. Let B be the eight-wall Gram matrix and p the prescribed eight inner products. Orthogonal decomposition gives

\[
r=r_{\parallel}+r_{\perp},\qquad
r_{\parallel}=\sum_i(B^{-1}p)_i u_i,\qquad
r_{\perp}^2=2-p^{\mathsf T}B^{-1}p.
\]

The complement of the edge span has signature (1,1). A positive residual norm defines a hyperbola, leaving a real parameter for the endpoint position. An additional angle condition with another wall generally reduces this to a quadratic equation. That is the mechanism in the current split-seed scout. Requiring such an extra finite angle is a search restriction: it can miss solutions whose other wall pairs are all disjoint.

Consequently, this route must combine diagram branching with exact metric constraints. It must check every retained wall pair, a strict interior point, actual vertices and facets, and eventual closure of the whole chamber. For finite-volume searches allowing cusps, a complete Euclidean cusp link is needed; a single parabolic extension is not enough. Without a field or other proved bound, a finite list of labels is also a search restriction, since dihedral spherical diagrams I2(m) allow unbounded m.

In particular, unrestricted spherical branching must allow the
noncrystallographic types H3, H4 and I2(m), rather than copying only the
crystallographic normed Dynkin extensions of Allcock's integer-lattice
algorithm. Nor does a chosen extension specify the dotted entries against
remote walls. A construction search must retain those metric freedoms,
propagate compatibility with the whole retained chamber, and backtrack when
the choices obstruct closure. Closing the complete edge queue of a *fixed,
verified chamber* is a different claim from exhausting a chosen list of
extension profiles or closing a few selected edges.

Combinatorial enumeration combined with Gram equations is an established separate approach; [Burcroff's classification](https://arxiv.org/abs/2201.03437) provides a concrete example. An unrestricted version is not supplied by Allcock's lattice theorem.

## What the repository currently implements

| Lane | Implemented | Still needed for a complete edge walk |
|---|---|---|
| [Exact Bugaenko calibration](../src/coxeter_polytopes/search/bugaenko_vinberg.py) | Distance-ordered arithmetic root enumeration with exact pruning at both real embeddings | This is the existing Vinberg calibration, not the new combinatorial search |
| [Quartic arithmetic H9](fresh_h9_quartic_next_frontier_20261013.md) | Five certified first-mirror steps; the first two also connect corners whose integral symmetries prove full-host nonreflectivity | This exact lattice is closed negatively; reuse the interval method on a different host |
| [Changed quartic H9 host](fresh_h9_new_host_next_20261013.md) | Complete first-edge censuses connect integral corner symmetries whose infinite-order product proves full-host nonreflectivity | This second F4 host is also closed negatively |
| [Quartic A9 host](fresh_h9_a9_outer_20261013.md) | Complete first-mirror paths connect two full corner sectors; an integral sector map with trace greater than ten preserves the full chamber and has infinite order | This third quartic host is closed negatively; the norm-orbit interface explicitly depends on certified PARI |
| [Quartic D9 host](fresh_h9_changed_spatial_next_20261013.md) | Eleven complete first-edge steps connect full B3+B6 corners with a D9-preserving transport of trace greater than ten | This fourth quartic host is closed negatively; an infinite-order full-chamber symmetry makes queue exhaustion unnecessary |
| [Tripled-coefficient D9 host](fresh_h9_tripled_d9_20261013.md) | Four new complete first-edge censuses, with all 26 norm classes dividing 24, connect full corners admitting an integral infinite-order transport | This different determinant squareclass is also closed negatively; the proof is recomputed for its actual lattice |
| [Quartic E8+A1 host](fresh_h9_e8_next_20261013.md) | Five complete first-edge steps reach six certified corners; actual lattice involutions at connected full sectors have an infinite-order product | This full host is nonreflective; separate Fraction replays and source review pass, with the certified PARI norm/unit dependency explicit |
| [Quartic A6+D3 host](fresh_h9_a6_d3_20261013.md) | Four complete first-edge calculations in a different determinant squareclass; an integral transport joins full A6+B3 sectors and fixes D3 | Full host nonreflective; separate Fraction edge and transport replays pass |
| [H22 π/7 incidence branch](fresh_h22_incidence_prism_20261013.md) | 835 exact facets and 23 certified finite vertices; conditional obstruction to quasi-arithmetic finite-volume completion | Actual finite vertices supply an edge frontier, but the current chamber has an exact strict escape; its vertex census is not exhaustive |
| [H22 exact frontier and coupled endpoints](fresh_h22_pi7_coupled_caps_20261013.md) | 506 directed edge classifications; two repaired ideal links; two coupled new endpoints give an actual compact quadrilateral in the 839-facet chamber | Exact global escapes and an untouched inherited cusp remain; the covering criterion requires global control of the old chamber |
| [Complete known H21 source walk](fresh_h21_source_corner_walk_20261013.md) | Exact traversal closes all 395 corner orbits and 7,528 directed edge orbits; separate rational replay verifies every endpoint and an empty queue | Complete geometry of a fixed known source, not a search for new mirrors or an H22 completion |
| [Multiquadratic π/6 H22 replacement](fresh_h22_multiquadratic_20261013.md) | 860 actual facets and a compact quadrilateral, independently verified by exact radical arithmetic | Strict escape persists; changing the angle field removes a conditional obstruction without proving finite volume |
| [H22 continuous bridge obstruction](fresh_h22_contact_bridge_20261013.md) | Exact all-real, all-label exclusion for positive-height caps perpendicular to 88 specified edge supports and directly contacting the 305-cap component | Scoped endpoint-family exclusion; nonright profiles, other vertices, indirect chains of new caps, and replacements remain outside it |
| [Coherent H5 edge closure inside H9](fresh_h9_h5_closure_20261013.md) | 123 actual facets, 26 finite H5 vertices and one newly closed edge | Both the H9 chamber and its H5 face have exact infinite-volume escapes |
| [Nonright H22 spherical profiles](fresh_h22_nonright_bridge_20261013.md) | Exact 6,480-profile census with 217 unbounded dihedral tails explicit; bounded numerical anchor search gives no survivor | Abstract profile completeness does not certify the numerical failures or exclude all nonright bridges |
| [Single-wall H22 cusp replacement](fresh_h22_core_replacement_20261013.md) | Delete one old wall and insert two perpendicular caps; exact cusp cut and actual finite vertex in an 813-facet chamber | Changes the fixed old core; separate radical replay passes, but a strict infinite-volume escape persists |
| [Coupled H9 hexagonal face](fresh_h9_h5_closure_next_20261013.md) | Two new caps give 125 facets and close an actual compact H2 hexagon | The containing H5 face has 53 open edges, up from 46; strict H5 and H9 escapes remain |
| [Full H22 replacement orbit](fresh_h22_replacement_orbit_20261013.md) | 12,331 actual facets and 5,760 certified finite vertices; exact compatibility for a whole parameter tail | A common escape excludes every acute wall in all 5,760 residual planes, with arbitrary real coefficients; an inherited cusp link also remains unbounded |
| [Complete replacement-vertex star](fresh_h22_replacement_edges_20261013.md) | All 22 edges are exactly open; two rational identities force every escape cutter to pair negatively with old walls 7 and 168 | All-label Schur bounds exclude a finite or ideal escape-cutting endpoint on any of these edges; preparatory additions and later vertices remain open |
| [Global H9 bank and opposite branches](fresh_h9_global_h5_20261013.md) | Exact common escape for all 38,209 original recipes; four further witnesses and 1,198 angle conflicts cover all 119 added branches | No finite-volume compatible selection from this dictionary, even after optional old-wall deletion; normals or seed must change |
| [Split-seed H9 scout](../scripts/fresh_h9_split_endpoint_walk_20261013.py) | Selected spherical extensions and extra angle equations; numerical proposals followed by exact reconstruction | Full stated branching scope, persistent edge queue, global compatibility and exact closure |

The earlier quartic search tested bounded powers of one explicit unit and square residual classes. Its successful caps were not certified nearest mirrors. The subsequent complete interval census includes nonsquare residuals and all bounded root solutions: it confirms the first two proposed endpoints and finds a nearer norm-(2+T) mirror before the third. At that stage, this established three complete local edge steps without exhausting the full chamber; the later global obstruction is described below. The split-seed scout still uses selected angle labels and anchors. Neither lane is currently a complete chamber implementation of Allcock's algorithm.

The [first explicit frontier audit](fresh_h9_split_edge_frontier_20261013.md) now enumerates all nine outgoing directions at each known simple finite vertex of the split-seed chamber. Its numerical traversal finds 4 closed and 32 open directed edges for the 72-wall chamber. Adding the exactly certified cap 2747 gives five finite vertices, 6 closed directions, and 39 open ones. Here an open direction reaches the light cone before another existing wall. The larger frontier records new edges exposed by a new vertex; it does not mean that adding a wall enlarges the chamber. This diagnostic does not enumerate missing mirrors or certify a complete lattice chamber.

The [next geometric frontier](fresh_h9_split_second_frontier_20261013.md) now has exact endpoint classifications for a selected 79-wall chamber: sixteen reached finite vertices, nineteen finite edges and 106 open ends. This uses exhaustive ninth-wall Gram checks, not numerical ray ordering. Its five closed finite face cycles do not eliminate the chamber's original escape.

The subsequent [global quartic certificate](fresh_h9_quartic_next_frontier_20261013.md) changes the arithmetic priority. Two integral involutions preserve local sectors connected by the complete first two edge steps. Their product preserves the same full chamber and has trace above the ambient dimension, proving infinite order and infinite chamber volume. No further wall enumeration in that exact lattice can yield a finite-volume reflection chamber. This is a full-host conclusion, unlike exhaustion of a bounded cap bank.

On the unrestricted geometric lane, [coupled remote insertions](fresh_h9_remote_caps_20261013.md) now retain the previous 79 walls and add six compatible caps, cutting two previous escapes and creating three actual finite vertices. The 85-wall chamber still has a strict escape. A separate [H19-to-H20 lift](fresh_dimension_bridge_20261013.md) creates another constructive lane with inherited cusp links and exact polar repairs. Further work must target the surviving gaps, retain metric parameters, and backtrack on incompatible choices. A persistent frontier remains useful, but local edge closure alone does not establish finite volume.

No new compact or finite-volume Coxeter dimension follows from this algorithm audit.

The [next construction round](new_host_first_mirror_positive_plane_bridges_2026_10_13.md)
also uses a different geometric move: two orthogonal caps in the positive
complement of a Lorentzian wall span. This creates a finite edge without
starting from a spherical edge extension. Both this H9 construction and
the targeted H20 extensions still have exact escapes. They supplement
edge traversal; they do not replace its global closure requirement.

The [following round](second_quartic_obstruction_coupled_caps_h22_orbit_2026_10_13.md)
closes the second arithmetic F4 host by the same complete-walk principle and
opens the A9 host above. On the H20 route, new parallel pairs reveal cusp
links of full quotient rank that are nevertheless unbounded because they
retain spherical components. Explicit highest-root walls complete those
components to affine diagrams with positive kernels. This repairs the two
links, but the 96-wall chamber's strict escape survives. It illustrates why
rank alone is not an adequate cusp-completion criterion.

The [next A9 walk](fresh_h9_a9_next_20261013.md) completes those eight
previously open central directions using four exact searches and reversal.
Every new endpoint exposes more directions: the ten known corners have
72 unresolved directions, or 36 reversal orbits. The central star is complete,
but the chamber frontier is not. The [94-wall geometric construction](fresh_h9_compact_next_20261013.md)
separately creates an actual compact H3 face with all six edges closed;
24 other directions at its four vertices remain open. The [H22 cusp-orbit construction](fresh_h22_cusp_orbits_20261013.md)
bounds over a million inherited links and still has an exact strict escape.
These examples distinguish local closure from the global volume requirement.

The [subsequent A9 walk](fresh_h9_a9_outer_20261013.md) closes that full host
without exhausting its outgoing queue. Eight additional complete first-mirror
calculations connect corners admitting an integral isometry with trace greater
than ten. The complete corner stabilizers and certified paths prove that the
map preserves the full reflection chamber, whose volume is therefore infinite.
The 92 reversal orbits remaining in the finite-prefix queue no longer warrant
further enumeration for this host.

The [following geometric checkpoint](a9_nonreflectivity_compact_prism_face_orbits_2026_10_13.md)
also sharpens the construction requirements. An actual compact H4 face does
not close its H9 chamber. In H20, a new facet with only 18 intersecting
neighbors cannot be repaired by any finite addition consisting solely of walls
dotted to that facet. A complete new H22 cap orbit admits every integer angle
label at least six, but one exact escape works throughout that family. These
are reasons to change the construction, not to enlarge the same wall banks.

The [latest checkpoint](changed_spatial_edges_full_orbit_obstructions_cusp_links_2026_10_13.md)
changes the arithmetic spatial lattice to D9 and carries out eleven complete
selected edge steps with its larger primitive norm pool. Certified relative
norm equations supply firstness over the quartic field; this is additional
work beyond the published integer-lattice theorem. The final step reaches a
corner that admits an integral oriented transport from another reached corner.
Its infinite order closes the full D9 host by nonreflectivity, without exploring
every remaining edge. On the geometric lane,
complete symmetry tests exclude the entire six-source H9 cap bank, while
eleven equivariant wall types repair all 241,920 specified new H22 cusp links.
Three compatible placements all retain exact global escapes. The distinction
between closing a local link and closing the whole chamber remains decisive.

The [following checkpoint](changed_core_trivalent_frontier_and_cusp_duals_2026_10_13.md)
adds exact trivalent endpoints outside the excluded H9 wall bank, then proves
a continuous all-label obstruction on their 16 open edges. All 940 small-label
spherical profiles and all 51 isolated dihedral tails are covered: a new wall
cannot both make a finite endpoint there and cut the specified remaining
escape while retaining the 104 walls. This is a scoped metric obstruction,
not a ban on intermediate noncutting walls or work in other components.
Separately, doubling the H19 source before a changed lift creates an escaped
179-facet H20 chamber, and quotient-space dual certificates exclude two
specified local H22 repair families from cutting their current escape.
