# Independent full-source review of the compact RA4 W>=16 manuscript

19 October 2026. Reviewer: finite57 research lane. **PASS within the stated
external geometric inputs.** This is a full written mathematical and source
review, not an execution, PDF-layout review, formal proof or priority survey.
No source edit or scientific, Lean or TeX job ran in this lane.

The reviewed primary source is
`manuscripts/compact_ra4_volume_constraints.tex`, all1073 lines,48233bytes,
SHA256
`20e5b4a3aff508520b7db3c4ae1a8904c0aa2d3299b4164a657bb6432433c49e`.
Its preparation/proof map is
`docs/compact_ra4_volume_constraints_preparation_20261019.md`, SHA256
`f5dd0c07e0881fc12bad110b38344475c6e24ae355e1a04c0d60b530a519cf7b`.
Both were read completely; hashes were checked by streaming. The source
and its reviewed result packets were left unchanged.

## Relation to the earlier reviews

The prior full-source review
`results/review_compact_ra4_volume_constraints_finite57_20261019.md`
binds the preserved900-line W>=15 source at9fe28180…d47738. That review
remains a review of those historical bytes. The present review separately
binds the W>=16 successor and includes the revised Section7 closure lemma,
the complete new Section9 endpoint proof, the changed statement/scope,
the literature credit and all four figure sources.

The underlying endpoint packet
`results/compact_ra4_fifteen_excess_equality_exclusion_research_20261019.md`
remains at SHA256
`237d9ccc7c91062b268124fc78c4fa82dbcf4bbf6b0b6b32c8eae013c713ff87`.
Its earlier independent review is
`results/review_compact_ra4_fifteen_excess_equality_exclusion_finite57_20261019.md`.
The complete manuscript now contains its argument rather than depending
on an omitted enumeration or referring readers to an agent message.

## Geometric and counting interfaces

Sections2--4 retain complete proofs of the compact normal gap, actual
intersection of orthogonal walls, common-neighbor clique and opposite-edge
distinctness, the three-dimensional polygon capacity/equality case, and
the four-dimensional singleton-footprint count. In particular a formal
normal intersection is not silently treated as an actual face: the
projection preserves every other wall inequality, using nonpositive wall
products. Compactness excludes parallel walls at a common ideal endpoint.
The dual-ear argument supplies the stated Loebell equality lattice; the
later metric step separately invokes three-dimensional uniqueness.

Section5 uses the exact incidences and Euler identities consistently:
`d(F)=(V(F)-20)/2`, `sum_F d(F)=2W`, `W=v-5f`, and
`Vol(P)=pi^2(f+W+16)/12`. The neighbor sum uses distinct facets, without
requiring their ridge contributions to be disjoint. Combining it with
`f>=2V(F)+5` gives the46/46/96 local vertex cutoffs with the correct strict
and nonstrict volume thresholds and parity. The120-cell and its double
are comparison examples, not assumed global extrema.

Section6 makes the external metric inputs explicit: recognition of a
12-face compact RA3 polyhedron as a regular dodecahedron, and compact
nonobtuse RA3 uniqueness. The whole-facet propagation for W=0 uses actual
regular dodecahedral facets, then connectedness; it does not infer a
universal120-cell decomposition from a local count.

## Revised closure and the W=15 endpoint

The revised `lem:three-edges` in Section7 is valid with **three** fixed
consecutive edges. Their oriented right-angle chain fixes the two next
supporting lines as the remaining sides of the regular pentagon. These
are actual polygon walls and are orthogonal. Orthogonal projection gives
their intersection as an actual finite vertex, with all other inequalities
strict, so an additional wall cannot truncate this closing point. This
last actual-wall step is essential and is now present. A zero-excess
third facet gives an entire regular-pentagon edge of length L, so the
carrier application has the required metric hypothesis.

The G/H distinction is maintained throughout. G records positive ridge
excess and H every adjacency between positive facets. An adjacency in H
is an actual edge/vertex incidence in the appropriate ridge by the
orthogonal-wall lemma. The forbidden induced square is excluded in H,
not in an uncompleted subgraph of G. The eligibility proof and four-node
graph lemma therefore give minimum H-degree5 whenever minimum G-degree
is at least3, including the W=15 equality case.

Section8 gives a complete written nine-node complement argument. In the
triangle case its three external degree slots and the forbidden induced
2K2 give the required independent five. In the trianglefree nonbipartite
case the shortest odd cycle is C5; zero-contact outside vertices are
isolated, and the last independent pair is chosen after deleting the
two nonadjacent contacts. The N<=8 cases are stated separately. No graph
census is imported.

I checked every case of Section9 against the reviewed endpoint packet:

* The degree-two G case at W=15 forces k=6 and equality in the entire
  distinct-facet sum. Exactly15 positive facets have excess2, hence each
  is a Loebell L6; each has14 positive neighbors. This contradicts the
  normal-space K5 obstruction.
* Otherwise the G minimum degree is3. The already proved small-graph
  argument forces N=10; the inequalities force E=15, cubic G and all
  weights1. Thus every positive facet has three hexagons and12 pentagons.
  If H had minimum degree6, its subcubic ten-node complement would give
  an independent five on any nine-node induced subgraph. Hence an actual
  positive facet F has exactly five positive neighbors.
* The pentagon lemma uses the unique common perpendicular of the two
  prescribed outer side lines. The hexagon lemma first excludes
  nonopposite exceptional positions by the three-edge closure, then
  explicitly derives the two right-angled hexagon equations. They force
  both exceptional sides to have length2L. The signs in the Lorentz
  normal-product calculation and the identity
  `phi^2/(phi-1)=phi^3=cosh(2L)` agree.
* The five-node graph J has exactly three marked hexagons. A marked node
  has degree>=2; at degree2 its two neighbors are nonadjacent and its
  two carrier edges are opposite and have length2L. At degree3 a
  three-node neighbor path is impossible: it makes the three remaining
  zero-carried edges consecutive. A neighbor triangle gives K4 in J.
* The universal-node case reduces to a four-node forest and allows at
  most two marks. The triangle cases distinguish absence/presence of the
  outside edge; in the latter case the two triangle-contact sets are
  both disjoint (maximum degree3) and nested (no induced C4), exhausting
  their possibilities. The trianglefree cases are a forest or C5.
  Three possible marks in the forest force P5. Its endpoint pentagon,
  and a pentagon next to a marked vertex in C5, have enough L sides to
  be regular, contradicting the shared2L side. Thus all cases close.

The warning that the abstract ten-node conditions alone are consistent
is appropriate: the proof does not claim an invalid graph-only exclusion
of C5[K2]. The resulting theorem is exactly W=0 or W>=16, together with
the separately proved zero-case recognition.

## Attribution, figures and verification boundary

I checked the changed pointers against the published primary PDFs.
Kellerhals--Perren, Proposition3.2(d), gives the cited volume formula;
its proof of3.2(b), printed p1311, contains the outgoing-edge/distinct
opposite-facet mechanism now explicitly credited. Inoue, Lemma3.3,
pp1529--1530, contains the3D excess identity now credited. These are
antecedents, not new repository results. The source does not turn a
limited literature comparison into a priority claim.
[Kellerhals--Perren](https://homeweb.unifr.ch/kellerha/pub/YEUJC1817.pdf#page=14),
[Inoue](https://msp.org/agt/2008/8-3/agt-v8-n3-p12-p.pdf#page=7).

All four TikZ figures were read as source. The dual-ear and incidence
figures keep their combinatorial roles. The chain figure now shows three
fixed L edges and the two prescribed closing walls, with an explicit
boundary-order rather than Euclidean-metric caption. In the nine-node
figure the U--V edge is retained; its highlighted set is independent,
and the example does not exhibit an unintended induced2K2. The figures
illustrate the written arguments and carry no independent proof claim.

No new issue or edit is requested. This PASS concerns the complete human
argument under its named acute-chamber, RA3 metric-uniqueness and
orbifold-volume inputs. It does not certify PDF rendering, end-to-end
Lean geometry, worldwide novelty, a minimum/second-volume theorem,
a universal120-cell decomposition, or a new dimension ceiling. The
updated source has not been built in this lane.
