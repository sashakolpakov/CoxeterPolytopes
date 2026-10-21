# Finite-edge triangles, forced bases and a better exact boundary

All included packages, including the coefficient-12 theorem and the three
new construction gates, are final and frozen. This checkpoint continues
the [previous archive](prism_base_capacities_and_boundary_replacement_2026_10_17.md)
without changing its frozen source, top-artifact or report sections.

This tranche supplies stronger finite-volume triangle accounting, closes
the union-of-constraints gap in the rank-13 arithmetic weighting diagnostic,
and improves an exact dimension-nine chamber boundary. It establishes no
new finite-volume chamber or ambient dimension ceiling. Existing dimension
proofs retain their previous scopes and validation status.

The unrestricted finite-sided, finite-volume dimension bound **163 is a
complete, internally checked computer-assisted proof with no identified
mathematical gap**. Its reviews are internal, not external refereeing, and
it is not formalized in a proof assistant; the
[proof-status note](dimension163_proof_status_20261017.md) makes that
distinction explicit. The repository's compact arithmetic bound21 is
separate from the published noncompact rational arithmetic case. Neither
scope should be substituted for the other, and this checkpoint changes
neither bound.

## Finite edges and ideal vertices

The [finite-ended-port theorem](fresh_finite_ended_simplex_ports_20261017.md)
assigns every edge with two finite endpoints to a unique maximal actual
simplex face having at least two finite vertices. Include isolated finite
edges as dimension-one blocks. At a finite vertex, two triangles sharing
a finite edge force a tetrahedral three-face by three-connectivity. A
larger clique closes to a simplex even when some displayed vertices are
ideal: deleting those ideal vertices cannot disconnect the remaining
finite vertices. This does not assert disjointness for arbitrary
ideal-ended ports.

The complementary [one-ideal cusp theorem](fresh_cusp_one_ideal_simplex_capacity_20261017.md)
works inside each simplex-factor fiber of a product cusp link. Within a
fiber, maximal simplexes whose sole ideal vertex is the cusp have disjoint
port sets. A port may occur once in each factor, so the factor multiplicity
cannot be discarded. The exact square pyramid saturates this local bound.

The [multi-ideal analysis](fresh_multi_ideal_triangle_scope_20261017.md)
keeps triangles with two or three ideal vertices separate. Its original-wall
star argument bounds the dimension of any actual simplex face by nine,
using a minimal indefinite subdiagram containing all of its black facet
roots. It does not assume that induced simplex-face angles are Coxeter
angles. Compact ambient chambers give the smaller Lannér case; a compact
face inside a noncompact chamber can still involve a quasi-Lannér diagram.

These statements yield actual incidence constraints, not a completed cusp
payment or an ambient dimension exclusion.

## Classical diagram census and finite-edge triangle coefficients

The [exact high-rank census](fresh_high_rank_quasilanner_census_20261017.md)
reproduces the classical quasi-Lannér diagrams of orders six through ten:
146,262 nonzero attachments to 24 affine source diagrams give 84 retained
attachments and 26 nonisomorphic diagrams. The counts by order are
`12, 3, 4, 4, 3`. Every proper deletion is checked positive semidefinite;
the complete finite/affine classifications, compact order-five bound and
general order-ten ceiling are explicit published inputs. This is a
reproduction with exact deletion data, not a new classification.

The [first finite-edge tradeoff](fresh_finite_edge_cusp_tradeoff_20261017.md)
uses the order-ten diagrams and additional-white-root restrictions to prove

\[
 2(T_0+T_1)\le14E_{ff},
\]

where T_s counts triangles with exactly s ideal vertices and E_ff counts
finite-to-finite edges. The ideal counts of an arbitrary induced simplex
are not simply read from the isolated Coxeter diagram: its supporting
white roots must also be controlled.

The [coefficient-12 successor](fresh_finite_edge_cusp_tradeoff12_20261017.md)
extends that argument to eight-dimensional simplex faces. It reconstructs
the relevant seven order-nine/ten diagrams using a separate Fraction-based
`Q(sqrt(2))` arithmetic engine, treats all 34 black/internal-white choices,
and handles the remaining affine-E8 black block with its full external
white attachments. Its proved conclusion is

\[
 2(T_0+T_1)\le12E_{ff}
 =6n v_f-6E_{fi}\le6n v_f.
\]

The relaxed extremizer is a seven-simplex with two finite and six ideal
vertices; no Coxeter realization or sharpness claim is inferred from that
numerical extremum. Triangles with two or three ideal vertices and the
quadrilateral demand remain unbounded by this finite-edge statement. In
particular, it gives no dimension ceiling for all-ideal chambers.

## Forced simplex bases and the combined arithmetic obstruction

The [mixed/shared certificate](fresh_sqrt5_mixed_shared_transfers_20261017.md)
first closes the previous scope gap between mixed three-port indicators and
represented prism/base transfers. Their joint exact lower bound is about
12.5349338525, above the dimension-13 target 12.

The [forced-base domain](fresh_sqrt5_forced_base_simplex_domain_20261017.md)
then adds every missing base type appearing in any provenance of a positive
quad row. There are 4,195 such stars: 3,460 maximal edges, 447 triangles,
205 tetrahedra and 83 four-simplexes. Their 9,496 finite corners use 694
already represented vertex forests, so no new ordinary pair variables or
forest budgets are needed. The patch adds 1,224 one-port variables and
2,848 shared-transfer columns.

A maximal edge is a singleton block in the disjoint-port theorem. Its row
has zero curvature demand, no pair-corner coefficients, and its two
endpoint one-port charges. Those charges are shared with existing simplex
rows in the same complete vertex budget. A square prism uses all four
edges as bases; no arbitrary opposite pair is selected.

The [exact analytic repair](fresh_sqrt5_forced_base_dual_repair_20261017.md)
retains every old quad and represented-simplex mass, adds the required
missing-base masses, and checks all 211,339 positive-cost variables and
5,693 signed zero-cost transfer columns. Its lower bound is

\[
 \frac{31154121494284588800}{2547155390157199561}
 \approx12.2309465746>12.
\]

The [full mixed/forced-base union](fresh_sqrt5_mixed_forced_base_transfers_20261017.md)
also satisfies every mixed three-port indicator with **no further increase
in vertex mass**. It checks 56,045 triple orbits on all 752 positive forests
and reconstructs all 4,195 new stars. Edges and triangles have zero mixed
indicator at each corner; the larger simplex footprints are included.
No new optimizer run was needed for either certificate.

Every merged quad provenance remains available. The transfer coefficients
use coordinatewise minima valid for all provenances; the certificate does
not assert a simultaneous realization of those minima by an actual split
of a merged row mass. This is a partial weighting-family obstruction,
not a compact realization or a dimension theorem. Other missing bases,
fixed arithmetic-form classes, unenumerated face types and 231 untouched
finite forests remain outside. Old zero-mass mixed features also need
reconstruction before these exports could be used as a complete primal
optimization matrix.

## Exact dimension-nine replacement

The [replacement47 chamber](fresh_nonarithmetic_replace47_close12_20261017.md)
has 64 actual facets and a complete ambient census of 198 finite vertices.
It replaces the sole old finite vertex incident to wall47 by two vertices.
All previously certified compact faces of dimension at least two persist.
The whole 12-face becomes a compact quadrilateral and the compact 8-face
remains a heptagon.

The [complete intrinsic S6 boundary](fresh_nonarithmetic_replace47_close12_boundary_20261017.md)
has 25 finite vertices, seven exterior rays, no ideal vertices, 32 compact
geodesic edges and 15 open ends. This improves the previous 16 open ends
while recording the newly opened 7–30 edge. The census exhausts all 816
projected side triples and every resulting edge. It is complete intrinsic
coverage, distinct from an ambient exterior-cone census.

A strict ambient escape survives, so the chamber still has infinite volume.
The replacement and its boundary both passed generation and fresh replay;
their physical arithmetic is inherited from the frozen exact construction
engine, without a separate independent arithmetic review of this package.

## Construction obstructions and bounded scouts

The four [region47 scouts](fresh_nonarithmetic_release54_region47_20261017.md)
are preserved as bounded numerical discovery. The initial screen tests 680
branches; the paired endpoint screen tests 292,896 cross-group pairs; the
third-side screen tests 100,548 branches; the relaxed replacement screen
tests 2,150 real branches. Their promoted simplest replacement is verified
exactly in the preceding package. Unexamined labels and continuous domains
are not excluded by those bounded banks.

The [single-cap gate](fresh_nonarithmetic_replace47_single_cap_gate_20261017.md)
excludes one wall making all four remaining endpoints of the targeted
boundary component finite, for every Coxeter label and coefficient field.
The H4-squared endpoint forces the necessary orthogonality; an exact norm
bound on a relaxed continuous rectangle is strictly below `219/125 < 2`.

The [two-edge gate](fresh_nonarithmetic_replace47_edges1252_3052_gate_20261017.md)
excludes one wall simultaneously terminating edges 12–52 and 30–52 while
retaining walls4 and7, even when47 is removed. All 14 finite profiles are
covered, including every nonzero A2 attachment. Four fail endpoint
positivity; all 20 real branches from the others are strictly obtuse to
one of the retained blockers.

The [paired-endpoint gate](fresh_nonarithmetic_replace47_b_pair_gate_20261017.md)
covers two right-S6 caps with specified duties on the7-side and52-side
endpoints. Thirty low-label norm branches leave three left and eleven right
caps, whose 33 mutual products all fail the required finite Coxeter contact.
Two infinite label tails are excluded by exact increasing norm polynomials.
Its preceding numerical pair scout remains a bounded diagnostic. Nonzero
A2 profiles on the right cap, different endpoint duties, larger cap
clusters and retained-wall replacements remain outside this paired theorem.

All three gates passed generation and fresh saved replay. They are
restrictions on stated cap duties and retained-wall families, not exclusions
of all coordinated replacements or finite-volume completions. No new wall
was appended in these obstruction packages.

## Archive status

The [auditor](../scripts/audit_finite_edge_forced_bases_checkpoint_20261017.py)
preserves all three frozen predecessor sections—sources, top artifacts and
reports—while allowing README and RESEARCH_PLAN to be updated and bound
as new entrypoints. It recursively checks the new source manifests,
record hashes, local links and frozen inputs.

The [audit output](finite_edge_triangles_forced_bases_and_boundary_2026_10_17_audit.json)
certifies archive integrity only. Exact arithmetic, inherited inputs,
written geometric proofs, internal reviews and bounded diagnostics retain
their separate scopes.
