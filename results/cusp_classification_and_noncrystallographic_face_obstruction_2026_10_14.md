# Cusp classification and noncrystallographic face obstruction

No new compact or finite-volume Coxeter polytope, and no new dimension,
has been established. This checkpoint replaces several sampled section
searches with uniform obstructions and changes the compact-H9 incidence
search beyond the previously excluded families.

## Literature constraints

The full texts were rechecked following the user's warning. Esselmann
excludes noncompact arithmetic reflection groups in dimension20 and every
dimension above21. Borcherds attains21. The bound does not exclude compact
arithmetic examples in22–29; Vinberg's general arithmetic and compact
bound is29. See the [literature audit](literature_frontier_recheck_20261013.md).
The July2026 primary literature still records compact examples only
through8 and noncompact examples through21. A nonarithmetic H20 example
would also fill a missing dimension in the known list.

Allcock's published algorithm operates on integral Lorentzian lattices.
Its edge-traversal organization can guide a nonarithmetic geometric
search, but does not supply metric realization or global closure for an
abstract diagram. No excluded high-dimensional rational reflective-host
enumeration was launched in this checkpoint.

## Every H3/H4 face in a Borcherds H20 section is excluded

The [coupled H3 proof](fresh_h20_coupled_h3_obstruction_20261014.md),
together with the preceding one-contact A3 proof, excludes every actual
H3 face in any hypothetical finite-volume Coxeter H20 section of the
saved210-wall source. It therefore excludes H4 faces as well. This does
not exclude sections without such faces.

The proof classifies the contacts on the three H3 walls, allowing
arbitrary contacts on every other source wall. Two-contact joins,
two-contact connected paths and three-contact paths all force unchanged
integral boundary roots on a compact face, while exact irrationality
prevents an H4 extension. That compact rational arithmetic face is
impossible by Meyer isotropy. The arithmetic check proves30 polynomial
irreducibility statements, two joined-field norm obstructions and a
nonzero D4 extension defect. A
[separate review](fresh_h20_coupled_h3_obstruction_20261014_review.md)
rederives the algebra and checks the local completeness argument.

## All-retained single-factor cusp sections are excluded

The [complete classification](fresh_h20_single_factor_classification_20261014.md)
reduces arbitrary real normal directions to rank-one graph cuts in all53
affine factors at the18 source cusp representatives. There are287 cuts:
119 fail a local Coxeter-angle condition and168 give valid product links.
Each of the168 surviving normal planes has an actual obtuse global ridge
for every interior shift. Seventy-seven were previous leaf-sum planes;
91 are additional directions.

The [independent SymPy checker](fresh_h20_single_factor_classification_20261014_root_review.md)
reconstructs the local Grams and all182 continuous ridge certificates,
including353 distinct rational endpoint points checked against all210
source walls. Fourteen planes require two overlapping half-interval
proofs. The theorem assumes all original facets of the active factor
remain distinct actual local facets. Higher-rank sections omitting one
local facet remain open.

## Dihedral edge sections

The [uniform dihedral proof](fresh_h20_dihedral_edge_sections_20261014.md)
excludes every one-contact A2-to-I2(m) section for all integers m>=7 and
all real positions on every source edge. An exact Gram table reduces
the only possible exceptional contact to
cos(2pi/m)cos(2pi/n)=1/2, whose sole integer solution here is m=n=8.
That case also forces an impossible compact rational face. The saved
source census has5,048 anchored families. The
[review](fresh_h20_dihedral_edge_sections_20261014_review.md) distinguishes
proof/source review and replay from independent arithmetic implementation.

The [label5 corollary](fresh_h20_dihedral_edge_sections_20261014_m5.md)
uses the new all-H3 obstruction to close both possible orientations of
an H3 extension. Labels4 and6 are not uniformly excluded.

For [label6](fresh_h20_dihedral_edge_sections_20261014_g2.md), a necessary additional-facet condition reduces each edge
family to finitely many candidate sections. The broad scout found481,082
numerically distinct candidates:481,016 proposed ridge failures and66
boundary remainders. Exact arithmetic identifies all66 remainders as
old-endpoint tangencies. Complete exact checks on three selected anchored
families reject131 sections,105 by obtuse ridges and26 by angles strictly
between pi/3 and pi/2. The remaining broad numerical rejections have not
been converted into an exact full-family theorem.

## Compact-H9 incidence changes

The [exact incidence-family checker](fresh_h9_incidence_repair_20261014_verify.json)
excludes arbitrary additional right-to-dotted mutations of72 earlier
15-facet patterns: their nerve would be a flag8-sphere with only15
vertices, whereas at least18 are required.

Adding finite contacts at originally right pairs requires at least seven
changes to hit the relevant surviving affine circuits, with eleven
minimum supports. All1,408 label3/4 assignments on those supports have
empty potential-vertex sets after exact unique-ridge peeling. Arbitrary
dotted distances and further right-to-dotted changes cannot repair them.
The [root review](fresh_h9_incidence_repair_20261014_review.md) includes a
fresh exact replay and mathematical/source review.

Separate searches that also release old parallel incidences produce
closed combinatorial proposals with216 potential vertices. Four bounded
proposal spaces were then rejected using55 exact modular rank obstructions
and145 numerical positive-definiteness failures; no model reached the
dotted-distance solver. This mixed exact/numerical scout is not an exact
exclusion theorem for the broader construction. See the
[full incidence report](fresh_h9_incidence_repair_20261014.md).

## Remaining constructive scope

For H20 sections, avoid H3/H4 faces and the excluded all-retained cusp
directions. One-omitted-facet cusp links, other crystallographic local
profiles, and changed ambient chambers remain open. The one-contact
label6 necessary candidate list offers a finite exact follow-up, but
the present numerical census is not a nonexistence proof.

For compact H9, further changes must leave the precisely described
flag or minimum-contact families. The four bounded216-vertex proposal
spaces did not produce a metric realization; broader incidence changes
remain open. Existing H22 constructions still have certified
infinite-volume escapes.

## Frozen verification

The H3 producer and separate identity reconstruction pass saved-result
replays. The all-retained cusp producer passes its replay, and independent
checks verify all353 endpoint points and the28 delicate midpoint recipes.
The dihedral theorem passes source/proof review and replay; the label5
corollary passes a fresh root replay. Both exact label6 certificates pass
their saved replays. The compact-H9 exact-family checker and the saved
combinatorial/modular-rank checker pass, with numerical limitations retained
explicitly in the report.

The final audit resolves65 local links in12 reports and checks62 source
hashes across13 authoritative manifests. No search job remains running
from this checkpoint. The goal of finding a new finite-volume example
remains open.
