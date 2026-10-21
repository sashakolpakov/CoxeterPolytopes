# Changed arithmetic host, trivalent frontier, replaced H20 core and cusp duals

**No new compact or finite-volume Coxeter example or dimension has been
established.** This continuation produces two exact geometric constructions
outside previously excluded choices, closes another full arithmetic lattice,
and proves continuous obstructions that identify what must change next.

## Arithmetic H9: the tripled D9 lattice is closed

For `K=Q(T)`, `T=2cos(pi/8)`, the form
`<−3(T−1)> + D9` has a different determinant squareclass from the preceding
D9 form. Its calculations cannot be inherited from the old lattice.

Four new complete first-edge censuses use all 26 safe primitive norm classes
dividing 24. They examine **136 finite profiles and 1,161 interval points**,
joining two full `B3+B6` corner sectors to the same initial chamber. A newly
constructed oriented transport preserves the actual D9 lattice in both
directions and has trace greater than ten. Hence the full chamber has an
infinite-order integral stabilizer and **infinite volume**. This is a
full-lattice obstruction, not just an escape from a finite wall prefix.

The [host report](fresh_h9_tripled_d9_20261013.md) and
[independent transport review](fresh_h9_tripled_d9_20261013_transport_review.md)
record the exact paths, lattice tests and explicit dependence on certified
PARI norm/unit computations. Separate Fraction implementations replay the
edge censuses and the transport.

A [parameterized matrix factorization](fresh_h9_tripled_d9_20261013_mechanism.md)
explains the repeated D9 transport under stated first-corner and congruence
hypotheses; it does not prove every coefficient fails. The next candidate
`<−3(T−1)>+E8+<1>` has a different `E8+A1` zero chamber and ten safe norm
classes dividing six. [Exact structural checks](fresh_h9_tripled_d9_20261013_e8_triage.md)
show the old transport is not integral there and the mirror arrangements
differ in the chosen coordinates. An unidentified conjugacy of the full
indefinite groups remains possible; no finite-volume or full-group novelty
claim is made for this candidate.

## Geometric H9: two exact trivalent endpoints, then an all-angle screen

Two new walls preserve all 102 preceding facets and cut both saved escapes.
One cuts the common escape of the entire old 13,824-image family, so this is
a construction outside that excluded family. The two walls close two old
open edges at actual finite vertices. Their 205 new pairs consist of ten
right angles, seven pi/3 angles, one pi/5 angle and 187 dotted pairs.

The resulting **104-facet chamber still has an exact strict spacelike escape**.
The [construction](fresh_h9_component_bridge_20261013.md) checks all new pairs,
the actual vertices, every relevant affine subset and the remaining escape.
The [root review](fresh_h9_component_bridge_20261013_root_review.md) combines a
full exact replay with an independent 110-digit linear-system reconstruction.

The two new corners expose 16 open edges. The
[continuous frontier proof](fresh_h9_component_bridge_20261013_continuous.md)
enumerates 940 spherical extension profiles with labels 2,3,4,5 and both
norm-hyperbola branches. Necessary rational quadratic intervals exclude
1,876 branches. On the four remaining full intervals, pairing with seed6
lies strictly between -1 and 0, which forbids every Coxeter label.

The accompanying
[dihedral-tail certificate](../data/enumerations/fresh_h9_component_bridge_20261013_dihedral_verify.json)
checks all 51 isolated `I2(m)`, `m>=6`, cases by positive two-wall dual
identities, even after discarding the norm constraint and angle discreteness.
Together these cover **all labels and all real positions**: no compatible
new wall can simultaneously retain all 104 walls, create a finite endpoint
on one of these 16 edges, and cut the fixed escape. This does not exclude
intermediate walls that do not cut it, other frontier components, ideal
endpoints, or a change of retained walls.

The [combined independent review](fresh_h9_component_bridge_20261013_continuous_review.md)
replays both exact proofs, independently enumerates all 1,048,576 finite
edge/profile assignments, reproduces the 940 admissible profiles, and checks
the two-wall dual cancellation used for all 51 dihedral tails.

## H20: replace the source core before lifting

Doubling the known finite-volume H19 source across an even wall gives 49
source facets. A different lift phase produces a 71-wall H20 core; retaining
102 compatible later walls and adding six exactly checked rational polars
gives **179 actual facets**. The old deficient facet is removed, so the
previous additive obstruction does not apply to this construction.

The new core cuts the preceding escape by `12171/4096>0`. All 10,780 changed
or new pairs pass; the other 5,151 pairs belong to the explicitly bound old
subset. The final strict escape has square
`10640902782393308787/20000000000000000000>0`, proving **infinite volume**.

An exact census of all 204 inherited source cusp rays checks every current
wall inequality. Seventeen links have only height-zero incident normals,
affine quotient rank 18 rather than 19, and an explicit free vertical
direction. This is a concrete repair target, not a census of every cusp.

The [construction and cusp report](fresh_h20_changed_assembly_20261013.md)
records both proofs. The [root review](fresh_h20_changed_assembly_20261013_root_review.md)
replays the complete certificate and independently reconstructs the rational
double, lift, six polars and escape. Retained radical normals are an explicit
shared dependency. No finite-volume conclusion is inferred from the source's
finite volume.

## H22: continuous dual certificates at the exposed cusp

At the lower placement's cusp `delta4`, a new proposal joins the five finite
chains into one affine A20 cycle while retaining the old affine A1. This
changes the bases used in the previous excluded repair family.

For each of the **32** possible chain-endpoint bridge types, an exact
certificate using at most two old walls proves that every compatible bridge
misses the fixed global escape. The proof covers every real shift and every
direction in the four unused tangent dimensions, assuming only acuteness
against the original 210 walls. An independent Fraction implementation
reconstructs all 32 quotient-space Cauchy inequalities. See the
[coupled repair report](fresh_h22_coupled_cusp_20261013.md) and
[independent audit](fresh_h22_coupled_cusp_20261013_independent.md).

A separate [five-wall conic identity](fresh_h22_incident_support_20261013.md)
handles the all-right incident profile: any normal perpendicular to `delta4`
and its 17 old incident walls also misses the escape, with no norm, field or
height restriction. Its independent checker verifies all 23 coordinate
identities and positive coefficients. Neither result excludes other incident
diagrams, walls away from this cusp, or repairs transported to other cusps.

## Consequences for the algorithms

The [Vinberg–Allcock audit](vinberg_allcock_search_strategy_20261013.md) remains
the relevant distinction: distance-ordered mirror enumeration versus edge
traversal. Allcock's published method is itself a lattice algorithm over Z,
combining corner diagrams with two-dimensional Lorentzian arithmetic. The
number-field first-edge proofs here supply additional certified norm-equation
work. The unrestricted construction route uses diagram extensions together
with continuous metric constraints and separate global checks.

This round shows both roles. Complete arithmetic edge paths can close a
lattice negatively without exhausting its whole edge queue. Combinatorial
extensions can propose genuinely different walls, while continuous interval
and dual proofs explain why certain complete local families cannot cut a
specific global escape. Closing a corner, edge or cusp link still does not
establish finite volume for the whole chamber.
