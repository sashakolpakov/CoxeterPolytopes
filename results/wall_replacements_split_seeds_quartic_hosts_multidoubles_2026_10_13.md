# Wall replacements, split seeds, a quartic host, and successive doubles

**No new compact or finite-volume Coxeter polytope in a new dimension has
been verified in this round.** Four constructive routes change the previous
starting configurations. Exact checks distinguish valid new mirrors and
angle repairs from successful finite-volume completion.

## Replacing the obstructed H9 core

The previous 46-wall configuration had a forced noncrystallographic cusp,
so adding walls while preserving it could not succeed. This search first
replaces walls. Removing old caps 31, 51, 197, 496, 1028, adaptive cap 19,
and either old cap 310 or 837 gives two different 39-wall configurations.
A complete exact connected-affine census finds no affine subdiagram in
either, and both retain the ten-wall seed. All wall pairings and a strict
interior are verified.

Targeted caps then cut actual escapes of these changed configurations.
The first two new walls give a 41-wall chamber. When the two-distance
polar construction misses its leading gap, the search changes the
constraint equations: it uses the first newly retained wall as an anchor
over `Q(phi,sqrt(132-60phi))`, prescribes a relation to polar wall 348,
and solves for three seed distances. Two of the resulting caps coexist,
giving 43 walls. A further polar coupling reaches 45 walls. Each extension
retains all preceding compatible walls.

The exact spacelike escapes remain, so these configurations still have
infinite volume. Absence of affine subdiagrams removes the old obstruction;
it does not establish compactness. The
[changed-core construction report](fresh_h9_replaced_core_20261013.md)
contains the final selected geometry, complete affine checks, exact escapes,
and independent replays. The final replay checks 1,980 pairs, 1,302 affine
inertias, and 24 strict escapes across two 45-wall configurations. A
[separate geometry audit](fresh_h9_replaced_core_independent_geometry_20261013.md)
also confirms every norm and the surviving strict escape of the preferred
configuration.

## Breaking the symmetry of a different H9 seed

The [split-leaf construction](fresh_h9_split_leaf_caps_20261013.md) starts
with sixteen walls, in two orthogonal sets of eight. Their cross pairings
are

```
(u_i,v_j) = -2 delta_ij - a - b epsilon_i epsilon_j,
a = 1/2+phi, b = -a/(12a+1),
epsilon = (1,1,1,1,-1,-1,-1,-1).
```

All cross pairs are strictly dotted. Breaking the leaf symmetry increases
the leaf span to Lorentz rank ten. A common strict timelike point and the
signature `(9,1,6 zero)` are exact.

Every real norm-two cap whose leaf angles are restricted to
`pi/2,pi/3,pi/5` is recovered by a complete 2,601-assignment census.
There are 192 such caps. Their compatibility graph has 8,478 maximal
cliques, of maximum size six; each has an exact strict spacelike escape.

The construction then leaves that finite-leaf restriction. The magnitude
`d=phi(1+sqrt(phi))` supplies 36 more exact caps with leaf pairings
`-1` and `-d`. The resulting **228-cap bank permits fourteen compatible
caps, hence an actual 30-wall Coxeter configuration**. Its complete graph
has 156,990 maximal cliques. Twelve exact rational-coefficient escapes
cover every one, excluding finite volume for every compatible selection
from this bank. Other caps are not excluded.

The [finite-leaf certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_exact.json)
and [combined-bank certificate](../data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_exact.json)
record all exact data. Finite-angle rejection in `Q(phi,sqrt(phi))`
uses the nonreal conjugates of elements with a nonzero radical coefficient,
so it is not limited to a numerically sampled angle alphabet.
An [independent Gram-inverse reconstruction](fresh_h9_split_leaf_caps_independent_20261013.md)
and a [separate Fraction/NetworkX replay](fresh_h9_split_leaf_caps_20261013_dotted_independent.md)
both confirm the complete graphs, clique censuses, and escape coverage.

The separate numerical scouts test 434 seed parameter points and a larger
bank of 3,040 numerical normal profiles at the selected seed. The latter gives 120 tested
configurations with 95–137 walls, all with numerical spacelike vertices.
Those are numerical failures, not exact exclusions of the entire larger
bank or a worldwide novelty claim.

The [subsequent exact reconstruction and endpoint proof](fresh_h9_split_leaf_endpoint_cover_20261013.md)
identify six duplicate profiles, leaving 3,034 distinct mirrors, and
exclude every finite-volume completion from that entire frozen bank.
This later proof is separate from the numerical scout described here.

## A different arithmetic ground field

The [conductor-sixteen host](fresh_h9_glued_host_20261013.md) uses

```
T = 2cos(pi/8), T^4-4T^2+2=0,
L = <-(T-1)> + F4 + F4 + <1>.
```

The coefficient has one positive and three negative real embeddings.
Consequently this is an admissible H9 form over a different
quartic field, rather than a scalar extension of the preceding quadratic
or cubic construction. Its ground-field isotropy is impossible, so a
finite-volume chamber made from its field-valued roots would be compact.

A complete bounded census produces fourteen walls through height 100.
The continuation frees all four coordinates of one F4 block together with
the scalar coordinate, checks both block placements through height 300,
and constructs **eight additional mutually compatible mirrors**. The
resulting 22-wall H9 cone and its 12-wall H5 face still have exact strict
spacelike escapes.

Independent enumeration reconstructs all 35 controls in those two slices,
all retained roots, and all 231 final H9 pairings. Every mirror in these
bounded slices compatible with the original fourteen walls is already
present in the final cone. This closes the declared bounded bank, not the
full host: roots using both F4 blocks, higher controls, and other
coefficient ideals remain open.

## H20 sections outside the old finite-star obstruction

The [successive-double construction](fresh_h20_multidouble_sections_20261013.md)
changes the ambient H21 chamber. Its first seams include parallel mirrors,
so their reflection subgroup is infinite. The resulting finite unions
therefore lie outside the previously excluded finite spherical stars.

Exact convexity checks give unions of 4, 8, 16 and 32 original Borcherds
chambers, with respectively 744, 1,456, 2,625 and 3,353 genuine boundary
facets. Six new actual right-angle repairs occur in the eight-chamber
sections. Further doublings remove the displayed offending seams but
create actual obtuse ridges in the tested family.

The next construction changes the cusp quotient direction itself. Four
previously unmaterialized directions give eight further finite-volume
H20 sections, all with complete rank-19 Coxeter cusp links. Six prescribed
angle equations are certified actual pi/3 or pi/5 ridges. Every section
also has a different actual forbidden ridge.

An exact interior-point construction solves the cusp linear equations and
then shifts toward the cusp. It succeeds where projection of the old
center fails. Both independent section replays pass. The failures concern
the literal sections and repairs checked here, not every section of these
larger unions or every chamber assembly.

## What remains open

The strongest changed H9 core still requires another compatible cutter
and complete finite-volume geometry. The quartic arithmetic host needs
roots beyond its exhausted two-slice bank. The split seed has a larger
numerical cap bank whose compatible completions remain undecided beyond
the tested selections. The H20 construction leaves other section
parameters and chamber assemblies available.

None of the exact finite-bank exclusions is a dimension-wide
nonexistence theorem. A positive result still requires a complete Gram,
actual facets, a strict interior, and an exact compactness or finite-volume
certificate, followed by an independent audit and a literature comparison.
