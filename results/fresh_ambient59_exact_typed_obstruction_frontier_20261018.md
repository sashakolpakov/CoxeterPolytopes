# Exactly eight actual finite forests remain above ambient59 curvature

The complete current payment rule has only **eight actual formal finite
forests** of rank59 with cost strictly greater than783/2. They occur
in two rank patterns. The previous third representative row mixed
incompatible type allowances: its true maximum is781/2, already below
curvature. The global finite maximum remains787/2, and the completed
unrestricted dimension bound remains59.

The [producer](../scripts/fresh_ambient59_exact_typed_obstruction_frontier_20261018.py)
and [certificate](fresh_ambient59_exact_typed_obstruction_frontier_20261018.json)
refine the frozen [complete A14/C/D16 certificate](fresh_a14_cd16_ambient59_frontier_20261018.md).
They do not change any local weight or coupon. The completed geometry,
complete Bellman bounds and positive ideal-surplus theorem are inputs.

## Exhaustive reduction to three small type searches

The frozen finite certificate covers961 fixed multiplicity vectors.
Expanding every rank5 count gives2574 vectors. Its current integer
sixteenth tables include the explicit rank15 interaction coordinate
and both contextual payments. This new checker reconstructs that
complete expanded domain and matches every saved vector bijectively.
Exactly three expanded vectors have representative upper cost greater
than6264 sixteenths, or783/2:

|Fixed components by rank|Residual rank|Representative upper|
|---|---:|---:|
|three rank15 and two rank5|4|392|
|one rank14 and three rank15|0|787/2|
|two rank14, one rank15 and one rank16|0|1569/4|

All other actual forests are already at most curvature by the frozen
simultaneous type-domination and complete Bellman certificate. This
new checker imports that proved bound; it does not rerun its millions
of last-item inequalities.

For the three remaining vectors, enumerate every actual finite type
of each fixed component. At ranks14--16 these are A, B and D; rank5
also has precisely A, B and D. For the residual rank4, enumerate all
five partitions4,3+1,2+2,2+1+1,1+1+1+1 and every actual type in them.
This includes A4/B4/D4/F4/H4, A3/B3/H3 and the inherited rank2 classes
A2/B2/I2*. The I2* class represents all other rank2 dihedral types,
whose current weights and coupon behavior are identical. Thus there
is no finite-label cutoff in this refinement.

Permutations of identical components are removed. The three searches
contain1080,30 and54 actual typed forests respectively,1164 in total.
Every forest has rank59 and its cost is computed with the exact actual
base, all ordinary cross-component payments, and both contextual
costs. The source checks each cost against its old representative
upper bound as an additional consistency check.

The first vector's maximum is781/2=390.5, including every residual
partition, so this entire vector is removed from the obstruction
list. The other two maxima are attained and remain unchanged.

## The eight actual types

|Actual whole-component forest|Cost|Excess over391.5|
|---|---:|---:|
|A14 + 3A15|392|1/2|
|D14 + 3A15|393.5|2|
|D14 + 2A15 + D15|393|3/2|
|D14 + A15 + 2D15|392.5|1|
|D14 + 3D15|392|1/2|
|2D14 + A15 + A16|392.25|3/4|
|2D14 + A15 + D16|391.75|1/4|
|2D14 + D15 + A16|391.75|1/4|

Every other actual formal forest is at most391.5. Equality is allowed
in that statement; it does not count as a positive obstruction.
These are actual type costs, not the representative bounds of the
displayed D forests. Formal forests here are potential finite vertex
diagrams, not asserted geometric realizations in a polytope.

There is a sharper form of the existing necessary vertex inequality.
Let n_F count finite vertices of each of the eight types, and let e_F
be the excess in the last column. All other finite vertices have
nonpositive excess. The unchanged local face theorem, ideal surplus
at least91/2 and the single4350 correction imply

    sum_F e_F*n_F >= 4350 + (91/2)*v_i.

In particular a hypothetical dimension59 polytope needs at least2175
vertices **among these eight types**, and at least2198 if it has an
ideal vertex. This refines the previous bound on the total number of
finite vertices. It does not bound that number from above or prove
dimension58. Excluding or further reducing these remaining types is
a specific next target.

## Verification

Generation chunk56e15f and separate fresh replay chunk2c4774 both
completed exit0/PASS with1164 actual forests and eight obstructions.
The checker verifies the frozen
source bindings and report hash, the complete expanded-vector
bijection, all residual rank4 partitions, actual-type products and
all1164 exact costs. No geometric census, Bellman optimization or
ideal convolution is repeated.

The ideal-lane and arithmetic-lane peers each completed a full
source/report mathematical read. The arithmetic reviewer also read
the complete frozen combined producer/report and the inherited finite
type/representative table. Both verified the complete expanded-domain
reduction, all residual partitions and actual-type choices, the I2*
equivalence, component-order deduplication, eight obstruction costs,
and direction of the weighted-excess inequality. Both returned PASS
with no mathematical correction. Neither executed the producer,
1164-cost enumeration, Bellman proof, geometry or ideal convolution.
Metadata checkdc3bdd passed four source bindings, three local links
and Python syntax. This packet is final and frozen; no external
refereeing or proof-assistant formalization is asserted.

    python3 scripts/fresh_ambient59_exact_typed_obstruction_frontier_20261018.py --verify
