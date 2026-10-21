# Prism-base capacity at a vertex–coface flag

For a maximal compact simplex face S, a vertex v of S, and an immediate
coface U of S, at most one counted maximal prism base in U contains v.
This is the pointwise statement already proved in the
[shared-base theorem](fresh_prism_simplex_shared_capacity_20261017.md).
It yields a refinement of the
[vertex-weighted capacity](fresh_vertex_weighted_prism_base_capacity_20261017.md):
retain U rather than summing its n−dim(S) possible values.

The exact flag constraints detect4,298 violations of each saved witness.
Monotone repair lowers the exported-only witnessed ratio from12.6142151747
to10.6043152578, and the full optimized-alpha ratio from12.9082390023 to
10.7310811081. These are repaired witness values, not optimized bounds or
primal feasibility results. No new dimension bound follows.

## Capacity and literal identification

Let I(S) be the actual incidences (Q,F) of the selected maximal prisms
along simplex bases F contained in S. A square contributes all four edge
bases, as in the inherited theorem. The selected objects are the proved
dotted full-wall-star maximal-prism aggregates, not arbitrary overlapping
products. Each incidence has a unique immediate
coface U(Q,F) containing S and Q. At any fixed pair (S,U), distinct counted
bases have disjoint vertex sets. Thus, for every flag (v,U),

\[
 \sum_{(Q,F)\in I(S)}
  [v\in F\text{ and }U(Q,F)=U]\le1.                 \tag{1}
\]

For any invariant class O of such flags, sum(1) over O. If c(S,O)=|O|
and m(rho,S,O) counts these flag occurrences in one actual prism
provenance rho, global frequencies obey

\[
 \sum_\rho m(\rho,S,O)f_\rho\le c(S,O)f_S.         \tag{2}
\]

There is no additional factor n−dim(S) in(2): the coface is already part
of the flag. Summing over all cofaces recovers the vertex capacity. This
also applies to compact counted faces in a finite-volume ambient chamber,
by the [finite-volume transfer](fresh_finite_volume_prism_capacity_20261017.md).
It makes no assertion at a nonsimple ideal vertex.

The saved prism template uses names that need care. Its `deleted_normal`
b is the facet **opposite** the selected base F; removing b exposes the
star of F and its unique maximal containing simplex S. Its
`opposite_base_facet` q is the other root of that two-element prism
factor, and **contains** F. It remains a supporting white root of S.
The vertex v is identified by the omitted black root p of S. Removing q
from S's supporting walls gives U(Q,F). In the full literal diagram, the
checker verifies

`supporting_whites(S) intersect prism_black_roots = {q}`.

Hence the ordered flag label is (p,q), with p black and q supporting
white. The order and the two different roles cannot be discarded.
Every actual prism vertex is taken from the inherited complete prism
partitions, not from arbitrary positive-definite two-deletions that may
extend a base to a larger simplex.

## Complete star orbits, including identical finite components

The complete star of S is the disjoint union of one connected Lorentz
diagram and a finite spectator forest. Its black roots all belong to the
Lorentz component. Full diagram automorphisms preserve this unique
component.

If q lies there too, the checker computes simultaneous ordered (p,q)
orbits. It uses exact edge labels, verifies every returned full graph map,
and preserves the set of black spherical deletions. Separate orbits of p
and q would not be sufficient.

If q lies in an orthogonal finite component, an orbit is the product of
the pointed Lorentz-core orbit of p and a rooted finite-component orbit
of q, allowing every copy of that same finite component in the full
spectator forest. Its size is

\[
 |O_p|\;|O_q|\;\#\{\text{copies of the finite component}\}.          \tag{3}
\]

Deleting b may have separated q into a newly finite component. The
saved Lorentz-core map then does not contain q. The checker reconstructs
that literal component, verifies a rooted isomorphism to the named
finite model, and includes identical components inherited from the
original spectator forest in(3). Counting only the newly separated
component would give an invalid capacity.

The exact sums of these flag capacities recover every prior pointed
capacity. The same summation of actual positive incidence coefficients
recovers every pointed incidence column. These checks include zero-mass
provenances and flags with zero selected positive coefficient.

## Exact diagnostic scope

The [source](../scripts/fresh_sqrt5_flag_refined_base_capacity_20261017.py)
uses the frozen [pointed-star domain](fresh_sqrt5_pointed_simplex_base_capacity_20261017.md),
the literal [base-vertex maps](fresh_sqrt5_vertex_refined_base_capacity_20261017.md),
and the [exported-only dual](fresh_sqrt5_labelled_packing_dual_support_20261017.md).
Its complete common provenance distribution is retained; there are no
independent coordinatewise minima or changes to prism masses.

For each new inequality, it first measures the exact initial signed load.
It then retains or increases each containing-simplex mass to the maximum
required load/capacity, and invokes the inherited exact dual checker to
restore and verify every existing signed and positive-cost constraint.
Every new flag inequality is checked again against the resulting masses.
No optimizer is invoked. Unselected positive incidences remain zero,
which supplies valid weaker columns rather than asserting that those
geometric incidences do not occur.

A repaired ratio is a witnessed lower bound for this partial payment
model, not its optimum. A failed old witness alone would not establish
primal feasibility, geometric realizability, or an ambient dimension
bound. All predecessors remain unchanged.

## Exact results and verification

The diagnostic has287 connected simplex-core types,4,969 complete simplex
profiles and77,391 flag columns, refining11,412 pointed columns. It
reconstructs all6,016 bases of the1,858 actual prism provenances and all
14,284 base-vertex flags. Of these,7,993 have q in the Lorentz core and
6,291 have q in a finite component. There are18,222 exact simultaneous
core-isomorphism tests and252 saved rooted finite-component maps.

Both initial witnesses require1,767 containing-simplex mass increases;
no prism mass or positive support changes. The exact repaired ratios are

\[
 \text{exported only: }\quad
 {9280034639332505\over875118705335371}
   =10.6043152577525\ldots,
\]
\[
 \text{full optimized alpha: }\quad
 {18979113911719330\over1768611542538437}
   =10.7310811081099\ldots.
\]

All inherited model costs and signed transfers are rechecked. The
separate740-forest local packing covers were constructed at the old face
masses; they are not silently reused after these increases. Any
combination with those covers must supply new cover marginals or an
explicit valid incremental repair.

The [JSON](fresh_sqrt5_flag_refined_base_capacity_20261017.json) records
every map, occurrence, capacity, initial signed load and repaired exact
dual. Generation and fresh saved replay both passed. The
construction_audit agent's independent narrow mathematical/source read
passed the flag theorem, literal q identification, full spectator
multiplicities and signed-repair logic. That read does not independently
replay the arithmetic. Syntax, direct source hashes and local links
passed. No optimizer was invoked. **FINAL/FROZEN.**
