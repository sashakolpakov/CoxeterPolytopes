# Direct Gale search and two-coordinate dimension raising

**14 October 2026. No new finite-volume or compact Coxeter polytope has
been established.** Both specified compact-H9 surgery catalogues are now
closed. Direct Gale configurations supply further realizable types, all
excluded in the bounded search. Two different H20 constructions pass
their exact Coxeter checks but fail global volume tests. A new structural
obstruction stops additions that cannot change the problematic facets.
This supersedes the frontier in the
[preceding checkpoint](changed_endpoints_and_complete_surgery_exclusion_2026_10_14.md).

## Literature determines the targets

The dimension-21 restriction applies to **noncompact arithmetic**
reflection lattices: the exact range is `n<=19` or `n=21`. Dimension 20
and every dimension above 21 are excluded in that category. Compact
arithmetic examples have the different bound `n<=29`; Vinberg also
proved that bound for all compact Coxeter polytopes. See
[Esselmann](https://doi.org/10.1006/jnth.1996.0141),
[Borcherds, corrected thesis §6.8](https://math.berkeley.edu/~reb/papers/thesis/thesis.pdf),
[Vinberg](https://www.mathnet.ru/eng/mmo444), and the
[detailed dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md).

The sources checked still give known compact examples only through H8,
and finite-volume examples through H19 and in H21. H20 remains absent
from the known dimensions. The
[July 2026 construction paper](https://arxiv.org/html/2607.14715v1)
does not raise these records, and the
[October 2026 revision of the d+4-facet classification](https://arxiv.org/html/2608.02894v2)
does not settle unrestricted existence. Compact H9 and nonarithmetic
finite-volume H20 remain the active targets. A noncompact arithmetic
lattice enumeration in H20 or H22+ is ruled out before computation.

[Allcock's published edgewalking](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf)
uses an integral Lorentzian lattice with a rational bilinear form.
It finds the next lattice mirrors from corners. Vinberg orders mirrors
by distance from a control point. Edge traversal is useful more broadly,
but an unrestricted combinatorial construction requires its own metric
realization and global completeness proof; it does not inherit those
guarantees from the published lattice algorithm.

## Compact H9: close the finite catalogues, then change the incidence

The [changed-order surgery report](fresh_h9_newfacet_frontier_20261014.md)
closes all **5,536** admissible types obtained by truncating the fixed
Bugaenko H8, wedging the new facet, and truncating again, with the two
cuts of codimension two through five. This includes all 3,244 cases
left open previously. New induced-nerve obstructions close 3,112 of
those cases; complete propagated label searches close the remaining
132 in 1,465 evaluated states. No necessary full labeling survives,
so no numerical metric fit contributes to the exclusion.

Separate checkers reconstruct all 5,536 incidences, verify the 44 local
obstructions and the 96-row overlap relation, and replay the remaining
search arithmetic. The report distinguishes shared classical tables
and shared DFS control flow from independent incidence and arithmetic.
The preceding 2,953-type wedge-first catalogue was already closed.
Neither result classifies all fourteen-facet compact H9 polytopes.

The [direct Gale search](fresh_h9_gale_frontier_20261014.md) changes the
source of incidence data. Integral fourteen-by-four Gale matrices give
actual simplicial nine-polytopes; their duals have fourteen facets.
Two bounded walks produce **770 distinct realizable nerves**. Exact
obstructions eliminate 715 by induced configurations and 32 by the
two-Lannér-block theorem. The remaining 23 close by complete label
searches in 189 evaluated states. No metric survivor is found.
Exactly 749 types lie outside both closed surgery catalogues. The two
walks make 80,000 and 100,000 proposals; this is a bounded search, not
exhaustion of the rank-four Gale space. The detailed report records
sampling limits and verification scope.

## H20: three repaired cusp links do not establish global closure

The [two-coordinate lift](fresh_dimension_construction_routes_20261014.md)
starts from Vinberg's nonarithmetic H18 hybrid. Twenty-six norm-two
roots remain unchanged; each of fourteen norm-one roots is replaced
by a signed pair in one of two added orthogonal directions. Every
axis assignment passes the pair gate, but a parallel triangle in the
source forces a strict escape for all 8,192 assignments modulo axis swap.

A uniquely optimal coloring for the specified parallel-pair objective
requires three cusp repairs. Six added roots give an exact **60-facet
H20 chamber**. At each repaired cusp, positive affine kernels and full
span certify a bounded link. All 1,770 wall pairs and all actual facets
are checked. Nevertheless, the report gives a literal strict spacelike
escape, so the volume is infinite. Retained Gram cycles also prove that
any finite-volume completion keeping those walls would have to be
non-quasi-arithmetic. The current infinite-volume group is not called a
nonarithmetic lattice.

The useful next test partitions the old H18 cone according to which
linear bound controls each new coordinate, and checks all lifted extreme
rays of those cells. This is a finite global closure criterion. It has
not been run to exhaustion for a successful modified construction.

## Weighted-box constructions: coverage and facet geometry are separate

The [pentagonal-radius theorem](fresh_h20_weighted_pentagonal_20261014.md)
extends the preceding weighted-box exclusion to a four-family union:
all seed-compatible radius-two caps, radius-`4/phi` caps, literal
half-vertex radius-one caps, and arbitrary-radius old-vertex caps.
Exact adjacent-probe relations and cube-edge counting leave at least
35,328 probes not strictly cut. Boundary cases admit an explicit outward
perturbation, proving infinite volume. A separate arithmetic and proof
review passes. Other centers are outside this theorem.

The [parity-controlled binary construction](fresh_h20_weighted_binary_20261014.md)
uses such different centers. Its radius-one orbit strictly cuts all
524,288 half-vertex probes, and is compatible with all 524,288
radius-`sqrt(2)` endpoint caps. Compressed exact pair formulas, actual
facets, and the old finite endpoints have independent arithmetic checks.
A different explicit escape survives, proving infinite volume.

More significantly, a complete local graph and 73 checked colorings
prove that **none of the binary facets contains a finite vertex**.
The graph and orbit transitivity have a separate exhaustive review.
Each binary facet also has actual `pi/5` ridges. The
[dotted-only completion theorem](fresh_h20_weighted_binary_20261014_dotted_completion_review.md)
therefore stops every acute-compatible completion whose new walls are
all dotted to even one retained binary facet. Such walls leave that
facet unchanged. A finite-volume Coxeter polytope cannot have an ideal
vertex on a `pi/5` ridge: the induced Euclidean reflection product would
have irrational trace while preserving a full translation lattice.
Consequently the ridge would need a finite vertex, contradicting the
unchanged facet geometry.

There is a simpler obstruction for the retained central facet as well.
Its only possible finite vertices in the original 39-wall core would
use one wall from every coordinate pair. Such intersections have
`x_i=+/-1`, `D=2`, `U=1`, and square `16+8phi>0`, so none is finite.
Its eight golden walls nevertheless give actual `pi/5` ridges. All the
binary, endpoint and outer-code caps are dotted to the central wall,
leaving this obstructed facet unchanged. Thus central-dotted additions
are already ruled out before any large orbit search. A useful retained-core
construction must produce finite corners on those central ridges as well
as repair any retained binary facets. The same trace argument applies to
every finite ridge label outside `2,3,4,6`.

The [outer unit-code construction](fresh_h20_outer_unit_code_20261014.md)
uses a shortened binary code and three new radii. It repairs the
preceding explicit escape, but another exact escape survives. Its added
walls are all dotted to every binary facet, so the structural theorem
already prevents finite-volume completion by this route. Further
expansion of this disconnected cap family is stopped. Any continuation
retaining binary facets must introduce walls that meet them; deletion
or replacement is a different construction.

## Verification and next work

Each linked report distinguishes exact certificates, independent
arithmetic, shared search control flow, and bounded discovery. Earlier
frozen artifacts are unchanged. The
[archive auditor](../scripts/audit_direct_gale_checkpoint_20261014.py) and
[saved audit](direct_gale_and_two_axis_constructions_2026_10_14_audit.json)
check current links and recursively declared source hashes; they do not
replace the mathematical replays.

```sh
python3 scripts/audit_direct_gale_checkpoint_20261014.py --verify
```

Do not revisit either closed surgery catalogue, the excluded weighted-box
cap dictionaries, or dotted-only additions to the binary facets. The
remaining constructive directions are Gale configurations avoiding the
recorded induced obstructions, changed or coupled height constraints in
the H18 lift, and cap replacements that create complete finite vertices
on the retained noncrystallographic ridges.
All still require an exact global volume certificate. The existence
goal remains active.
