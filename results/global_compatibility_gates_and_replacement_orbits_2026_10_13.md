# Global compatibility gates and the full H22 replacement orbit

No new compact or finite-volume Coxeter polytope is established. This checkpoint
tests whether the recent local constructions can close globally, and proves
that several large candidate families cannot. It follows the
[dimension-gate and cusp-replacement checkpoint](arithmetic_dimension_gate_cusp_replacement_2026_10_13.md).

The literature restriction remains essential: Esselmann excludes noncompact
arithmetic examples above dimension 21, while Vinberg's bound for all arithmetic
and all compact Coxeter polytopes is 29. Compact arithmetic examples in dimensions
22–29 are not excluded by the noncompact bound. The present H22 construction is
a geometric, conditionally nonarithmetic route; compact H9 is a separate target.
See the [dimension audit](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
and [Vinberg–Allcock algorithm audit](vinberg_allcock_search_strategy_20261013.md).

## H22: the full replacement orbit is coherent but cannot close its escape

The [full stabilizer construction](fresh_h22_replacement_orbit_20261013.md)
extends the one-wall replacement to all 5,760 symmetry images. It gives 12,331
actual facets and 5,760 certified finite vertices, with exact compatibility
for every common real parameter at least 350. These are explicit chamber
configurations, but they have infinite volume. Their saved strict spacelike
escape survives, and a checked inherited cusp retains an unbounded link.

The stronger obstruction is continuous: on each of the 5,760 residual
three-dimensional planes, the escape functional is a strictly positive
combination of three independent retained-wall functionals. Thus every acute
nonzero normal in any such plane misses the escape, regardless of its real
coefficients, height, norm or angle labels. Increasing parameter cutoffs or
changing the two caps' mutual angle cannot repair this family.

Two [literal rational vector identities](fresh_h22_replacement_escape_dual_20261013.md)
give a useful condition outside those planes: every acute wall cutting the
saved escape must pair strictly negatively with both old walls 7 and 168.
The [complete star calculation](fresh_h22_replacement_edges_20261013.md)
finds all 22 outgoing edges of the replacement vertex open. On each edge,
the required two contacts imply an exact projected-norm lower bound at
least three for an added normal of square two. Its endpoint Gram therefore
cannot be positive semidefinite. No escape-cutting addition can terminate
one of these edges at either a finite or an ideal vertex, with no label bound.

This does not exclude intermediate walls that leave this escape feasible,
new vertices reached afterward, normals outside those planes, or further
old-wall replacements. It is not a nonexistence theorem for all completions
of the retained core. The [scoped source review and three fresh exact replays](fresh_h22_replacement_global_20261013_review.md)
pass; the old group and source reconstruction remain explicit dependencies.

## H9: global enclosure must survive compatibility

The [global-bank calculation](fresh_h9_global_h5_20261013.md) starts from the
previous 125-wall chamber, which contains a compact hexagonal H2 face but whose
containing H5 face still has 53 open edges. The union of twelve preceding
proposal archives contains 38,209 literal normal recipes across 121 defining
contact spans. A common strict spacelike escape is verified against all of
their exact formula branches. Thus no subset of that finite bank can complete
the chamber, even when mutual compatibility is ignored and optional old walls
are removed. This is a finite-bank result, not an exclusion of other profiles
or all completions of the seed.

Adding 119 opposite quadratic branches initially gave a useful geometric
relaxation: with all old walls and all proposals included, numerical H5
halfspace enumeration found 122 vertices, all timelike, with largest square
about −0.00786453. Some of those walls conflict, so this is not a Coxeter
polytope. Conditioning separately on each required new branch and retaining
every wall numerically compatible with it produced a strict spacelike escape
in every one of the 119 relaxed cases. This test permits old optional-wall
deletions and ignores all remaining mutual conflicts. The enlarged-bank
diagnostic therefore supplies no compatible bounded H5 construction.

That numerical diagnosis is now promoted to an exact finite-bank theorem.
Four exact ambient escapes have complete cutter lists of sizes 15, 54, 1 and
2. All 38,453 dictionary entries (125 old walls, 38,209 original cap recipes
and 119 opposite branches) are checked against each witness: 153,812 exact
sign classifications. Every added branch is assigned one witness; 1,198 exact
incompatibility inequalities show that every possible cutter of that witness
conflicts with the assigned branch. Each inequality is either positive or
strictly between consecutive Coxeter values, so the proof has no angle-label
cutoff. A selection without any added branch has the original common escape.
Every pairwise Coxeter-compatible selection from the entire dictionary
therefore has infinite volume, even after arbitrary optional old-wall deletion.

The [original-bank review](fresh_h9_global_h5_20261013_review.md) and
[enlarged-bank review](fresh_h9_global_h5_20261013_enlarged_review.md) pass.
Separate rational checks verify escape geometry, projection, cover assignments
and angle gaps; fresh complete primary replays also pass. Neither theorem
excludes other normals or another seed. The authoritative chamber remains
the previous 125-wall construction; no rejected trial was appended to it.

## A different H9 seed with a compact H5 face

The [new construction](fresh_h9_makarov_face_extension_20261013.md) starts
with the known compact arithmetic Makarov P0
prism in H5, obtained by polar truncation of the `[5,3,3,3,3]` orthoscheme,
and a perpendicular spherical H4 root system. The base form is defined over
Q(sqrt(5)) and has signature (9,1). An exhaustive old-face ray check and a
positive dual-center identity certify the compact H5 face directly. The
11-wall H9 seed has ten certified actual finite vertices, but infinite volume.

Before expanding its cap bank, the new checker excludes all forty families
centered at those ten vertices and directed along the four H4 fundamental
weights. A single strict integer escape survives even all forty label-three
caps together. The [monotonicity argument](fresh_h9_h4_transverse_monotonicity_20261013.md)
then covers every larger label and every real transverse magnitude at least
one. This conclusion does not depend on mutual cap compatibility. A bounded
200-entry screen separately finds 152 individually compatible caps; that
count does not overcome the all-label obstruction.

Releasing the vertex-center restriction gives an exact cap centered on an
old A3 face. It has seven right, two label-three and two dotted relations
with the 11 seed walls, cuts the forty-family escape, and preserves the
compact H5 face and all ten old finite vertices. Its coordinates lie in
Q(sqrt(5),sqrt(11)). This is a coherent 12-wall H9 chamber with an actual
new facet, but a new strict integer escape proves that its volume is still
infinite. Its purpose is to leave the excluded family while retaining a
compact lower-dimensional face from the outset. The
[scoped source review](fresh_h9_makarov_face_extension_20261013_review.md)
records the exact seed, all-label obstruction and changed-center construction.
The [separate sparse-radical implementation](fresh_h9_makarov_face_extension_20261013_independent.md)
also passes: all twelve norms, all 66 pair gates, the strict center, ten
rank-nine finite vertices, the new facet, both escape checks and all forty
family inequalities. This second implementation uses Q(sqrt(5),sqrt(11))
directly instead of the producer's nested-radical arithmetic.

## Consequence for the search

Search globally before expanding another local cap bank. First test whether
the entire proposed bank can enclose the required slice at all. Then condition
on required walls, remove conflicts, and repeat the enclosure test. Exact
reconstruction and a full chamber closure certificate are still required
for any surviving selection. Local polygons, new finite vertices and larger
facet counts do not substitute for that final step.

Further H22 work must leave the proved residual planes or change retained
walls. Further H9 work must change the excluded finite bank and retain global
compatibility during selection. A different seed begins with a known compact
H5 prism and a transverse H4 diagram. The standalone
[H4 Schur gate](fresh_h9_h4_transverse_gate_20261013.md) and
[all-label vertex-centered cap envelope](fresh_h9_h4_transverse_monotonicity_20261013.md)
give exact restrictions before another bank is expanded. The repository's independent global
[H9/14 combinatorial search](h9_14_facet_completion.md) also remains open;
ordinary polytopal incidence candidates are not Coxeter realizations.
