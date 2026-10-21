# Coupled E8 caps and direct rank-five Gale configurations

**2026-10-14. No new finite-volume hyperbolic Coxeter polytope or dimension is established.**

The H20 construction now has 87 actual facets. Eight new walls close the
previously saved escape and create eight finite vertices joined by four
finite edges. A different exact spacelike escape proves that this chamber
still has infinite volume. The compact H9 search moves outside the completed
truncation families, but all 1,063 saved direct Gale incidence patterns are
excluded by exact necessary Coxeter constraints.

## Literature restrictions applied before search

Noncompact arithmetic reflection groups can occur only in dimensions at most
19 or exactly 21. Esselmann's result excludes dimension 20 and every dimension
above 21 in that category; Borcherds realizes dimension 21. This is not a
ceiling of 21 for compact arithmetic groups: Vinberg's general arithmetic and
compact upper bound is 29. The current targets are a nonarithmetic finite-volume
H20 example and a compact H9 example. H20 would fill a missing dimension;
compact H9 would raise the known compact dimension record.

Primary references checked again this turn:

- [Esselmann, original paper and publisher abstract](https://www.sciencedirect.com/science/article/pii/S0022314X96901419), maximal integral Lorentzian lattice rank 22.
- [Borcherds, Section 8, Example 5 and the 1998 remark](https://math.berkeley.edu/~reb/papers/al/al.pdf), the rank-22 lattice and Esselmann's uniqueness in rank at least 21. Lattice rank is one more than hyperbolic dimension.
- [Vinberg, original nonexistence theorem](https://www.mathnet.ru/rus/mmo444), ruling out compact or arithmetic reflection groups in hyperbolic dimension at least 30.
- [Allcock, *An alternative to Vinberg's algorithm*, introduction](https://web.ma.utexas.edu/users/allcock/research/edgewalking3.pdf), expressly restricted to lattices over the integers. Walking edges does not make the published algorithm an unrestricted nonarithmetic generator.

The [existing literature note](arithmetic_dimension_bounds_and_effective_high_dimension_search.md)
records Vinberg's original source, the precise noncompact dimension range,
and the recent checks on known examples. No rational-lattice search is run
in the excluded dimensions. Our nonarithmetic edge extensions are separate
geometric constructions, with explicit pair, facet, corner and volume checks.

## H20: coupled interior caps and four finite edges

The [construction report](fresh_h20_e8_interior_caps_20261014.md),
[exact construction certificate](fresh_h20_e8_interior_caps_20261014_exact.json)
and its [producer](../scripts/fresh_h20_e8_interior_caps_20261014_exact.py)
retain all 79 predecessor walls. Two norm-two caps have two nonzero pi/3
contacts with the E8 support, at positions `{3,6}` and `{4,5}`. Their residual
vectors lie inside the plane of an edge of the completed H12 face. This
differs from caps directed along a single source vertex.

The original pair is orthogonal. Applying all four independent sign changes
in the two added coordinate directions gives eight distinct compatible caps.
Exact arithmetic in Q(sqrt(2),sqrt(79),sqrt(11)) verifies all 660 new wall pairs:
144 right angles, 20 pi/3 angles, and 496 disjoint pairs. All 87 inequalities
are actual facets with a common strict timelike interior point.

The complete finite-vertex census on the eight new facets gives eight vertices:
four of type `A11+A7+A2` and four of type `A11+E7+A2`. Four actual finite edges
join the two types. The 35 previously certified vertices of the H12 face
and both completed ambient cusp links are preserved. This is not a census
of every vertex or cusp of the full chamber.

The original pair cuts the saved 79-wall escape. The full eight-cap orbit
also cuts the next numerical escape, but the global witness search finds
another point. Its exact squared norm is

    1336028658157091827 / 250000000000000000 > 0.

It pairs strictly negatively with all 87 wall normals and lies in the future
affine chart. Hence the chamber has infinite volume. Local corners and bounded
individual cusp links do not replace a global volume proof. The retained
non-quasi-arithmetic obstruction remains conditional on achieving finite
covolume; this chamber is not claimed to be a lattice.

The [independent 87-wall review](fresh_h20_e8_interior_caps_20261014_exact_review.md)
passes using a separate nested-quadratic arithmetic implementation with
relative-norm sign tests and no project or producer imports. It reconstructs
the eight normals, full Gram matrix, all 1,680 candidate corner subsets,
eight actual vertices, four edges, facet witnesses, preservation signs and
the exact remaining escape.

The new E7-supported H13 face has exactly thirteen facets, of type
`A11+A2`; the other 67 nonincident wall projections are future timelike.
It is an infinite-volume cone with one finite apex. Exact normalization
gives the rational intrinsic form `diag(C_A11,C_A2,-2)`, determinant -72.
The [face report](fresh_h20_e7_face_20261014.md) and
[certificate](fresh_h20_e7_face_20261014.json) also restrict every
new actual face facet: its ambient wall must be orthogonal to E7, or have
exactly one pi/3 contact at wall 33 and be orthogonal to the other six E7
walls. In the second case its projected squared norm is one quarter of
its ambient squared norm. This extra option was absent for E8.
The [independent face review](fresh_h20_e7_face_20261014_review.md) passes,
including every projection and the rational normalization.

The known reflective form `-2t²+x1²+...+x13²` is a possible source suggested
by [McLeod's thesis, Section 3.1.4 and Table E.2](https://etheses.durham.ac.uk/id/eprint/7743/1/thesis31072013.pdf).
However, the [exact source-corner probe](fresh_h20_e7_face_20261014_source_corner_probe.json)
checks all 2,380 relevant subsets of its published nineteen-wall chamber
and finds no `A11+A2` corner. Thus that specific direct corner transfer
does not work. A [separate rational-form calculation](fresh_h20_e7_face_20261014_rational_form.json)
does establish rational similarity: the face form is rationally isometric
to three times `diag(1^13,-2)`. It has Hasse invariant -1 at primes 2 and 3
and +1 at every other place; the unscaled diagonal form instead has +1
everywhere. Exact Hilbert-symbol calculations and independent PARI/GP
products agree. This supplies an existence certificate for a rational map
after scaling, not an explicit map or an integral reflection-host
identification. Other lattices, reflection subgroups, and face completions
remain separate possibilities; no
claim of having identified a suitable complete source is made.

## Exact restrictions on further caps

The [all-residual-vector gate](fresh_h20_e8_completion_gate_20261014.md)
uses eleven exact Lorentz-Cauchy dual certificates, proposed by a numerical
optimization and subsequently checked without an optimizer. A
[separate arithmetic replay](fresh_h20_e8_completion_gate_20261014_review.json)
reconstructs every certificate using independent field multiplication and
rational interval signs. It also checks two null-ray bounds.

With the 79-wall core retained, no escape-cutting cap can have a single E8
contact at positions 0, 1 or 7, for any Coxeter label or disjoint-wall distance.
A single contact at position 2 must be pi/3. Six two-position supports are
excluded for arbitrary independent magnitudes at least one; `{0,7}` is also
excluded for equal magnitudes. These are restrictions on the listed supports,
not a noncompletion theorem. The successful new patterns are outside them.

The [transverse-edge gate](fresh_h20_e8_transverse_edges_20261014.md)
separately covers all 280 families formed from 35 H12 vertices and eight E8
fundamental weights, with a cap orthogonal to the 19 retained edge supports.
It proves that every compatible cap in these families retains the old escape,
for all real parameters and all Coxeter labels. Independent rational interval
replay passes. General residual directions and multiple E8 contacts remain
outside this statement.
The [separate mathematical proof review](fresh_h20_e8_transverse_edges_20261014_proof_review.md)
also checks orientation, monotonicity, the full Coxeter coefficient gap and
the finite-extension escape argument.

## Compact H9: a different realizable incidence family

The [direct rank-five Gale report](fresh_h9_rank5_gale_20261014.md) and
[ledger](fresh_h9_rank5_gale_20261014.json)
records a bounded search on fifteen facets: 3,000 initial proposals and
four walks of 10,000 proposals each. The 1,063 saved distinct labeled nerves
have minimal nonfaces of size at most five and no simplex-boundary join
factor in any vertex link. The latter is a sufficient obstruction to being
a standard single-face truncation, independently of the previously saved
parent catalogues.

Exact necessary Coxeter conditions exclude every saved nerve: 902 literal
local-pattern applications and 161 two-Lanner-block applications. No case
reaches a complete label survivor or a metric solve. This is a bounded
search result, not a classification of all fifteen-facet H9 polytopes or
of all Gale configurations. Separate integer incidence verification checks
3,192,189 five-row minors, 5,320,315 circuit signs, and 15,945 vertex links.
The individual audit describes which local obstruction proofs it inherits.
The [final certificate](fresh_h9_rank5_gale_20261014_certificate.json) also
gives an explicit rational primal realization of one saved type with all
310 supporting affine functionals. It freshly replays successfully. This
is a Euclidean realization of the incidence type, not a Coxeter metric.
An additional [proof and primal-realization review](fresh_h9_rank5_gale_20261014_review.md)
passes with a separate Fraction replay of the fifteen points and all 310
supporting functionals; it does not duplicate the full 1,063-case census.

## Integrity and next search boundary

The [archive auditor](../scripts/audit_e8_caps_and_rank5_checkpoint_20261014.py)
checks recursive source hashes, saved case references and local links; its
[result](e8_coupled_caps_and_rank5_gale_2026_10_14_audit.json) binds this checkpoint.
It also verifies that the previous checkpoint's frozen mathematical sources,
artifacts and reports remain unchanged. README and research-plan hashes in
the previous archive are historical entrypoint snapshots; this new archive
binds their updated contents.

Continue from the actual new corners and their supported faces, or use
contact profiles outside the proved gates. The new E7 face supplies a
rational 13-dimensional cone and a proved restriction on the walls that
can complete it. Do not repeat the closed
single-vertex transverse caps or treat the bounded Gale batch as all possible
compact H9 incidence types.
