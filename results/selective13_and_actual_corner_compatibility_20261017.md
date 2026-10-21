# Selective corner payments and actual finite-vertex compatibility

This checkpoint proves an unrestricted dimension bound **146** for
finite-sided, finite-volume real hyperbolic Coxeter polytopes. It also
improves the ambient bound for chambers having a compact facet to111.
The compact-arithmetic ceiling21 and the minimum-cusp-factor-rank9
ceiling29 are unchanged. No new finite-volume construction dimension
has been established.

The new dimension proof has a separate internal mathematical/source
review and independent Fraction computations. It is a complete internally
checked computer-assisted proof, with no identified gap. External
refereeing, proof-assistant formalization and publication priority are
not asserted.

## Unrestricted156 to146

The [selective distance13 proof](fresh_selective_distance13_dimension146_20261017.md)
starts with a counterexample to the simpler idea. Two minimal L25
diagrams sharing their farthest node give the known finite-volume H17
pyramid. Its quadrilaterals need distance14. The new calculation checks
the complete chamber cone:80 finite vertices and2 ideal vertices.
This reconstructs an existing example rather than a new one.

The improvement comes from classifying that exception. Under the
assumption that all adjacent marked white paths exceed13, every other
core overlap or joining bond contradicts positive definiteness or
negative index one. For the surviving shared core, its actual radical
relation and nine preferred-corner inverse calculations force every
additional supporting white to be orthogonal to the whole core. A
nonzero attachment would instead create a path of length at most4.

The relevant corner component therefore has only17 nodes. Paying all
pairs in such a small component, along with all pairs of distance at
most13 elsewhere, costs at most13n at a finite vertex and
sum min(13,r_i/2) per cusp-link vertex. It pays every triangle and
quadrilateral, including corners at nonsimple cusps.

Convex-order concentration at rank26 reduces all cusp-rank partitions
in dimensions147–156 to1255 profiles. Their exact curvature margins
are positive, with minimum830407/9565938. The inherited156 theorem
covers larger dimensions. The same vertex payment on compact facets
gives ambient dimension at most111. The
[separate review](fresh_selective_distance13_dimension146_20261017_review.md)
reconstructs the geometry and budgets independently, including the
proper cusp links of the known pyramid.

The local profile(26,26,26,26,26,15) has negative margin at ambient146.
Thus repeating this scalar comparison does not give145. This is a
limitation of these weights, not an example of an ambient chamber.

## A new geometric constraint in the compact arithmetic model

The [first finite-corner compatibility theorem](fresh_dominant_cap_compatibility_20261017.md)
uses complete contacts of omitted walls with the same finite vertex.
Those contacts determine both the positive projection and the common
future timelike component of each actual normal. Two complete prism
stars that abstract port packing allows can consequently force a wall
product strictly between−1 and0, impossible for normalized Coxeter
walls.

At the selected forest A1+A3+H3+H4+I2(5), six literal orientations of
two prism classes are mutually incompatible. Seven pairs share an
actual quadrilateral; the other eight force product
−sqrt((5+sqrt5)/16). The resulting whole-face inequality invalidates
the preceding fixed witness. Increasing just its allowance at this
vertex repairs it to11.986190572394824, below the target12.

This removes the old fixed witness's obstruction. It does not establish
optimized feasibility or a new arithmetic dimension theorem. The
[successor family](fresh_dominant_cap_compatibility_family_20261017.md)
extends the same actual-normal method across all57 positive raw23
classes:570 finite corners,299 forests and6010 literal orientations.
Of134445 tested pairs,134435 have exact exclusions. Whole-orbit clique
constraints give294 columns;120 allowance increases repair the saved
witness further to8.804856730635848. Unexcluded pairs are conservatively
treated as compatible.
The [independent family review](fresh_dominant_cap_compatibility_family_20261017_review.md)
checks the saved corner multiplicities, clique capacities and exact radical
signs with separate arithmetic; generation and fresh replay pass.
No new optimization was performed. Full field, face and forest coverage
remain necessary for a dimension exclusion.

## The remaining fine-payment defects have nonincident carriers

The [small-cusp localization lemma](fresh_small_cusp_defects_20261017.md)
shows that a face whose original wall star contains no affine diagram
of order at most9 receives the finer Vinberg payment. Assign weight2
to every small varying-factor ideal corner. Any face still deficient
then has a small affine carrier whose actual cusp is outside that face.
The lemma gives an exact lower demand for such defects; bounding the
number assigned to each nonincident cusp is an unresolved global step.

Its separate one-cusp155 consequence was obtained from the earlier156
weights and is now superseded by the unrestricted146 result. Its
localization statement remains useful. The square-pyramid example
distinguishes a carrier cusp from an incident vertex; its finite A3
corners are paid, so it is not offered as an unpaid-face example.

## Construction: two changed cavities tested

The [three-cap gate](fresh_three_cap_pentagon_20261017.md) excludes the
specified central triangle with boundary arcs1+2+2 after deleting
walls12,47,52. Every rotation needs a forbidden mandatory cap. Exact
checks cover18 low-label norm branches and4 complete high-label
interval branches, without a label cutoff.

The [changed-seed scout](fresh_three_cap_pentagon_changed_seed_20261017.md)
also releases64 and65, replacing the boundary by[7,8,44,63,30]. The
exact deletion cost is11 finite-vertex certificates and6 compact
two-faces, with no higher-dimensional compact face touched. Its94
coherent low-label caps leave11 open ends and supply no full patch.
These are numerical diagnostics, not accepted new chamber walls.

The [exact changed-seed gate](fresh_three_cap_pentagon_changed_seed_gate_20261017.md)
then resolves all mandatory four-contact high-label families, including
the one two-variable norm locus. Exact bounds reduce its first label
to16 or17, and the forced second label lies strictly between consecutive
Coxeter values. Together with three mandatory three-contact exclusions,
this rules out the stated one-/two-center and positive-arc three-center
central-triangle closures. Other incidences, additional deletions and
larger patches remain open. Enlarging the same label grid cannot repair
these excluded configurations.

The exact construction control remains63 actual walls,199 finite
vertices and10 target-face open ends, with certified infinite volume.
The search for a new frontier-dimensional finite-volume example remains
unfinished.

## Archive and continuation

The [incremental auditor](../scripts/audit_selective13_and_actual_corner_compatibility_20261017.py)
and [saved audit](selective13_and_actual_corner_compatibility_20261017_audit.json)
bind the new reports and certificates and rehash all frozen predecessor
files. They compare all1255 producer/reviewer margins,21 inverse ports,
all294 new clique loads and120 allowance raises. Historical censuses and large
optimization models are not rerun by this archive check.

The next nonexistence steps are to classify the next genuine long-path
exceptions and to bound nonincident small-cusp carrier multiplicities.
For compact arithmetic, optimize only after adding the new actual-corner
compatibility constraints and keep feasibility separate from repaired
fixed witnesses. The next finite model can retain all exported earlier
columns and add the cap, forced-incidence, joint-cut and294 new clique
columns, for at most430666 variables. The larger implicit packing and
joint families require renewed separation for any new witness; their
old fixed-witness certificates alone do not cover a new solution.
Construction needs a different retained-wall geometry
or boundary incidence, with a complete finite-volume certificate before
promoting any chamber.
