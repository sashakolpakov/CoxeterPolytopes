# A forest-sensitive partial rule reaches the rank29 vertex budget

The [checker](../scripts/fresh_sqrt5_forest_sensitive_weights_20261015.py) and
[certificate](fresh_sqrt5_forest_sensitive_weights_20261015.json) give a
nonnegative weight rule with **total weight at most84 at every rank29 finite
vertex**. It pays all112 previously sealed quadrilateral classes, including
arbitrary external white diagrams, and preserves the established triangle,
short-path and H4 payments. Complete quadrilateral coverage remains open;
this is not yet a new dimension bound.

Start with the saved fixed-baseline pair weights from the
[sealed-core calculation](fresh_sqrt5_sealed_weights_20261015.md). Change only
the bonuses in an entire D7 component. If the full vertex diagram has an
entire A1 component, retain the old D7 bonuses. If it has no A1 component,
retain bonus1 on the two marked-pair orbits forced by classes535/563,
bonus1/2 on the two orbits appearing twice in classes785/860, and zero on
the remaining D7 orbits. All other component bonuses stay unchanged.

The D7 totals are therefore16 without A1 and21 with A1. The latter condition
is essential: four sparse D7 forcing faces have an isolated white spectator,
and the reviewed unpaid-sealing argument forces it to remain an entire A1.
The two additional half-weight orbits are also essential; using only the
two sparse pure-D7 witnesses would leave classes785/860 unpaid.

For each of the112 sealed cores, the existing complete one-white contact
argument says that an external white node touching the core either fails a
necessary test or creates an existing short-path/H4 payment. In the unpaid
case every external white node is orthogonal to the entire core. Every
corner diagram is consequently its displayed core diagram plus the same
external white diagram. Testing whether that external diagram has an A1
component gives exactly two cases per core. All224 exact face sums are at
least1. Additional pair bonuses are nonnegative, so all existing baseline
payments remain valid, including the reviewed arbitrary-extension triangle
rule.

The checker enumerates all92 allowed irreducible finite component types up
to rank29 and uses exact dynamic programming over their ranks. The maximum
without an A1 is159/2, attained by D6+E7+E8+E8. With an A1 it is84, attained
by A1+D7+D7+D7+D7. This dependence on the entire vertex diagram bypasses
the old component-only lower bound86.

The remaining obligation is substantial and explicit: the130 nonsingular
finite unpaid classes with458 accepted one-white profiles, and24 symbolic
dotted classes with369 profiles, still need coverage under all external
white extensions. The44 unpaid-sealed nonsingular cores and two certified
affine dotted cores remain excluded by the inherited argument. No extension
cutoff or sampled rank calculation is being promoted to a proof.

[Independent review](fresh_sqrt5_forest_sensitive_weights_20261015_review.md)
passes: it reconstructs all face cases and vertex totals separately. The
historical component-only
obstruction remains valid within its stated weighting class; it does not
apply to weights depending on whether the full vertex diagram contains A1.
