# Whole-region cap on the C(3) companion triples

None of the100 [coherent companion triples](fresh_c3_companion_continuation_20261017.md)
admits the specified single cap closing the five-edge B boundary in this
numerical calculation. The finite label domain is complete for these fixed
prefixes, but the numerical root calculations are not an exact exclusion
certificate and do not cover continuously deformed prefixes.

The boundary is30–52–R–A–L–30. A closing cap D must give finite corners on
all five edges. Its mandatory three-contact frame is30/52/R. The retained
30–52 bond is4 and52–R is at least4. Each assigned D contact therefore
shares a required finite corner with another nonright bond of label at
least4. A label at least6 would give a norm2 acute3x3 Gram determinant

    8 - 2(a^2+b^2+c^2) - 2abc <= 8 - 2(3+2) < 0.

Thus all three assigned cap labels lie in2..5; there is no omitted high-label
tail in this frame. All seven A2 profiles are then tested against the full
finite endpoint diagrams. The100 fixed prefixes give1,010 recipes and2,020
norm branches. Every branch is checked against all retained native contacts
and the unchanged strict common center.

There are340 native-coherent caps. Every complete numerical target cone has
35 finite and7 exterior rays with11 open ends, agreeing with independent
positive-definite incidence. The two frame corners are finite, but the cap
does not meet A or L: across these candidates D.A is at most-4.7738598 and
D.L is at most-30.1780837. These are measured dotted products, not near-miss
finite angles. No candidate has the required five finite corners; none
improves the ten-open-end control. There are no unresolved numerical cases.

This failure is driven by the earlier companion geometry, not a terminal
label cutoff. A useful successor must deform those companions, use more
than one final connector, or change the incidence. Merely increasing the
three terminal labels cannot add a permitted finite-corner recipe.

No original wall is deleted. All199 old finite vertices and old compact
faces would persist under any exact coherent realization by acute-center
projection. No new ambient finite-vertex count, exact chamber or finite-volume
claim is made. The frozen63-wall/199-finite/10-open-end control remains
unchanged and infinite-volume.

[Source](../scripts/fresh_c3_companion_region_closure_20261017.py) and
[saved data](fresh_c3_companion_region_closure_20261017.json) contain all
domains, native products and complete accepted target cones. Primary
generation and fresh saved replay both pass. A separate saved-data check
confirms the finite-corner pattern true/true/false/false/false in all340 caps;
it does not independently recompute their numerical coordinates or cones.

    python3 scripts/fresh_c3_companion_region_closure_20261017.py
    python3 scripts/fresh_c3_companion_region_closure_20261017.py --verify
