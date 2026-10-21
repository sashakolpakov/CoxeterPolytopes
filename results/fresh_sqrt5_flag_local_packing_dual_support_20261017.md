# The new obstruction survives removing unexported quadrilaterals

The exact dual of the
[flag/local-packing successor](fresh_sqrt5_flag_local_packing_full_model_20261017.md)
has ratio12.7939560363. Removing all1,820 quadrilateral rows without
exported actual flag provenance still leaves an exact dual of ratio
12.5077183952. Thus the current obstruction does not depend on those
unexported quadrilateral rows.

In particular, merely adding flag-positive coefficients on those omitted
rows cannot defeat this restricted witness: their masses are zero, while
the corresponding negative simplex-capacity terms are nonpositive. This
does not rule out stronger incidence relations or new payments on the
retained classes.

## Literal support and feature domains

The full new dual has4,428 positive face masses and2,104 positive forest
masses;2,093 face indices were absent from the preceding optimized dual.
Its positive support splits as follows:

| Class | Positive rows |
|---|---:|
| Selected actual quadrilateral provenance, with flag features |355|
| Unselected merged quadrilateral rows, without flag features |1,820|
| Simplex rows with legacy mixed/integer and full-block features |1,467|
| Simplex rows with full-block features only |153|
| Simplex rows with legacy mixed/integer features only |397|
| Simplex rows outside both feature domains |236|

The selected flag domain contains1,858 actual prism provenances. The
unselected quadrilateral rows remain inherited merged rows. Extending
their flag coefficients would require a consistent expansion into actual
provenances; assigning independent coordinatewise minima does not provide
that expansion.

The full-block domain contains3,574 simplex rows with reconstructed
literal block features. The legacy mixed/integer domain contains5,462
simplex rows; their intersection has3,135. Having full-block features
available is distinct from adding every possible local-packing charge:
the successor contains only its specified finite175 positive cuts.
The old740-forest complete local covers were certified at the preceding
face masses. They are not certificates for this new dual or these
restricted masses. Testing the full local packing cones on the new
restricted support remains a separate task; it was not duplicated here.

## Exact restrictions, including forced simplex closure

The [checker](../scripts/fresh_sqrt5_flag_local_packing_dual_support_20261017.py)
builds the entire frozen successor once per run. Each test retains the
original exact masses on a named subset, then applies the successor's
exact dual repair. All82,806 face obligations,6,558 vertex budgets and
94,481 signed columns remain, including all77,391 flag transfers.
Quadrilateral masses never change during repair. Required containing
simplex masses may be restored or increased.

| Retained seed | Seed rows | Restored simplex rows outside seed domain | Exact repaired ratio, rounded |
|---|---:|---:|---:|
| All simplex rows and selected actual quadrilaterals |2,608|0|12.5077183952|
| Legacy exported mixed/integer classes |2,219|0|12.3842565597|
| Full-block simplex classes and selected actual quadrilaterals |1,975|113|12.4610172878|
| Simplex classes with all feature families, plus selected actual quadrilaterals |1,822|113|12.3743445180|
| Complement of the full-block/selected-quadrilateral domain |2,453|482|2.2309796126|

The first two rows are genuine restrictions to their stated domains.
Their exact ratios are respectively

\[
 {2704367297060214240\over216215876598467731}>12,
 \qquad
 {2676324107975525280\over216106965733328291}>12.
\]

The next two rows must not be called exported-only certificates: their
signed constraints force113 simplex rows outside the seed domains.
Those rows are explicitly listed in the
[JSON](fresh_sqrt5_flag_local_packing_dual_support_20261017.json). They
identify a remaining feature-coverage issue; this diagnostic does not
assert that completing those features would remove the obstruction.
The last row's low ratio is only one restricted witness, not a feasibility
certificate for that domain.

All original and restricted duals are verified against the same exact
model. The report does not infer realizability from a dual, a dimension
bound from a partial domain, or primal feasibility from a restriction
falling below12. There is no optimizer call or arbitrary provenance split.
Generation and fresh saved replay both passed. Syntax, direct source
hashes and local links passed. All predecessors remain unchanged.
**FINAL/FROZEN.**
