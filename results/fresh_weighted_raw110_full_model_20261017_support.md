# Support of the weighted/raw110 finite-model dual

The new exact finite-model dual gives **12.368567257872337 >12**. Restricting
it to all positive simplex/edge rows and the209 positive Q rows in the named
248-class actual-corner union still gives **12.089357918752444 >12**. Thus
unrepresented Q rows alone do not explain the remaining finite-model obstruction.

The [source](../scripts/fresh_weighted_raw110_full_model_20261017_support.py)
and [saved diagnostic](fresh_weighted_raw110_full_model_20261017_support.json)
rebuild the [enlarged model](fresh_weighted_raw110_full_model_20261017.md),
verify both exact certificates, recount every349 appended load and compute
two exact restrictions. No optimizer or new geometry census is invoked.

## Current support and coverage

|Statistic|Count|
|---|---:|
|Positive face masses|5053|
|Positive maximal-simplex/edge rows|2873|
|Positive Q aggregate rows|2180|
|Positive vertex allowances|1025|
|Newly positive face rows|1450|
|Previously positive rows now zero|1110|
|Changed common positive face masses|3600|
|S rows inside old complete feature domain|1994|
|S rows outside that domain|879|
|Q rows with explicit actual provenance|433|
|Q rows with one unique saved context|1321|
|Ambiguous merged Q rows left unexpanded|426|

The S rows comprise1643 maximal-simplex aggregates and1230 maximal-edge
incidences. The named actual-corner domain is the union of the old243
class identifiers with the189 raw110-package identifiers, totaling248.
The raw110 package retains184 of the old classes and adds five. Classes
that formerly had zero mass are not silently removed from this union.

There are206 currently positive Q rows in the old243 domain and174 in the
raw110 package, with209 in their union and1971 outside it. Only339 positive
Q rows lie in the old joint-port export;1841 lie outside it. The86 weighted
columns all have positive current load,12 exactly tight. Of263 raw110
columns,250 have positive load and one is exactly tight. Every appended
inequality passes exactly.

An unchanged class identifier does not validate a packing mixture at changed
mass. The preceding complete weighted graph-cone covers, local simplex
packing covers and implicit joint scans are not reused. Every earlier
exported column remains in the full finite matrix. The new checks concern
those fixed coefficients, not an unexported complete charge family.

## Exact restrictions

Both restrictions retain all2873 positive S masses, discard only the stated
Q rows and recompute vertex allowances against the entire matrix. The96176
zero-budget signed columns and52 positive-cost signed columns remain in
verification. No outside S row is restored, no retained mass changes and
every resulting vertex allowance is at most its original value.

|Retained Q domain|Kept Q|Dropped Q|Positive faces|Vertex allowances|Bound|
|---|---:|---:|---:|---:|---:|
|Named248-class actual-corner union|209|1971|3082|790|12.089357918752444|
|All explicit or unique-context Q|1754|426|4627|996|12.325780504848609|

The exact ratios are respectively

    4667687600305587 / 386098883966806,
    44417119684509572 / 3603594893405505.

These are finite-model dual obstructions. Neither is an ambient realization
or a certificate of complete implicit-family coverage. Dropping every
ambiguous merged Q row still leaves an exact bound above12.

## Concrete next constraints

The largest uncovered cohorts with explicit or unique saved provenance are:

|Raw family|Classes|Actual corners|Demand-weighted mass|Dotted square|
|---|---:|---:|---:|---|
|4|83|498|0.0295293495|(7+sqrt5)/2|
|76|18|72|0.0203302832|(7+sqrt5)/2|
|205|11|66|0.0197689524|6+2sqrt5|
|125|9|54|0.0165768582|13+5sqrt5|
|57|34|136|0.0147571629|6+2sqrt5|
|14|74|592|0.0147412536|(7+sqrt5)/2|
|112|12|72|0.0144577452|19+8sqrt5|

Raw125 is a concentrated next geometric target: nine unique-context classes
give54 actual corners. Face51173 alone has mass41941921/10000000000 at
source context70783. Its dotted square13+5sqrt5 requires a correct native
square-class embedding; no new physical compatibility is asserted here.
Raw205 offers11 classes and66 corners with all dotted contacts already in
K=Q(sqrt5). Its leading explicit row61785 has provenance1398, source89825.
Raw76 has18 classes and72 corners and uses the existing raw23 dotted square.

A broader common-field opportunity consists of312 uncovered actual classes
across60 raw families with finite labels only2,3,5 and all dotted contact
magnitudes in K. They give1580 corners and demand-weighted mass
12927357/250000000. This is saved-metadata selection, not a projection or
pair-capacity theorem. Another147 positive class identifiers belong to
previously treated raw families but fall outside the named geometry export.

The represented-Q restriction is itself above12. A meaningful first check
would rebuild complete weighted graph-cone capacities at these changed
masses, including newly positive classes and conservative unknown pairs.
The old fixed-mass mixtures cannot answer that question. New mixed actual
corner constraints for a concentrated cohort are a separate geometric
route. Repeating the unchanged global model is not justified by these data.

## Verification

Generation68781 and separate fresh saved replay84551 both terminated
exit0/PASS. Each rebuilt the whole model, checked both parent certificates,
recounted349 appended loads and verified both exact restrictions. A separate
checker importing no project producer passes all saved349 loads, support
partitions, restriction objectives and normalizers, retained-mass identities
and vertex-allowance monotonicity. It does not duplicate the full matrix
reconstruction. All16 direct input bindings, syntax and3 local links pass.

Root completed full source and mathematical report reads, finding no gap
in the explicit248-class union, provenance distinction, changed-mass load
checks, signed-column-aware restrictions or coverage limitations. This was
a source/composition review, not a second matrix or physical-geometry
implementation. Source, data and report are final/frozen, with no change
to frozen predecessors.

The scope is a support and exact coefficient-load diagnostic. It adds no
arithmetic-field coverage, global realization claim or dimension theorem.
