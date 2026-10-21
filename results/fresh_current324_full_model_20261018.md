# A changed global arithmetic model from the current324 local charges

This packet concerns the finite exported compact arithmetic H13 model over
Q(sqrt5). It integrates newly proved local charge columns and changes the
optimized face masses. One changed-model solve gives the exact enclosure

    12.35590648805293 <= optimum <= 12.355919609844.

Both endpoints remain above 12. The preceding optimized enclosure began at
12.35984031684756, so the objective improves modestly while the optimized
face masses change substantially. This is not a dimension theorem or a
geometric realization.

The preceding global model contained439619 columns, including209 additions
from the earlier incidence and geometry packages. Since that solve, three
complete local-cover packets produced47 further exported columns:16 from
the290-class reconstruction at its new5121 masses,27 from the32-class
historical extension, and4 from the two actual raw4 classes. Every one of
these47 columns strictly violates the preceding optimized dual and is
absent from its entire matrix, not only absent from the209 later additions.
This is therefore a genuinely changed optimization problem.

The frozen all89 raw4 barrier is a different fixed-mass diagnostic. It is
not inserted as a packing column, and its retained unknown is not reclassified.
No raw4 pair census is performed here.

## Saved charge audit and exact matrix append

Instead of retaining only the47 violated exports, the producer checks every
saved charge certificate in the three input packets. There are738 records:
571,161 and6 respectively. It verifies3273 complete-domain charge inequalities
and their attaining cover certificates at the old 5121 face masses. Fifty-nine
records have zero charge. The679 nonzero records give647 distinct primitive
columns after rational scaling and deduplication.

For a saved nonnegative orbit charge w_i, the exported coefficient of face f
is sum_i w_i c_if, with c_if counting actual corner occurrences. Its vertex
cost is1. Charge feasibility against the complete count-vector domain makes
this inequality valid independently of the old face masses. The old attaining
mixtures certify the old local optima only; they are never assigned to the
newly optimized face masses.

The producer reconstructs the full frozen model and compares each new charge
with every positive inherited column after joint primitive normalization of
both face and vertex-cost coefficients. Of647 distinct charges,444 exactly
match inherited columns and203 are appended. Every679 original nonzero
certificate origin is retained, including all47 violated exports.

| Quantity | Count |
|---|---:|
| Inherited columns |439619|
| Appended columns |203|
| Total columns |439822|
| Active columns |173996|
| Face requirements |82806|
| Vertex budgets |6558|
| Inherited zero-budget signed columns |96176|
| Inherited positive-cost signed columns |52|
| Preserved two-owner columns |14|
| Face-matrix nonzeros |829446|

The exact matrix digest is
`0ddd1b3a76ca6041463b6f1e149efa44710eb12d0261472e5f55be0ddde7ee98`.
Every old variable record, face requirement and coefficient, budget
coefficient, zero-budget record and provenance mapping is preserved
literally. All14 multi-owner columns retain both actual vertex costs.
Every new charge is also checked against the current fixed-mass repair
12.119932746175158 before optimization.

## One solve and exact recovery

The producer guards the solver invocation count. It permits one300-second
HiGHS interior-point solve of this changed matrix. Primal recovery and exact
verification use the inherited full budget matrix and strictly positive old
repair columns. Dual recovery uses every actual vertex-budget cost, applies
the inherited zero-budget simplex correction, and makes only nonnegative
vertex-multiplier repairs. It then verifies every column exactly.

Preparation 8551/chunk642393 and the sole solve 89378/chunk5ee521 both
exited0/PASS. HiGHS terminated optimally after89 iterations in70.022 seconds,
at numerical objective 12.355919609276977. Exact recovery gives

    lower = 222628992286411620/18018021785910581,
    upper = 3088979902461/250000000000.

The displayed lower decimal is rounded outward below the exact fraction.
Their gap is
59107179334521447839841/4504505446477645250000000000, approximately0.00001312.
The primal certificate has23662 positive rounded variable numerators and
repair delta3/500000000000;30 appended columns have positive primal weight.
The dual has5073 positive face masses,1095 vertex allowances and1312 exact
allowance-repair records. No second optimizer call is made.

Separate full-matrix replay 35853/chunkcd25ad exited0/PASS. It rebuilt the
same matrix and checked all exact primal face/budget inequalities and dual
column inequalities without an optimizer.

## Changed face masses

The old 5121-face support becomes 5073: 1002 previously zero face masses become
positive and 1050 previously positive masses become zero. A direct saved
Fraction comparison (chunka45716) gives normalized face-mass L1 change
approximately0.1796247365. Here each certificate is normalized by its own
sum of vertex multipliers before comparison.

| Raw4 face | Former mass | Newly optimized mass |
|---|---:|---:|
|14792|15906303/5000000000|548451/10000000000|
|14927|31081609/10000000000|14837/5000000000|

The new optimized certificate assigns much less mass to these two
strongly constrained faces. Of the 324 geometrically represented classes, 240 have positive new
mass, compared with 261 previously. Of the 1002 newly positive faces, 982 lie
outside those 324 classes. These support counts are saved-data comparisons,
not new geometric provenance or coverage claims.

The frozen all89 raw4 barrier concerned the previous5121 face masses and
cannot be transferred to these new masses. Likewise the complete local
covering optima of the three input packages have not been recomputed at the
new 5073 masses. Their fixed charge inequalities remain valid and are the
only information from them used in this global matrix.

## Independent saved checks and limits

The standard-library saved check passed with738 local certificates,
3273 charge inequalities,2016 mixture terms,679 source records,647 distinct
charges,203 appended columns and444 matches, plus all14 multi-owner loads.
Final helper execution 40377/chunks3a453c and431c88 exited0/PASS; the earlier
identical arithmetic check a282c2 also passed before adding the support-count
display. The final display independently confirms5073 positive faces,1002
added/1050 removed, and30 appended columns used by the primal.

The standard-library helper imports no producer or recovery engine. It
reconstructs all738 saved local charge/cover checks, every primitive new or
reused charge identity and all679 origin records. It checks the final new
dual-column loads, every inherited multi-owner load, saved primal
nonnegativity and new-column budget contributions, and exact interval
arithmetic. Matching reused charges to the old full matrix, old-column
inequalities, all primal face obligations and the face-weighted numerator
remain the scope of the separate full-matrix replay. No historical native
pair, graph or clique census is repeated.

The parent completed a full producer, helper and draft-report mathematical
accounting review, also reading the inherited primitive-scaling, cost-row,
primal/dual recovery and full-verification implementations. It checked all
679 origins,203 appends/444 matches,47 old-dual violations, preservation of
the old439619 columns and all multi-owner costs, exact recovery and the
old-mass-only cover scope. That review did not repeat matrix reconstruction,
arithmetic, optimization, native geometry or clique enumeration. The parent subsequently completed the full final report and helper-delta
read (terminal read bb5831), with the same source/math scope and no
computation replay.

These647 selected charges do not exhaust the complete implicit current324
local packing cones. After face masses change, the saved local optima and
covers must be recomputed before claiming complete cone saturation at those
masses. Missing face classes, finite-forest and field coverage, and global
geometric consistency also remain outside this finite matrix. Neither a
remaining above 12 dual nor an eventual target12 primal alone settles a new
dimension theorem.

Files: [producer](../scripts/fresh_current324_full_model_20261018.py),
[exact model and certificates](fresh_current324_full_model_20261018.json),
[independent saved checker](../scripts/check_current324_model_saved_20261018.py).
Inputs: [previous global model](fresh_cross_vertex_current_full_model_20261018.md),
[current290 covers](fresh_cross_vertex_current290_local_covers_20261018.md),
[historical32 extension](fresh_cross_vertex_historical32_local_covers_20261018.md),
[current324 covers](fresh_new5121_raw4_two_face_local_covers_20261018.md).
