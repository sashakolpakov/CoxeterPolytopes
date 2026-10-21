# A new arithmetic solve with cross-vertex incidence costs

The enlarged compact arithmetic H13 model over Q(sqrt5) still has optimum
above12. One bounded solve gives the exact enclosure

    12.35984031684756 <= optimum <= 12.359855171582.

The preceding optimized interval began at12.3685672578. The new constraints
improve that obstruction modestly, while the optimized face masses change.
The recent fixed-mass repair11.985461837279875 is a different witness and
does not establish target12 feasibility. This remains a finite-model result,
with no new compact arithmetic dimension theorem.

The [producer](../scripts/fresh_cross_vertex_current_full_model_20261018.py),
[exact certificate](fresh_cross_vertex_current_full_model_20261018.json) and
[separate saved checker](../scripts/check_cross_vertex_current_model_saved_20261018.py)
extend the [last optimized model](fresh_weighted_raw110_full_model_20261017.md).
No frozen predecessor is changed.

## The changed finite matrix

Twelve proved packages export213 selected actual columns since the previous
solve. They include historical weighted corners, new physical-pair exclusions,
raw125/raw205 corners, forest1325 incidence, the raw23/raw25 family, shared
tetrahedra and both raw39 target packages. Hypothetical zero-event columns
are excluded. Their positive face coefficients and nonnegative vertex costs
are taken from their exact saved exports, with the actual face and forest IDs.

A common positive denominator clears every coefficient of each entire
column, including every vertex cost. Dividing by their common integer gcd
gives a primitive column. Four exact positive-scaling duplicates are merged,
retaining all213 origin records. The result has209 distinct new columns:
195 charge one vertex type and14 charge two. No several-vertex cost is
replaced by a zero-budget transfer or assigned to just one vertex.

|Quantity|Count|
|---|---:|
|Inherited variables|439410|
|Distinct added variables|209|
|Total variables|439619|
|Active variables|173793|
|Face obligations|82806|
|Vertex budgets|6558|
|Inherited zero-budget signed variables|96176|
|Inherited positive-cost signed variables|52|
|Face-matrix nonzeros|828714|

Every old variable record, face coefficient, face requirement, vertex-budget
coefficient and zero-budget record is compared with the frozen reconstruction.
The provenance-to-face map stays fixed. The matrix digest is
`c2196971b65d5ff951cb36297076353be4fd059fa7ea5416b4f509d1d2b28758`.

All input bindings are checked. All twelve packages use the same former
5053 positive face masses. Their213 columns also satisfy the latest exact
full fixed-mass repair, which has ratio below12. That check explains why a
new optimization is warranted; those old local covering mixtures are not
assigned to the newly optimized face masses.

## Exact recovery with several-vertex costs

For a column j, let A_fj be its face coefficients and B_vj its actual
vertex-budget coefficients. Every variable is nonnegative, including
variables with signed face entries. A dual certificate must satisfy

    sum_f alpha_f A_fj <= sum_v mu_v B_vj

for every column, with nonnegative alpha and mu. The recovery and verifier
construct B directly from all6558 saved budget rows. They do not infer it
from a single-owner metadata field. Exactly the96176 recorded zero-budget
columns have empty cost rows; all14 new several-vertex columns have both
costs present.

Recovery rounds the solver's face and vertex multipliers to rational numbers.
It applies the inherited simplex-mass correction for zero-budget signed
transfers. If a remaining column has positive deficit, it increases one
of its actual vertex allowances by the exact required amount, choosing a
largest positive coefficient. Because every vertex cost is nonnegative,
this cannot invalidate any previously satisfied column. The verifier checks
all439619 final column inequalities, exact numerator and normalizer, and
the target flag. The saved certificate has1250 allowance-repair records.

Primal recovery uses the inherited strictly positive old repair columns.
It rationalizes the proposed weights and adds a common nonnegative amount
only on those old columns until every face obligation holds. It then
recomputes every full vertex budget, including all costs of the new columns.
The unchanged inherited recovery already reads the budget matrix, so no
single-owner assumption enters this part either.

The guarded solve counter permits one HiGHS interior-point solve, with
a300-second limit. It terminates optimally in83.7995 seconds after106
iterations, at numerical objective12.35985517100107. Exact recovery yields

    lower =205569148850452800/16632023034330289,
    upper =6179927585791/500000000000.

The upper certificate has23645 nonzero rounded variable numerators and
repair delta3/500000000000. The lower certificate has5121 positive face
masses and1050 positive vertex allowances. Both endpoints exceed12.

## Verification and limits

Preparation9636da and the single solve872707 terminate exit0/PASS. Separate
full-matrix replay749374 also terminates exit0/PASS. It rebuilds the entire
changed matrix and checks every exact primal face/budget inequality and
dual column inequality without invoking the optimizer.

The independent stdlib saved checker3389d5 terminates exit0/PASS. It imports
neither producer nor recovery engine. It checks41 direct bindings, all213
export origins, common positive scaling and primitive deduplication, every
new dual load against all its vertex costs, nonnegative saved weights and
the new primal contributions, exact ratios, interval order and target flags.
It does not reconstruct the old matrix or independently verify its face
requirements and old-column inequalities. Those remain the full replay's
explicit scope.

A separate peer completed a full producer/helper source read, including the
inherited primal recovery and actual-versus-hypothetical export interfaces.
It checked matrix preservation, cost scaling, all-owner budget treatment,
monotone exact dual repair, full verification and the helper's limited scope.
That source review reran no matrix construction, optimization or certificate
arithmetic. The peer's full final report read also passed. Its small saved-data
display check4caa94 confirms the stated support counts and that the displayed
lower decimal rounds outward below the exact fraction; the upper decimal is
exact. Source, data, report and separate checker are final/frozen.

The matrix includes the213 stated exports, not every charge in the complete
implicit local packing cones. Old fixed-alpha covers cannot be reused at
these new face masses without reconstruction. Missing face, finite-forest,
ground-field and charge coverage also remain. Thus an eventual target12
primal in this finite model would still require further mathematical work
before a dimension theorem. Its present dual is an obstruction to this
finite charge family, not a simultaneously realizable hyperbolic polytope.
No external refereeing or proof-assistant verification is claimed.
