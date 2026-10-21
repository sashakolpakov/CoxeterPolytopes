# Enlarged finite model with weighted and raw110 corner capacities

The enlarged finite dimension13 charge model still has optimum above12.
One bounded solve produced exact primal and dual certificates after349
additional necessary capacities were appended. This is a partial-model
obstruction, not a dimension theorem.

The [source](../scripts/fresh_weighted_raw110_full_model_20261017.py) and
[saved model certificate](fresh_weighted_raw110_full_model_20261017.json)
append the86 exported [weighted orbit capacities](fresh_weighted_actual_corner_capacities_20261017.md)
and263 [raw110 mixed capacities](fresh_raw110_mixed_actual_corners_20261017.md)
to the [frozen preceding model](fresh_raw33_295_mixed_full_model_20261017.md).
The latter package had reduced its fixed witness to11.9202817799. That fixed
repair justifies testing this different finite matrix; it does not imply a
target12 primal certificate or an ambient dimension theorem.

## Matrix and provenance preservation

|Quantity|Count|
|---|---:|
|Inherited variables|439061|
|Distinct added variables|349|
|Exact duplicate new exports|0|
|Total variables|439410|
|Active variables|173584|
|Face obligations|82806|
|Vertex-forest budgets|6558|
|Zero-budget signed variables|96176|
|Positive vertex-cost signed variables|52|
|Face-matrix nonzeros|827823|

The exact prepared matrix digest is
`4cc01a7a85e739a52d537200929b442847ca1a3f42c06d3e3f1aacc8be144968`.
Every inherited face coefficient, requirement, variable, vertex-budget
coefficient and zero-budget record is asserted unchanged. All52 signed
vertex-cost incidence columns remain ordinary positive-cost budget columns.
Every variable is nonnegative; signed face coefficients do not make its
variable free.

New columns use their actual saved face IDs, vertex forest, positive cost
and occurrence coefficients. Exact rational tuples are deduplicated.
A common positive denominator scales the face coefficients and vertex
cost together, preserving the inequality. No orientation multiplicity
replaces an actual whole-face occurrence coefficient. No merged face row
receives a newly chosen provenance.

All recorded input hashes are checked before use. The two source packages
are also checked against the preceding dual's same4,713 face masses.
The full weighted repair starts from that dual, and the raw110 repair starts
from precisely its saved weighted vertex masses. Both stages only increase
vertex allowances, and the final fixed repair is below12. These checks
explain the input chain; no old fixed-mass cover is assigned to newly
optimized masses.

## Solve and exact recovery

The inherited engine uses one HiGHS interior-point solve with a300-second
solver limit. A guarded invocation counter prevents repeating that solve
through this package. The prepared model is rebuilt and compared before
optimization; large duplicate matrix files are not saved.

Primal recovery rationalizes the proposed nonnegative weights and repairs
deficits using old columns that have strictly positive coefficients in every
affected face row. Its final certificate checks every face obligation and
vertex budget. Dual recovery rationalizes nonnegative face masses, raises
simplex masses only where needed by zero-budget signed transfers and then
recomputes every positive vertex-cost load. It includes the52 signed-cost
columns. The new349 loads are additionally recounted directly. A separate
fresh verification rebuilds the whole matrix and checks all saved exact
certificates.

The single solve terminated optimally in90.2202 seconds after60 iterations.
Its numerical objective was12.368578902229354. Exact rational recovery gives

    lower = 55714268997452065 / 4504504671872248,
    upper = 1236857890281 / 100000000000.

Both endpoints exceed12. The dual has5,053 positive face masses and1,025
positive vertex allowances. The upper endpoint is the largest exact saved
primal vertex budget, after nonnegative repair; the lower endpoint is the
exact demand divided by the sum of vertex allowances. The target flags
therefore certify that this finite model does not achieve12.

The preceding optimized lower bound was12.3822483175. Its subsequent fixed
weighted/raw110 repair11.9202817799 was a different dual witness; the new
optimization changes face masses and finds the stronger obstruction above.
No claim transfers the old complete weighted-cone or implicit-family
certificates to these changed masses.

## Scope

This is a finite exported charge model. It preserves every inherited
obligation and includes every stated new export, but it does not contain
all implicit local simplex-packing or joint-port constraints. The86 weighted
exports also do not encode every charge of the complete conservative graph
cone at arbitrary new masses. Its old exact covers depend on old face masses.

An exact dual above12 obstructs this finite charge family. It is not a
simultaneously realizable ambient polytope. A target12 primal, if obtained,
would still require the missing face, forest, field and charge coverage
before any dimension theorem. The fixed ground field remains Q(sqrt5).
No unrestricted or compact dimension bound follows from this package.

Preparation53376, the single solve88251 and separate fresh full-matrix
replay29055 all terminated exit0/PASS. A separate saved-data checker imports
no producer: it verifies all349 new exact loads, nonnegativity, ratio and
target flags, rational scaling and the largest saved primal budget. This
small check does not duplicate the full matrix reconstruction.

Root completed a full append-source and mathematical report read, finding
no gap in sparse preservation, rational scaling/deduplication, the exact
weighted-to-raw110 witness chain or the scope of inherited recovery. That
read did not duplicate the solver or full matrix replay. The recovery
engine remains an explicit frozen input. All10 direct bindings, syntax and
5 local report links pass. Source, data and report are final/frozen. All
frozen predecessors remain unchanged.
