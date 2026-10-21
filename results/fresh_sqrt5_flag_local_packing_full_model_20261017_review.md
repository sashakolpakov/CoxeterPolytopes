# Independent source review of the flag and local-packing successor

The mathematical/source review of the
[combined model](../scripts/fresh_sqrt5_flag_local_packing_full_model_20261017.py)
passes. No gap was found in the column construction, preservation of the
inherited domain, or exact certificate recovery. This is a narrow source
review: it does not repeat the large LP, independently rebuild the full
matrix, or certify a particular numerical solver outcome.

The reviewed source SHA256 is
`8cad3c6a3a688171cc7cb42804c8bdf2fcb22b5c6bd233fc980217676701950a`.

The model retains all82,806 inherited face obligations,350,744 old
variables,6,558 vertex budgets and17,090 signed columns. Its explicit
submatrix checks preserve the old coefficients, face kinds, demands and
indices, including every actual-provenance copy of a split quadrilateral
row. Adding features does not discard the old unexported rows.

The77,391 new [flag transfers](fresh_sqrt5_flag_refined_base_capacity_20261017.md)
have positive actual prism incidences and one negative containing-simplex
capacity. The inherited provenance map identifies their1,858 possible
actual quadrilateral rows. The signs agree with the geometric inequality
`sum(prism flag incidences) <= capacity * simplex count`; their global
cost is nonpositive. They have zero vertex cost. The underlying ordered
flag and spectator multiplicities were separately reviewed in the flag
package.

The175 positive columns are exactly148 saved all-subset charges,4
intersecting-orbit clique charges and23 local-packing charges. Each has
nonnegative simplex coefficients and a positive vertex cost. Multiplying
all coefficients and its cost by the same positive common denominator
only reparametrizes its nonnegative variable. The source asserts the
scaling identities exactly. The subset package's inherited model hash
matches this model, and the three cut packages use the same original
face masses and row indices. These are finite selected cuts; the source
does not substitute the entire local packing cone.

New positive coefficients are supplied only on the exported simplex
classes, with zero quadrilateral coefficients. Flag positive coefficients
are supplied only on the selected actual prism provenances. Omitted
positive payments weaken the rule, which is valid; this does not assert
that the omitted geometric incidences vanish. The negative simplex
capacity terms are retained. This distinction is explicit in the saved
scope.

The active-column restriction is safe. Any omitted column has no positive
face coefficient; setting it to zero cannot reduce a face payment and
cannot increase a nonnegative vertex cost. Every face has at least one
retained inherited positive repair variable.

For a primal certificate, rounding leaves nonnegative rational variables.
The exact common repair increment is chosen from every face deficit and
added only to those inherited positive variables. The verifier checks
every face and vertex budget, the exact maximum, valid variable indices,
and the target flag. Signed-column rounding is therefore covered by the
same exact deficit check, rather than presumed harmless.

For a dual certificate, all face and vertex masses are checked positive
and correctly indexed. Each signed constraint is repaired by increasing
its simplex mass. All signed positive entries lie on quadrilateral rows,
so increasing simplex masses cannot worsen another signed inequality;
one pass suffices. The exact routine then recomputes all positive-column
loads and vertex masses. It also checks the claimed ratio and target
flag. Saved replay reconstructs the model hash and exact certificates
without a new optimizer invocation.

The preparation/solve/replay protocol records the sole bounded optimizer
call separately from the exact certificates. A solver success alone is
not a mathematical certificate. Conversely, even an exact primal at the
target would apply only to this preserved partial geometric domain; it
would not by itself be an ambient dimension theorem.

The source was read independently of its author; inherited geometric
lemmas and frozen charge certificates remain explicit inputs. No source
changes were requested. **FINAL/FROZEN.**
