# Forced missing base simplexes, including maximal edges

The missing bases of the positive prism masses form a manageable exact
successor domain. Reconstructing **4,195 missing full-star types requires no
new finite vertex forests and no new ordinary pair variables**. The additive
patch needs 1,224 further one-port variables, of which 475 occur in the new
rows. This package does not run an optimizer or establish a new dimension
bound.

The [source](../scripts/fresh_sqrt5_forced_base_simplex_domain_20261017.py)
and [saved exact domain](fresh_sqrt5_forced_base_simplex_domain_20261017.json)
bind the frozen [base-star census](fresh_sqrt5_prism_simplex_base_incidence_20261017.md),
[shared-capacity dual repair](fresh_sqrt5_prism_shared_transfers_20261017.md),
and [combined disjoint-port model](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md).
Every old index, variable, face row and budget is preserved.

## Selection and complete local reconstruction

The shared-repaired dual has positive mass on 1,582 quad inequalities,
representing 1,858 literal full-star provenances. The selection takes every
missing base type occurring in **any** of those provenances. It does not
replace a merged row by its first representative.

| Maximal simplex dimension | Missing full stars | Positive conservative forced mass |
|---|---:|---:|
| 1: maximal edges | 3,460 | 1,904 |
| 2: triangles | 447 | 402 |
| 3: tetrahedra | 205 | 150 |
| 4: four-simplexes | 83 | 81 |
| Total | 4,195 | 2,537 |

The other 1,658 selected stars are retained even though the coefficientwise
minimum over the merged provenances gives them zero load at the current
dual. Their presence in a possible provenance is not confused with a forced
positive lower multiplicity.

There are 275 connected base cores. Exactly 228 are new relative to the
earlier 310-core triangle catalogue, and all 228 have maximal simplex
dimension one. The 735 higher-dimensional missing stars instead use 47
already known connected cores with previously absent spectator forests.
The 147 selected cores containing dotted pairs are all maximal-edge cores.
The other finite edge labels occurring in the selected cores are 3, 5 and 10.

For every selected core the source tests every one-root deletion against the
complete finite component catalogue, then checks equality with its saved
spherical-deletion list. It adjoins the actual orthogonal spectator forest
and reconstructs every finite rank-13 corner. There are 9,496 such corner
occurrences, using 694 existing whole-forest types. At each corner it records
the literal omitted root, the entire simplex port block, all pair-orbit
coefficients and all one-port-orbit coefficients. Complete vertex budgets
are independently recounted as 78 unordered pairs and 13 ports.

The parent prism stars already have exact locally admissible
`Q(sqrt(5))` gauges. Deleting a normal gives the indicated nonsingular
base star in the same admissible quadratic space; the exact parent
arithmetic is inherited, while the finite corner reconstruction is repeated
here. No new ground field is introduced. Each full rank-14 star lies in the
parity-flexible branch of the previous determinant necessary-condition
filter. This is not a proof that every global arithmetic-form class or every
simultaneous incidence pattern is realizable.

## Why maximal-edge rows are legitimate

The [disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md) proves
that triangular adjacency is an equivalence relation on the edge ports at a
simple vertex. Its components generate the maximal simplex faces. Including
isolated singleton components makes these blocks a partition of all ports.
An isolated edge cannot lie in a larger simplex, since every edge of a
simplex of dimension at least two lies in a triangle. Conversely, a
nonisolated port lies in its unique larger simplex block.

Thus the same nonnegative one-port charge can be shared by maximal edges and
higher-dimensional maximal simplexes, with the same orbit-size vertex
budget. For a maximal d-simplex, the curvature demand from its triangular
two-faces is

\[
 2\binom{d+1}{3}.
\]

The pair-corner coefficient sum is `3 C(d+1,3)`, and the one-port coefficient
sum is `d(d+1)`. In particular a maximal edge has demand zero, no pair-corner
coefficients, and its two endpoint port charges. It is an incidence
bookkeeping row, not an extra positive curvature demand.

For a shared transfer parameter `tau >= 0`, an edge row imposes

\[
 \text{two endpoint port payments}-24\tau\ge0,
\]

with the appropriate sum when several transfers use the same edge type.
The coefficient 24 is `(d+1)(13-d)` at `d=1`. The corresponding positive
quad coefficients count all four edges of a square prism, as in the frozen
invariant base convention. This makes a forced edge carry genuine incidence
cost without inventing a triangle or an arbitrary choice of two square sides.

Construction_audit separately confirmed this isolated-edge extension and
zero-demand row semantics in a narrow conceptual review. This is distinct
from an independent arithmetic implementation.

## Additive model interface

The old positive-cost domain has 210,115 variables, 78,335 rows and 6,558
vertex budgets. The patch appends:

- 4,195 face/incidence rows, with 16,276 nonzero positive-cost coefficients;
- 1,224 one-port variables and their additions to the existing budgets;
- 2,848 zero-vertex-cost shared-transfer columns, with 6,127 positive quad
  entries and one negative simplex/edge entry each.

There are no new ordinary pair variables or vertex budgets. The combined
positive-cost domain would have 211,339 variables and 82,530 rows, with
85,841 variables active in a face row. Besides the 475 active new variables,
the new rows activate 187 old variables. The remaining 749 new port variables
are saved as zero-only columns. The full finite-forest universe still has
231 untouched forests; this package does not silently add their budgets.

One-port variables are looked up by the complete whole forest and its
one-marked-root code. An existing variable is reused, and a new variable is
charged to the same complete vertex budget. There is no independent second
port budget for the added stars. All quad one-port coefficients remain zero.
The 22 old clique charges deliberately receive zero on the new rows; this
restriction is valid and avoids borrowing uncomputed new clique features.

For each selected base type, the transfer coefficient is computed from
**every one of the 81,183 provenances of all 68,150 quad rows**, not merely
from positive rows. It is `k` times the minimum multiplicity over every
provenance of that quad row. The negative entry is
`-(d+1)(13-d)` on the corresponding new base row. The frozen represented-base
shared transfers remain unchanged. All other missing base types retain
unknown masses and remain outside this patch.

The main JSON interfaces are:

- `new_variables`, `new_vertices`, `existing_vertex_budget_additions` and
  `touched_vertex_budgets` for index-preserving budget updates;
- `new_faces`, including literal `corners`, separate pair and port
  coefficients, positive source provenances, and the current conservative
  forced mass;
- `new_shared_transfer_columns`, with explicit positive quad and negative
  simplex/edge entries;
- `base_to_new_face`, `selected_base_types` and `core_profiles`, referring to
  the complete frozen core and provenance data through hashed inputs.

Ordinary pair variables in a future patched domain should be identified by
the inherited pair interval together with `new_pair_variable_indices`, not
by assuming that future additions must precede the existing port variables.
This particular patch adds no pair variables.

## Scope and validation

This is one prepared partial domain, with no optimization and no claimed
successor primal or dual. The current shared-dual masses determine the
selection and the saved necessary base masses; they are not asserted to
remain a valid dual after adding these rows and port coefficients.

Fixed determinant classes, other connected cores and spectator choices,
all-finite quadrilateral aggregates, the remaining missing bases, and other
uncovered faces still require separate treatment. Mixed-port or clique
features on the new positive rows must be reconstructed before combining
their earlier diagnostics. Local arithmetic admissibility is not a claim
of a compact global chamber.

**FINAL/FROZEN:** generation and a fresh saved reconstruction replay passed.
The direct source manifests, Python syntax and local report links were also
checked. No optimizer or independent second arithmetic implementation was
used.
