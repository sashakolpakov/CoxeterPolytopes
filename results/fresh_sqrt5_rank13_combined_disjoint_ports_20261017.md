# Combined rank-13 model with disjoint simplex ports

This successor adds the newly proved one-port simplex incidence charges to
the complete saved combined model. It is a mathematical enlargement of the
payment family, not a repeat of the earlier timed-out optimization.

The [source](../scripts/fresh_sqrt5_rank13_combined_disjoint_ports_20261017.py)
and [saved artifact](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.json)
bind the frozen
[combined model](fresh_sqrt5_rank13_simplex_quad_combined_20261017.md) and
[disjoint-port theorem/export](fresh_simplex_disjoint_ports_20261017.md).

## New columns and preserved domain

Every inherited variable, face coefficient and budget entry remains. The
model still contains 6,558 complete 78-pair vertex budgets, all 10,185
flexible simplex aggregate rows, and 68,150 distinct inequalities representing
all 81,183 canonical flexible dotted quadrilateral stars. The 22 older
simplex packing charges also remain.

Distinct maximal simplex faces through a simple vertex use disjoint port
sets. For a labelled diagram port orbit O, one new nonnegative coefficient
pays the number of ports in O used at a simplex corner. Its vertex cost is
|O|. This is the integer normalization of an averaged budget-one charge.
It pays an actual simplex flag once; it is not repeated for each triangle
inside that simplex. Every quadrilateral coefficient is zero.

The implementation uses the verified maps from the original 1,870 forests
and 10,185 rows into the combined model. It independently reconstructs all
13 one-marked-root codes at each old forest, checks each exported literal
port orbit, and checks that every simplex row's new coefficients sum to
d(d+1). All 6,528 active exported columns are added, with 54,749 nonzero
face coefficients. The 8,591 exported columns absent from every face row
are omitted: their only effect would be to consume nonnegative budget.

There are 210,115 variables in the exact domain, including 203,565 pair
variables, 22 earlier packing charges and 6,528 new port charges. Numerical
optimization uses the 85,179 columns occurring in at least one face row;
the 124,936 other inherited columns remain explicit and are fixed to zero
numerically. No mixed simplex/quadrilateral incidence columns are added.

## Outcome and scope

The sole `highs-ipm` solve completed within its 300-second limit. Exact
rational recovery and verification give the following interval for the
minimum budget of this partial model:

`46073614230483840 / 3603600104368189 <= optimum <= 95890805370 / 7499999999`.

The endpoints are approximately 12.7854403641 and 12.7854407177. In particular,
the exact lower bound exceeds the target 12. The stronger port charges do
not suffice for this combined domain. The recovered full dual is preserved
without support minimization: 2,397 positive face rows and 620 positive vertex
budgets. Its face support comprises 815 simplex aggregates (212 triangles,
259 tetrahedra and 344 four-dimensional simplexes) and 1,582 quad aggregates.

The largest vertex masses include `2A1+A4+H3+H4`,
`A1+A3+H3+H4+I2(5)`, and `2A1+A3+A4+H4`. The leading simplex row is 599
(saved star 681). One high-mass quad row is 26448, with demand 10 and
source context 27469 / canonical core 19. These are exact diagnostic
identifiers, not a statement that the displayed configurations can coexist.

Every merged quad inequality retains all its source context/core pairs.
Equal pair-orbit coefficients do not imply equal coefficients for a new
mixed-incidence charge. Such a successor must check every represented star,
split the row by the new coefficients, or use a coefficient-wise lower bound
valid for every provenance. The first provenance alone is insufficient.

Since the new columns only enlarge the old payment family while retaining
every old constraint, the same exact lower bound also obstructs the preceding
combined family. Its frozen timeout artifact remains unchanged. Generation
and a fresh saved replay both passed, including complete domain reconstruction
and exact verification of both rational certificates. No further optimization
is being run.

The same 231 untouched finite forests remain outside the budget domain.
Fixed-determinant stars, other connected cores, all-finite quad stars and
other missing face domains are not covered. Exact feasibility would therefore
remain a partial weighting result; a dual would obstruct only this enlarged
family on the listed local configurations. Neither result is a dimension
theorem or an assertion that the local configurations coexist in a compact
ambient polytope.
