# Combined flexible simplex and dotted quadrilateral model in rank 13

This is one bounded optimization of the combined saved face domain. It retains
the 22 proved simplex packing charges and adds every flexible canonical dotted
quadrilateral aggregate from the completed separation. Mixed simplex/quad
incidence columns are deliberately absent from this first test.

The [source](../scripts/fresh_sqrt5_rank13_simplex_quad_combined_20261017.py)
reconstructs its exact domain from the frozen
[simplex rule](fresh_sqrt5_rank13_simplex_packing_repair_20261017.md),
[quad separation](fresh_sqrt5_rank13_packing_quad_separation_20261017.md), and
[spectator catalogue](fresh_sqrt5_rank13_dotted_b1_spectators_20261017.md).
The [saved artifact](fresh_sqrt5_rank13_simplex_quad_combined_20261017.json)
uses compact JSON and is the sole domain export.

## Exact model

The model has 6,558 finite rank-13 forests, including all 1,870 previous
forests and 4,688 newly touched forests. It reconstructs every unordered
root pair at each forest: exactly 78 pairs, grouped into 203,565 labelled
diagram orbits. The 22 simplex packing variables are additional budget
charges; they are not counted as root pairs.

Every old variable has an explicit injective map into the new domain.
Its forest, orbit multiplicity and literal pair list are checked against the
new enumeration. Every old face row and every old packing coefficient is
preserved. Packing columns retain their denominator scaling in both the face
coefficients and the destination vertex budget.

All 10,185 previous flexible simplex aggregate rows remain. The 81,183
canonical flexible dotted quad stars produce 68,150 distinct aggregate
inequalities after deduplicating equal demands and coefficients. Source
context/core indices preserve every occurrence. Thus the model has 78,335
face rows and 203,587 variables in total.

For a quad star, each certified `K2,2` contributes its four literal corner
pairs, and the demand is the number of quadrilaterals. Orthogonal spectators
are added to the whole corner forest and its marked code before looking up
the corresponding variable. Every quad coefficient refers to an ordinary
pair variable: the simplex-specific charges contribute zero to quads.

The complete budget domain retains all variables. The sole numerical solve
fixes the 124,936 columns absent from every face row to zero, leaving 78,651
active columns. This cannot weaken the optimum: these nonnegative columns
only consume vertex budget. The optimizer has a 300-second limit and is not
restarted on a timeout.

## Outcome and limits

The sole 300-second optimization reached its time limit. HiGHS returned no
primal solution, objective or dual multipliers, so there is no rational
certificate to recover. The saved status is `UNRESOLVED_COMBINED_PARTIAL_MODEL`.
The run was not restarted. A fresh saved replay rebuilt the whole domain
without invoking the optimizer and passed exact saved-data equality.

The 231 finite rank-13 forests outside the touched domain remain explicit in
the saved artifact. Fixed-determinant stars, missing connected simplex cores,
all-finite quadrilateral stars and other missing face domains are not covered.
A feasible rule would still be a partial model result. A dual obstruction
would exclude only this pair-orbit plus 22-simplex-charge family on the
listed local arithmetic configurations. Neither outcome is a dimension
bound or a claim that these configurations occur together in a compact
ambient polytope.
