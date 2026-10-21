# Sharp finite-edge payments, coface flags and local packing

This successor preserves the [preceding checkpoint](finite_edge_coefficient10_and_refined_arithmetic_20261017.md).
The new unrestricted incidence coefficient is sharp. Complete local packing
certificates identify a limitation of the arithmetic model, and retaining
prism coface flags removes its previous fixed dual obstructions. No new
frontier-dimensional example or lower ambient dimension ceiling is claimed.

The unrestricted finite-volume bound163 retains its complete internally
checked computer-assisted proof, with no identified gap; see the
[proof-status statement](dimension163_proof_status_20261017.md). External
refereeing and proof-assistant formalization have not been performed.
The compact arithmetic21 result is unchanged.

## A sharp unrestricted finite-edge inequality

The [539-case six-simplex frontier](fresh_six_simplex_ideal_frontier_20261017.md)
and [completed theorem](fresh_finite_edge_cusp_tradeoff_28over3_20261017.md)
exclude an actual six-simplex having two finite and five ideal vertices.
Every actual six-simplex with at least two finite vertices has at most
four ideal vertices. The argument chooses a minimum-order indefinite
original-wall diagram, uses the complete quasi-Lanner catalogue, and
handles arbitrarily many further supporting whites by exact Schur
corrections. No cutoff on their rank or coefficient field is assumed.

The last74 configurations reduce to111 single-profile models,31 pairs
and three triples. All651 five-vertex determinant systems for the pairs
are checked. Their equations are inconsistent, force at least six ideals,
or force a forbidden positive offdiagonal product. Nonnegative kernels
resolve the triples. Together with the previous higher-dimensional
simplex bounds this proves

\[
 2(T_0+T_1)\le {28\over3}E_{ff}
   ={14\over3}n v_f-{14\over3}E_{fi}.
\]

A classical H9 Coxeter simplex with seven finite and three ideal vertices
attains equality. Thus28/3 is sharp for this inequality; it is not a new
example. The [independent Fraction review](review_finite_edge_28over3_20261017.md)
rechecks all651 pair systems,111 single-profile envelopes,three triples,
and the sharpness Gram and finite-volume vertex census without importing
producer arithmetic or SymPy. The written minimum-order and arbitrary-white
argument has also received an internal mathematical/source review. These
are complete internally checked lemmas. They do not alone control triangles
with two or three ideal vertices and quadrilaterals sufficiently for a new
ambient ceiling.

The [finite-vertex counterexample](fresh_finite_vertex_t2_port_counterexample_20261017.md)
addresses one tempting extension: two triangles with two ideal vertices
can share a finite-ideal edge at a finite vertex without lying in one
simplex face. The exact six-wall finite-volume H3 example has a pentagonal
face separating the other two ports. Thus the triangular-port equivalence
used for finite-finite edges cannot simply be extended to these triangles.
This refutes that proposed proof step, not a numerical bound on T2.

## Exhausting nonnegative local simplex-port charges

The [all-subset package](fresh_sqrt5_all_subset_port_packing_20261017.md)
reconstructs all corners of the3,574 positive simplex rows in the saved
full optimized arithmetic dual. It separates17,828,080 integer packing
inequalities over the13-port forests and repairs the ratio from12.9082390023
to12.8186145350. Its [independent source review](fresh_sqrt5_all_subset_port_packing_20261017_review.md)
checks the integer identities, complete marked-forest orbits, common
symmetry averaging and exact separation logic.

The [intersecting-orbit clique package](fresh_sqrt5_intersecting_block_cliques_20261017.md)
finds four more violated capacity-one charges, leaving12.8175545624.
The [pointed-simplex refinement](fresh_sqrt5_pointed_simplex_base_capacity_20261017.md)
is entirely redundant on the prepared selected domain: its only new splits
are dominated by existing columns with the same capacity.

After a [fourteen-forest pilot](fresh_sqrt5_local_simplex_packing_cones_20261017.md),
the [complete local packing calculation](fresh_sqrt5_all_local_simplex_packing_cones_20261017.md)
covers all740 contributing forests and6,011 positive block orbits. It
finds23 additional exact cuts. Explicit rational mixtures of disjoint
literal blocks cover every required marginal.623 local values have matching
exact lower/upper certificates;117 have honest exact brackets. No cases
remain open. The repaired ratio is

\[
 {69774201720926976000000000000000
  \over5448173768070467817598921440139}
 =12.80689726342\ldots>12.
\]

The upper covering mixtures suffice to prove this obstruction against
**every nonnegative charge depending only on the finite forest and the
maximal-simplex block orbit**, not merely the finitely many charges found
by optimization. A [separate standard-library checker](review_sqrt5_all_local_packing_cover_20261017.md)
verifies all5,255 mixture terms,23,324 literal blocks and6,011 marginal
inequalities without an optimizer or packing enumeration, and reproduces
the ratio. Neither the abstract mixtures nor the dual are claimed
realizable in one hyperbolic chamber. Additional coface geometry lies
outside this exhausted local charge family.

## Retaining the prism coface removes the old fixed witnesses

The [flag-capacity theorem and exact diagnostic](fresh_sqrt5_flag_refined_base_capacity_20261017.md)
retain both a vertex v of a maximal simplex S and its immediate coface U.
For each fixed flag(v,U), at most one counted maximal prism base contains
v inside U. An invariant flag orbit has capacity its literal cardinality.
The older vertex capacity had summed over all n-d possible cofaces.

The implementation uses joint ordered core automorphism orbits. If the
coface root belongs to a finite component, its capacity includes every
identical spectator component. It reconstructs all14,284 actual base-vertex
flags and proves that summing the77,391 flag columns recovers every prior
pointed column. This uses the same1,858 actual prism provenances throughout.
Unselected positive incidence coefficients remain zero, giving valid weaker
columns on the explicit selected face classes.

There are4,298 violations in each saved witness. Increasing1,767 containing
simplex masses and restoring all inherited constraints gives repaired
ratios10.6043152578 for the exported-only witness and10.7310811081 for the
full optimized-alpha witness. Both are below12. These computations remove
the previous fixed obstructions, not all possible dual obstructions; they
are neither feasible primal rules nor dimension bounds. The theorem,
literal coface identification, spectator multiplicities and source logic
have an independent internal review, and both exact diagnostics passed
fresh saved replay.

The [changed-mass local-cover companion](fresh_sqrt5_flag_incremental_local_cover_20261017.md)
revalidates all inherited packing mixtures and supplies2,922 explicit
singleton-block increments for the changed marginals. Its deliberately
conservative ratio is about7.7991. This only verifies a simultaneous repair;
it is not the minimum repair, an optimized bound or a feasible primal rule.
The old fixed-alpha packing certificates are never silently reused after
the flag correction changes simplex masses.

## The enlarged full model has a different obstruction

The [new full-model solve](fresh_sqrt5_flag_local_packing_full_model_20261017.md)
preserves all82,806 inherited face obligations,350,744 inherited variables,
6,558 vertex budgets and17,090 signed columns. It adds all77,391 flag
transfers and175 explicitly selected positive packing cuts:148 maximizing
subset columns,four clique columns and23 complete-local-cone cuts. The
resulting model has428,310 variables,162,484 active variables and801,967
face-matrix nonzeros. It includes those finite cuts, not the entire
all-subset or local-packing cone.

One bounded solve succeeds in about52 seconds. Exact recovery places this
partial model's optimum in the bracket

\[
 {345782670398463480\over27027032875363841}
 =12.793956036279\ldots
 \ \le\ \mathrm{optimum}\ \le\
 {3198492461931\over250000000000}
 =12.793969847724.
\]

Both endpoints are above12. The new inequalities therefore remove the
previous fixed witnesses while a different valid dual obstruction remains.
No second optimization of the unchanged model is warranted. The
[independent source review](fresh_sqrt5_flag_local_packing_full_model_20261017_review.md)
checks preservation, transfer signs, simultaneous positive-column scaling,
class restriction and exact recovery logic. A fresh saved replay verifies
the exact certificates without an optimizer. This is a complete conclusion
for the explicitly partial model, not an ambient dimension exclusion.

The [new support restrictions](fresh_sqrt5_flag_local_packing_dual_support_20261017.md)
show that completing the unexported flag entries alone cannot remove this
obstruction. The new dual uses2,253 simplex rows,355 selected actual quad
rows and1,820 unsplit quad rows outside the flag-positive domain. Removing
all1,820 outside quad masses and exactly repairing the remaining witness
still gives12.5077183952>12. Restricting further to the legacy exported
mixed/integer classes leaves12.3842565597>12, without restoring any outside
face row. These are simultaneous exact duals for the new finite model.

The next test must strengthen already represented configurations. In
particular, the complete local packing covers were proved at the earlier
fixed face masses, while the new model imports only175 selected cuts.
They do not settle the full local cone at the new restricted witness.
That cone can be tested on its actual support before pursuing additional
joint simplex/prism compatibility. Missing provenance coverage remains
necessary for completeness, but its completion alone is insufficient.

## A more connected construction patch

The [disk bridge investigation](fresh_nonarithmetic_b_disk_bridge_20261017.md)
keeps the frozen64-wall H9 control with198 finite vertices and15 target-face
open ends. After releasing wall52, a conditional complete seven-duty disk
filling would reduce18 open ends to11. The exact combinatorial identity
requires the stated retained incidences, simple boundary and absence of
ideal vertices; it does not establish geometric realizability.

An explicit six-cap abstract disk passes exact Euler, link and incidence
checks. A bounded numerical search finds188 bridge candidates with one old
neighbor and two new framing contacts. None completes a compatible disk;
the tested necessary opposite-side contacts fail. Missing labels, seed
families and other disk topologies remain open. No candidate is promoted
to an exact improved chamber, and the control still has certified infinite
volume.

## Reproduction and archive integrity

All listed producer packages are final and frozen. Their generation and
fresh saved replays passed. The exact full model used one numerical solve;
its saved replay uses no optimizer. Independent arithmetic and source
reviews retain their explicitly stated scopes. These internal checks do
not establish publication priority or constitute external refereeing.

The [checkpoint auditor](../scripts/audit_sharp_finite_edge_flags_and_packing_20261017.py)
and [saved audit](sharp_finite_edge_flags_and_packing_20261017_audit.json)
preserve every frozen source, top artifact and report bound by the preceding
archive, recursively check current manifests including prepared model
sources, and check local report links. They also independently check the
new LP certificate signs, indices, exact ratios and target flags. Archive
integrity does not replace the separate geometric proofs or exact replays.

```sh
python3 scripts/audit_sharp_finite_edge_flags_and_packing_20261017.py --verify
```
