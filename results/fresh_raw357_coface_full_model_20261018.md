# Sixteen coface transfers leave an exact lower bound above12

The sixteen existing-row transfers from the frozen
[raw357 coface theorem](fresh_raw357_simplex_coface_transfer_20261018.md)
are now installed in a new copy of the complete current exported matrix.
An exact feasible dual survives with lower bound

    4272824056360160 / 346500030107667 = 12.331381486554207…

The inherited exact primal remains feasible, giving upper bound

    1233971157349 / 100000000000 = 12.33971157349.

Therefore these sixteen transfers cannot lower the optimum of this finite
model to12, even when all its face masses are free. No numerical optimization
was needed to establish that obstruction. The lower certificate is a chosen
feasible assignment, not an optimized assignment or an exact optimum.

The [producer](../scripts/fresh_raw357_coface_full_model_20261018.py),
[certificate](fresh_raw357_coface_full_model_20261018.json), and
[saved checker](../scripts/check_raw357_coface_full_model_saved_20261018.py)
keep this model conclusion separate from geometric completeness. A feasible
dual above12 is not a hyperbolic chamber, and no compact-arithmetic H13
exclusion or new dimension bound follows.

## Exact changed matrix

The producer reconstructs the frozen
[all-marker model](fresh_all_marker_families_full_model_20261018.md), including
all prior ordinary charges and raw56/57/76/108/112 marker exports. It compares
the complete prepared record with the frozen normalized-dual continuation.
It then appends the sixteen coface columns, each of the form

    +1 on the raw357 source Q; −2b on its actual simplex S; vertex cost0.

These are nonnegative transfer variables with signed face coefficients,
not equalities or unrestricted variables. The simplex requirement remains20.
The remaining31 simplex types absent from the inherited matrix are not
installed, given artificial zero multipliers, or silently treated as present.

The new matrix has442226 variables,82806 face rows,6558 vertex-budget rows,
176400 active variables,96192 zero-budget signed variables and833711 face
coefficients. Every old442210 variable,96176 zero-budget record, face
coefficient and requirement, and vertex-budget coefficient is preserved
literally. The52 old positive-budget signed variables and all14 multi-owner
columns retain their complete costs. No unknown physical compatibility or
historical owner assignment changes.

The exact matrix hash is
`367da8a42e7b8f8cfd309f345486a7b7429f668657fcfdd60a57591422139948`.
Appending all sixteen valid columns does not assert that every one is new
relative to every old constraint. The frozen coface certificate already
establishes six strict violations by the old exact dual, so this is a
genuinely strengthened model.

## Cheap exact certificate before optimization

The old certificate has positive mass on nine raw357 source rows. Their
required-payment-weighted mass is83284683/10000000000; each source has
requirement1. Set these nine multipliers to zero and retain every other
face multiplier and every vertex multiplier exactly.

This candidate satisfies all sixteen new inequalities: their positive
source terms vanish and their simplex terms are nonpositive. Deleting face
mass could in principle invalidate old signed constraints, so the producer
does not infer feasibility solely from this observation. It reconstructs
and checks every one of the442226 exact column inequalities against the
full vertex-cost matrix. All pass. It also records that the nine changed
face rows have no negative old coefficients, consistent with that result.

The objective numerator decreases by exactly83284683/10000000000, while
the normalizer remains115500010035889/115500000000. The resulting rational
ratio is the lower bound above. The old primal is verified on the same
changed matrix with every new variable zero, giving the stated upper bound.

The model allows arbitrary nonnegative face and vertex multipliers; the
chosen deletion is merely one feasible point. Weak duality makes its value
a lower bound for the optimum over all such assignments. Thus changing
other masses cannot yield a primal bound at most12 for this matrix. This
does not assert that the optimizer would choose the deleted masses, nor
that the improvement from the old optimum is known. The manual lower
certificate is weaker than the previous lower certificate and is used
because it remains feasible after the new constraints.

This gives a concrete next requirement: a further constraint must exclude
this new feasible witness. A larger effect on the old, now infeasible
raw357 masses alone is insufficient. Broader coface restrictions or other
structural constraints should be tested against this surviving witness
before another global solve.

## Verification scope

Generation97599/94d717, separate full-matrix replay87103/8cef97, and the
standalone saved checker b1ed83 all terminated exit0/PASS. No optimizer
has been invoked.

Root completed a full producer/report mathematical accounting read. The
construction reviewer completed a full producer/report/helper read and
reread the inherited loads, cost-row and dual-verification interfaces.
Both reviews passed: literal preservation, signed additions, exact mass
deletion, normalization, complete feasibility obligations and the free-mass
weak-duality conclusion are coherent. Neither reviewer executed the
producer/helper, reconstructed the matrix, or repeated optimizer, geometry
or historical census work. This report is final and frozen.

The producer rebuilds the exact inherited matrix and checks both complete
certificates. It rechecks the bound saved charge arithmetic through the
inherited builders; it does not repeat physical pair tests, native embedding
censuses, clique enumeration or geometry proofs.

The independent checker imports no producer. It checks the complete sixteen
saved transfers, the exact nine-row deletion, unchanged remaining
multipliers, objective and normalizer arithmetic, inherited primal identity,
model counts,14-owner preservation and interval arithmetic. It does not
independently reconstruct the old matrix, establish old signed feasibility,
or verify the primal against all old rows. Those obligations belong to the
separate full producer replay. Historical geometry, catalogue completeness
and the actual coface injection remain the frozen theorem inputs.

    python3 scripts/fresh_raw357_coface_full_model_20261018.py --verify
    python3 scripts/check_raw357_coface_full_model_saved_20261018.py
