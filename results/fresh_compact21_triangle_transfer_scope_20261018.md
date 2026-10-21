# Why the new E7 triangle constraint alone does not close compact21

The new exceptional triangle argument transfers to a compact H21 vertex
of type E7³: at most three of its nine designated faces are triangles.
However, even an arbitrary improvement confined to the **vertex allowance
at E7³** leaves one frozen exact golden-field model dual strictly above40:

\[
\frac{4119593261}{102300000}
=40.2697288465\ldots>40.
\]

The dual has zero multiplier on the E7³ forest. This is a precise
limitation on that proposed use of the new theorem. It does not rule out
new incidence columns, changed face-payment rules, or geometric arguments
that also affect other rows. No improved arithmetic bound is claimed.

The [producer](../scripts/fresh_compact21_triangle_transfer_scope_20261018.py)
and [certificate](fresh_compact21_triangle_transfer_scope_20261018.json)
recheck the saved rational dual directly, without an optimizer or a new
geometry census. They bind the final
[universal triangle packet](fresh_compact29_three_triangle_star_constraint_20261018.md)
and the frozen
[E7/E8-context golden-field model](fresh_sqrt5_dimension21_e7e8_context_20261016.md).

## What transfers geometrically

The triangle proof needs a positive vertex span made entirely of E7 and
E8 components. Its untouched-E7 argument replaces an entire E7 by E8,
which would increase the positive rank by one; it does not depend on the
specific ambient number29. The complete seven-contact domain and native
gap are likewise unchanged. The only rank21 forest consisting entirely
of E7/E8 is E7³, since 7a+8b=21 has the sole nonnegative solution (3,0).
Thus the same proof gives the stated nine-face restriction at that type.

This observation does not extend the proof to E7 components beside D5,
A1, or other finite factors. A third side may have additional contacts
with such factors, changing its projection and height. No exclusion of
those contacts has been supplied here.

The general triangle theorem makes no arithmetic assumption. Its unique
triangle root has a label5 contact, whose cyclic product is
(3+√5)/2. Thus the same triangle is not available when the actual
arithmetic ground field has degree3 and does not contain √5. This is
consistent with the already proved
[cubic-seven bound20](fresh_compact_global_bound_20261015.md); it gives no
new bound for that field. The
[complete compact21 field reduction](fresh_compact21_single_field_20261016.md)
already leaves only Q(√5). Its historical “candidate” wording records
external-validation status; no new unresolved step in that completed
field reduction is asserted here.

## Exact surviving dual

The frozen model has nonnegative bonus variables and a common upper
variable z. Its certificate lists inequalities A_i x≤b_i and nonpositive
dual coefficients y_i. The new check recomputes

\[
\sum_i y_i A_{ij}\le c_j,\qquad
\sum_i y_i b_i=\frac{4119593261}{102300000}>40,
\]

where c is zero on bonus coordinates and one on z. The inequalities give
weak duality for nonnegative primal variables. Every saved term and
coefficient is included: 2,490 terms comprise 1,973 face terms and 517
forest terms, touching 5,441 columns. The total forest mass is
999997/1000000. The finite geometry and complete reconstruction
of the original face/forest coefficients remain inherited from the
frozen model and its recorded replay.

All forest terms are checked to have rank21. None is E7³. Therefore any
change solely to the E7³ vertex inequality, including replacing its
allowance or removing that row entirely, leaves every positive-support
dual term and the displayed lower bound unchanged. This is not a claim
that every conceivable use of the triangle theorem changes only that
inequality. In particular, incidence couplings may introduce new columns
on face rows and require a separate dual feasibility test.

The largest supported forests include A1+D5+E7+E8 and D5+2E8, as in the
earlier [nine-model diagnostic](fresh_sqrt5_dimension21_diagnostic_synthesis_20261016.md).
These are concrete structural targets for a further compact21 argument.
The current triangle theorem has not classified their stars.

## Reproduction and next work

```sh
python3 scripts/fresh_compact21_triangle_transfer_scope_20261018.py --verify
```

The source uses only the Python standard library. It rechecks source
identity and literal rational weak duality; it does not rebuild the full
model, optimize, or rerun the old arithmetic queues. Generation a84e58
and separate replay da28d3 both terminated exit0/PASS. The root agent
completed a full mathematical source/report read, covering the all-E7/E8
rank21 transfer and its limited forest domain, the cubic-field cyclic
product observation, literal negative-multiplier weak duality and all
touched coordinates, and zero E7³ support. The review did not execute
the producer, reconstruct the old model coefficients, or replay the
native/arithmetic queues. It passed without correction. This packet is
final and frozen after that review and the completed computations.

A useful next direction is to choose a recurring supported forest such
as D5+2E8 and derive an actual neighboring-star or coface constraint that
changes this dual. E8 isolation is still available, but a finite endpoint
classification alone is insufficient: every old endpoint type has an
allowed finite reextension. Any new condition must use simultaneous
Lorentzian/native compatibility or global incidence and be checked on
the surviving exact witness before another optimizer run.
