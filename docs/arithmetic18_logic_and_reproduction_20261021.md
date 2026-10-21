# Compact arithmetic dimension eighteen: proof and verification

The new result is dimension at most18 for compact arithmetic first-type
hyperbolic Coxeter polytopes, over every actual admissible field. The
[written dimension proof](../results/compact_arithmetic_dimension18_global_composition_20261021.md)
and its [full internal composition review](../results/review_compact_arithmetic_dimension18_global_composition_full_scope_outside_author_20261021.md)
state the exact inherited and computed inputs. This is a computer-assisted
written proof with full internal reviews; external refereeing and geometric
Lean verification remain separate.

| Class | Current repository ceiling | Change in this batch |
|---|---:|---|
| Compact arithmetic, first type |18|19 to18, all actual admissible fields|
| Arbitrary compact Coxeter |20|unchanged|
| Arbitrary finite-volume Coxeter |57|unchanged|
| Noncompact arithmetic |21|published; unchanged|
| Arithmetic finite-volume union |21|unchanged|
| Right-angled, finite volume |8|internally proved sharp bound; unchanged|
| Right-angled, entirely ideal |4|internally proved sharp bound; unchanged|
| Right-angled, compact |4|published sharp bound; unchanged|

## How the new proof fits together

The earlier all-field theorem already gives ceiling19. At dimension19,
the complete farthest-vertex analysis leaves only a whole D19 diagram.
The new proof rules out its actual outgoing fork edge. It first treats
orthogonal, finite, parallel and isolated dotted opposites. A nonisolated
opposite has one of six possible white contacts. The
[finite-contact proof](../results/arithmetic_d19_nonisolated_finite_fork_exclusion_20261021.md)
removes every finite own label without a height assumption.

For dotted contacts, the
[field reduction](../results/arithmetic_d19_actual_field_is_quadratic_20261021.md)
proves that the actual field must be quadratic. This uses the complete
changed-threshold prefix inventory and the centre-height floor; it is
not a restriction to a selected list of fields. The
[two near-contact proofs](../results/arithmetic_d19_near_dposition_and_four_dotted_fork_exclusion_20261021.md)
then reuse the earlier exact polynomial compatibility data. For the
[other four contacts](../results/arithmetic_d19_remaining_dotted_fork_physical_bounds_and_census_20261021.md),
actual feasible cores give a physical size bound and a complete finite
domain of290 integer rows.

The exact census leaves133 necessary fork models. Testing each against
the complete257-row weak-prefix inventory requires717 same-field square
tests. Three anchors and five pairs survive. The
[final four-contact proof](../results/arithmetic_d19_four_remaining_dotted_fork_exclusions_20261021.md)
handles every survivor with one positive reference vector: average a
few neighboring slopes, so the coefficients of the possibly weak
directions become zero. The energy remains greater than1, contradicting
the original normalization. The output is deliberately not called empty,
and no surviving formal model is asserted to be a polytope.

The dimension composition checks all fork cases, rather than extrapolating
a restricted lemma. It gives18 only for the stated compact first-type
class. The quadratic-field theorem proved along the way concerns the
hypothetical whole-D19 configuration; it must not be imported unchanged
into a new dimension18 analysis.

## Reproduce the new computed proposition

The [mathematical specification](../results/arithmetic_d19_remaining_dotted_fork_physical_bounds_and_census_20261021.md)
separates the domain, algorithm, termination, soundness/completeness and
actual output. The implementation uses rational pairs in quadratic fields,
exact embedding signs and the full field-square test, with a multiplied
witness for every positive result. A square rational norm alone does not
pass the test. The arithmetic and prefix generator are reused read-only
from the earlier fully reviewed worker.

From the repository root, choose NEW output names:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock \
  --ledger results/NEW_arithmetic18_census_bounded.json -- \
  scripts/check_d19_nonisolated_quadratic_fork_time_compatibility_20261021.py \
  --proof results/arithmetic_d19_remaining_dotted_fork_physical_bounds_and_census_20261021.md \
  --proof-sha 80d67e3b48b33534f7e716e2011323123d4b05a4f1aa6fe9fd770919bb9b930f \
  --out results/NEW_arithmetic18_census.json
```

The [recorded output](../results/d19_nonisolated_quadratic_fork_time_compatibility_20261021.json)
and [guarded execution](../results/d19_nonisolated_quadratic_fork_time_compatibility_bounded_20261021.json)
are preserved. This run finished in0.123 seconds with peak process-tree
RSS56,016,896 bytes. The prior isolated-fork census remains a separate
input, with its original reproduction command and one-survivor proof.
Replaying only the new worker is not a full transitive proof regeneration.

## What Lean now checks

The new [OmittedPrefixReferenceKernel](../formal/CoxeterProofs/OmittedPrefixReferenceKernel.lean)
compiled and received an actual axiom-output audit. Its eight theorems
prove the shared internal mechanism:

| Logic | Inputs supplied by an application | What is derived |
|---|---|---|
| Supported reference pressure and contradiction | Nonnegative finite weights; coordinate directions; pressure floors only at positive weights; energy>1 and original norm<1 | Weighted pressure, coordinate energy lower bound and contradiction; zero-weight directions need no floor |
| Two- and three-slope averaging | The equation defining the arithmetic mean | Exact cleared energy-debit identities and nonincrease, without an assumed variance identity |
| Nonflat graft | Coordinate vectors and a scalar coefficient | Complete square expansion, retaining the overlap with the original reference |

This does NOT formalize actual hyperbolic feasibility, the finite-white
classification, field normalization, the computed census, radical bounds,
the concrete knot/reference construction or the complete dimension18
theorem. Those remain explicit written or computed mathematical inputs
outside this selected kernel. There is no new global-bound axiom wrapper.

Immutable copies of the [compile record](../results/arithmetic18_omitted_reference_compile_20261021.json),
[axiom audit](../results/arithmetic18_omitted_reference_axiom_audit_20261021.json)
and [actual output](../results/d19_omitted_prefix_reference_guarded_lean_stdout_20261021.txt)
bind the current source, dependency objects and compiler. All eight actual
axiom lists contain only propext, Classical.choice and Quot.sound. Two
unused-section-variable linter warnings are retained; compilation succeeded.

For a sequential cached check and fresh axiom audit:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py \
  CoxeterProofs.OmittedPrefixReferenceKernel --audit-axioms
```

The runner uses one module/thread and a2048-MiB Lean cap, preserving
the required dependency objects and cleaning disposable products. Both
the scientific and Lean runners retain disk/memory gates. Below20GiB
free disk, start cleanup before another job. No unrestricted build is needed.

## Papers and evidence

There are still [six paper families](../manuscripts/README.md). The frozen
arithmetic PDF prints20; the stronger18 proof is a successor in that SAME
family and has not been silently substituted into its PDF or older Lean
kernel. Parked refinements, original certificates and source/PDF pairs
remain preserved. The
[new acceptance record](../results/compact_arithmetic18_written_acceptance_20261021.json)
binds the source proofs, full reviews, exact outputs, selected Lean records
and predecessor archives at their precise verification scopes.
