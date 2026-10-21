# Lorentz obstructions: internal logic and application boundaries

`formal/CoxeterProofs/LorentzObstructionKernel.lean` has an actual bounded
compile PASS and a separate actual axiom-audit PASS for final source
`5c0f8e57c9a1fb2dd1ca740d2424a5359169cde199c351d782188fe328c0a6d1`.
It imports only `Std`. All seven requested declarations were actually
printed; only `propext`, `Classical.choice` and `Quot.sound` occur, with
no `sorryAx`. The immutable
[verification manifest](../results/lean_lorentz_obstruction_verification_20261021/manifest.json)
binds the source, compiler, object ledger, actual output and two separate
successful runs. This verifies conditional internal algebra, not the
geometric application inputs below.

## A second-orbit obstruction

For a positive squared root norm c, let q be the squared coefficient
along a fixed future unit timelike vector and W=q+c the Euclidean squared norm
of its orthogonal projection. If an isometry fixes that timelike vector
and the first-orbit product is zero, Euclidean Cauchy gives

```
2*q² <= W*(W + second + q).
```

The kernel derives `second > -3*c` from this inequality. If instead
`second <= -3*c`, multiplication by positive W bounds the right side by
`(q+c)*(2*q-2*c) = 2*q² - 2*c²`, a contradiction. It specializes to
the norm-two cutoff `second > -6` and proves that the explicit golden
pair `second = -80 - 132*phi`, with `phi > 1`, lies below that cutoff.

The finite-order application has a named, substantive external bridge:
averaging a future timelike orbit supplies a fixed future unit vector;
the Lorentz/Euclidean decomposition and Cauchy supply q,W and the displayed
inequality. The kernel assumes that bridge as a theorem parameter, then
derives that finite order is impossible. It does not assume an infinite
order conclusion or declare a custom axiom.

Actual original-lattice integrality, the map of full chamber germs and
the index/covolume consequence remain separate written inputs. The new
norm79 D9 application has separate [written acceptance and exact-map
checks](../results/a14_a18_and_prime79_period_acceptance_20261021.json).
The verified conditional kernel alone excludes no lattice or polytope.

## The marked-D quadrilateral obstruction

The kernel expands the cleared two-by-two determinant numerator and proves
it cannot vanish when `m >= 3`, `h >= 0`, and `g = 0` or `g >= 1`.
At zero g the numerator is negative. Otherwise it is

```
h*((m+1)*g-2) + 2*m*g - 4 > 0.
```

It also checks the natural-number contradiction between rank n+2 and
the ambient upper rank n+1. The complete white support, the actual
long-end/center role assignment, native gap, Schur elimination and rank
interpretation remain written mathematics in
[the restricted theorem](../results/finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md).
This is not a universal transfer rule or a reduction of finite57's charge372.

## Actual verification and reproduction

The root-owned compile took0.52s and the separate axiom audit0.51s, with
initial free disk134.56GiB and available RAM7.54/7.57GiB respectively.
These are admission observations, not peak-memory measurements. One
unused-section-variable warning in `schur_numerator_identity` is retained.
The actual seven declaration names are:

* `second_orbit_gap`;
* `norm_two_second_orbit_gap`;
* `golden_pair_beyond_cutoff`;
* `no_finite_orbit_data`;
* `schur_numerator_identity`;
* `large_d_schur_numerator_nonzero`;
* `full_white_rank_contradiction`.

The initial tactic failure and the subsequent unqualified-name failure
are byte-preserved in
[the first preservation receipt](../results/lorentz_obstruction_initial_failure_preservation_20261021.json)
and [the second preservation receipt](../results/lorentz_obstruction_transitivity_name_failure_preservation_20261021.json).
Both receive ZERO verification credit, including any recovery prints with
`sorryAx`. The final source explicitly proves square positivity, rewrites
the Cauchy premise and composes inequalities through the ordered-preorder
class field. The seven statements and their external inputs are unchanged.

The successful commands, run separately and serially from the repository
root under the existing resource gates, were:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.LorentzObstructionKernel
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.LorentzObstructionKernel --audit-axioms
```

Do not rerun the successful module or dependencies without a substantive
source change or unresolved issue. The runner retained one module,
2,048 MiB, one thread, the shared exclusive lock and all existing
admission/abort gates; no successful dependency was rerun. Actual printed
axioms are distinct from the explicit geometric and arithmetic theorem
parameters above. Source and delta reviews are
[the full static review](../results/review_lorentz_obstruction_kernel_source_foundations_20261021.md)
and [the final source binding](../results/review_lorentz_obstruction_preorder_name_precision_delta_foundations_20261021.md);
the latter corrects the intermediate transitivity-name identification.
