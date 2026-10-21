# Multi-anchor payment for one paid-foot energy bound

Source: [MultiAnchorPaidFootKernel.lean](../formal/CoxeterProofs/MultiAnchorPaidFootKernel.lean),
SHA256 `801ccfc1dd19afee9ab65451ddc63640d95efc7a4c0d802d6433dfe4ff897d4f`.
**Bounded compilation and six actual axiom outputs pass.** The unchanged
module compiled in2.19s and audited in.38s, one thread/2048MiB, with
all five dependencies retaining capped provenance. The immutable
[compile](../results/multi_anchor_paid_foot_compile_20261021.json) and
[audit](../results/multi_anchor_paid_foot_axiom_audit_20261021.json) records
both bind the same772-byte [actual output](../results/multi_anchor_paid_foot_actual_axioms_20261021.txt).
Its [parsed summary](../results/multi_anchor_paid_foot_axiom_summary_20261021.json)
contains exactly the six declarations below, each using only propext,
Classical.choice and Quot.sound, with no sorryAx. This adds six selected
conditional declarations: twelve selected kernels/99 declarations overall,
not a completed formal geometric dimension theorem.

The only direct import is the already verified
[AnchorAndPaidFootKernel](../formal/CoxeterProofs/AnchorAndPaidFootKernel.lean),
SHA256 `ac811bf2cd70cc975ccbbcd446d0a5d575311e5506e7141a967d150e7a639457`.
Its finite-coordinate machinery is inherited from
[SummedActualFootKernel](../formal/CoxeterProofs/SummedActualFootKernel.lean),
SHA256 `41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28`.
Neither existing source was changed.

The written target is the coefficient-payment step in the
[A20 two-weak-end supporting argument](../results/compact_whole_a20_two_weak_ends_paid_feet_candidate_20261021.md),
source `bf61078630da718816a668c5df4164b7210dd489ce19c39558a926bdd1f712e0`.
It also supplies the payment mechanism for the new
[one-weak-end/whole-A20 argument](../results/compact_whole_a20_one_weak_end_and_global_family_exclusion_20261021.md).
There n20,d21 and the two nonzero a coefficients are5gamma/6 at Omega6
and1/10 at Omega16, with mu7/50. Numbering the Lean array from zero
puts those entries at5 and15. The written geometry proves capacities,
flatness and residual square; no concrete Lean instantiation of those
geometric facts is claimed by compiling the general module.
For any finite family of old anchor vectors Omega_i, set

    w = sum beta_i Omega_i, P = sum a_i Omega_i,
    b_i = beta_i - mu*a_i, H = sum a_i*(w.Omega_i).

The primitive payment inputs are mu>=0, a_i>=0 and mu*a_i<=beta_i
on the finite range. They imply both beta_i>=0 and b_i>=0, and positive
residual support implies positive original support. Thus original floors
`x.Omega_i>=w.Omega_i` are requested only where beta_i>0. No lower bound
is subtracted from an unknown height. Finite sums prove the exact
coordinate identity `w=sum b_i Omega_i+mu*P` and the pairing `w.P=H`.
Neither identity is supplied as a theorem hypothesis.

| Declaration | Derived content |
|---|---|
| `payment_coefficients` | Old/residual nonnegativity and residual-support inclusion from primitive capacities |
| `weighted_payment_decomposition` | Exact payment identity at every coordinate |
| `joint_anchor_height` | Composite payment height from the individual anchor pairings |
| `multi_anchor_paid_pressure` | Paid-vector pressure from original supported floors and one foot inequality |
| `multi_anchor_paid_energy` | One norm gain by instantiating existing `paid_reference_energy`, with scalar multiplier1 and the derived composite P/H |
| `multi_anchor_paid_contradiction` | That derived lower bound above1 contradicts the SAME x having squared norm below1 |

The scalar type is an arbitrary linearly ordered commutative ring, as in
the imported modules. Nat appears only in finite index/coordinate bounds;
it does not replace real heights, coefficients, pairings or norms.

Applications still supply actual coordinate realization and available
old floors, `g<=x.V`, flatness `w.(V-P)=0`, and the exact residual square
`||V-P||^2=S`. The reference w is DEFINED by the primitive beta weights;
its energy and H are actual coordinate expressions, not abstract supplied
pressure fields. A concrete strict scalar comparison and `||x||^2<1`
remain inputs to the final contradiction. If V is a sum of several actual
feet, its joint pressure, flatness and residual-square bound must first
be established for that SAME vector. The module does not independently
pay several copies of the same old coefficient or add separate norm gains.

Complete cap/contact inventories, feasible auxiliary hyperbolic points,
native scalar estimates, census correctness, actual farthest vertices and
the A20 or global dimension conclusions remain written application
inputs. This is a conditional payment mechanism, not a new geometric or
dimension theorem. The complete [pre-execution source review](../results/review_multi_anchor_paid_foot_source_replay_reviewer_20261021.md)
binds the source and original preparation map. That map's exactb909940f
bytes are [preserved](../results/source_archive/multi_anchor_paid_foot_map_before_actual_compilation_20261021.md.gz)
with gzip round-trip verification. Source preparation and actual execution
retain their separate evidence. No successful dependency or exact search
was repeated; the batch runner removed its temporary products.
