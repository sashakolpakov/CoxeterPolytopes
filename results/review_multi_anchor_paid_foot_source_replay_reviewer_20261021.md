# Multi-anchor paid-foot kernel: outside-author source review

21 October 2026. **PASS for internal mathematical logic and the documented
assumption boundary; precompilation source review only.** No compilation,
Lean import, axiom audit, scientific job or geometric computation ran.
This review does not describe the six new declarations as kernel verified.

| Reviewed file | Bytes | SHA256 |
|---|---:|---|
| formal/CoxeterProofs/MultiAnchorPaidFootKernel.lean |7989|801ccfc1dd19afee9ab65451ddc63640d95efc7a4c0d802d6433dfe4ff897d4f|
| docs/multi_anchor_paid_foot_logic_20261021.md |3779|b909940fff002ed7712d866f65b1fba9131c84b497edc7720bf41f9f4bb4495e|
| formal/CoxeterProofs/AnchorAndPaidFootKernel.lean |7362|ac811bf2cd70cc975ccbbcd446d0a5d575311e5506e7141a967d150e7a639457|
| formal/CoxeterProofs/SummedActualFootKernel.lean |14269|41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28|

All four source pins match. I read the complete new module/map and the
relevant imported pressure, energy, finite-sum and dot-product proofs.
Contribution disclosure: I contributed inherited actual-foot/paid-reference
foundations and previously reviewed the imported kernel. I did not author
this new module or map. This is an outside-author review of the new
derivation and its use of those explicit dependencies, not an independent
authorship claim for shared prior work.

## Primitive payment and the six declarations

Write w=sum beta_i*Omega_i, P=sum a_i*Omega_i, b_i=beta_i-mu*a_i and
H=sum a_i*(w.Omega_i), on the stated finite index and coordinate ranges.
The new source uses an arbitrary linearly ordered commutative ring. Nat
controls finite ranges only; heights, coefficients, pairings and squared
norms remain ring elements. It does not replace real geometric quantities
by rational or natural-number surrogates.

1. `payment_coefficients` starts from mu>=0, a_i>=0 and
   mu*a_i<=beta_i. Their product is nonnegative, so beta_i>=0,
   b_i>=0, and b_i>0 implies beta_i>0. In particular nonnegativity
   and support inclusion are conclusions, not extra assumed records.
2. `weighted_payment_decomposition` proves the coordinate identity
   w=sum b_i*Omega_i+mu*P. Its entrywise ring expansion and the imported
   finite-sum addition/scaling laws establish every coordinate, not just
   the coordinates eventually used by the dot product.
3. `joint_anchor_height` proves w.P=H by finite dot-product linearity.
   No independently supplied composite height or payment identity occurs.
4. `multi_anchor_paid_pressure` uses the derived b_i and support
   inclusion. Original floors are requested only where beta_i>0; hence
   every active residual coefficient has an available floor. The
   specialization of `paid_reference_pressure` has composite anchor P,
   scalar multiplier1, actual reference w, actual energy w.w and the
   derived H. All argument directions match those requirements.
5. `multi_anchor_paid_energy` makes the corresponding specialization of
   `paid_reference_energy`, with Y=V-P. The supplied flatness and exact
   residual-square hypotheses are precisely the imported hypotheses;
   scalar1 simplification introduces no new assertion. The imported
   proof derives the paid-vector identity and its norm from finite
   coordinate algebra and one completed square.
6. `multi_anchor_paid_contradiction` applies that derived lower bound
   against the supplied strict scalar comparison and x.x<1. It does not
   assume the desired lower bound on x.x.

For clarity, the mathematical chain checked in the source is

    E=w.w=w.(sum b_i*Omega_i)+mu*H,
    W=sum b_i*Omega_i+mu*V=w+mu*(V-P),
    x.W>=E+mu*(g-H),
    W.W=E+mu^2*S,
    x.x>=2*x.W-W.W>=E+2*mu*(g-H)-mu^2*S.

The pressure inequality follows by applying the old floors only AFTER
forming the nonnegative remaining coefficients, and by multiplying
g<=x.V by mu>=0. It never subtracts an unknown actual anchor height
from a lower bound. The energy inequality is for ONE paid vector on
the SAME x. Multiple anchors do not authorize summing overlapping
norm gains or spending an old coefficient several times.

## Application boundary and verification status

The map accurately identifies the payment step of the written A20
two-weak-end proof bf610786. In a multiple-foot application, V may be
their sum, but the application must establish its combined pressure and
the flatness and squared norm of its full residual, including cross
terms. These facts are not inferred from individual foot names or
separate norm inequalities by this module.

Actual cap classification, actual coordinate realization, geometric
feasibility and farthest-point comparison, construction of old floors,
reference flatness, the exact residual square S and any subsequent bound
on S remain external application proofs. The final comparison
1<E+2*mu*(g-H)-mu^2*S is also an explicit scalar input. This is a
conditional finite-coordinate payment theorem, not a formal theorem
about A20, a polytope, a field, a census or a global dimension ceiling.

The complete new source contains six corresponding `#print axioms`
commands and no `sorry` or newly declared axiom. Source inspection is
not elaboration or an actual axiom-output audit. A separately admitted
bounded compilation and subsequent audit must establish their execution
status and bind source/object/import/compiler/output evidence. I found
no mathematical or documentary blocker requiring an author-source edit.
