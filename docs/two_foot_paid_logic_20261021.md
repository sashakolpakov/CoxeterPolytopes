# Two actual feet with joint coefficient payment

**Bounded compilation and all five actual axiom outputs PASS.**
[TwoFootPaidKernel.lean](../formal/CoxeterProofs/TwoFootPaidKernel.lean),
SHA256 `ff40c49823d399777d3f93ae869d1ee3beaa33aca616d85487b2fd27ba5f1e65`,
contains five checked declarations and five `#print axioms` commands.
The [compile snapshot](../results/two_foot_paid_compile_20261021.json)
and [separate audit snapshot](../results/two_foot_paid_axiom_audit_20261021.json)
bind the same source, compiler and six local dependency objects. Both
use one thread and a2048MiB cap; the compile took0.81s with6.01GiB
available memory, the audit0.61s with6.45GiB. Their ledgers record the
memory setting and admission resources, not a measured peak RSS.

The retained [actual588-byte stdout](../results/two_foot_paid_actual_axioms_20261021.txt)
matches BOTH ledger hashes. Its [parsed summary](../results/two_foot_paid_axiom_summary_20261021.json)
contains exactly the five requested declarations, each with only
`propext`, `Classical.choice` and `Quot.sound`, and no `sorryAx`.
No existing successful dependency was recompiled and no catalogue changed.

The first compile failed at the final bare tactic; its source and full
failure are [byte-preserved](../results/two_foot_initial_failure_preservation_20261021.json).
The repaired final step explicitly doubles the derived pressure and
uses the algebraic square cancellation. Its [independent delta review](../results/review_two_foot_explicit_pressure_repair_20261021.md)
checks all five statements are unchanged. The earlier failed energy
print earns no verification credit. A later [audit admission refusal](../results/two_foot_paid_axiom_audit_gate_refusal_20261021.json)
at5.7GiB is also preserved; no child started on that attempt.
The original preparation map is preserved in
`results/source_archive/two_foot_paid_logic_before_actual_compilation_20261021.md.gz`.

Its only direct import is the unchanged
[MultiAnchorPaidFootKernel](../formal/CoxeterProofs/MultiAnchorPaidFootKernel.lean),
SHA256 `801ccfc1dd19afee9ab65451ddc63640d95efc7a4c0d802d6433dfe4ff897d4f`.
The proof uses the inherited finite-coordinate dot product, sum linearity,
supported-floor pressure and completed-square inequality. It does not
reuse the imported one-foot energy theorem's flatness premise.

For finitely many coordinate vectors Omega_i, define

    w = sum beta_i Omega_i, E = w.w,
    P = sum a_i Omega_i, Q = sum b_i Omega_i,
    r1 = Z-P, r2 = U-Q, W = w+mu*r1+nu*r2.

The primitive scalar inputs are mu,nu>=0, a_i,b_i>=0 and the JOINT
capacity `mu*a_i+nu*b_i<=beta_i`, only for i<n. The actual original
floors are `x.Omega_i>=w.Omega_i` where beta_i>0. The two supplied
actual-foot comparisons are `x.Z>=h` and `x.U>=t` on that SAME x.
The scalar type is an arbitrary linearly ordered commutative ring;
real heights and coordinates are an application instance. All sums
and coordinate hypotheses have explicit finite bounds n,d.

Joint capacity derives beta_i>=0, the remaining coefficients
`c_i=beta_i-mu*a_i-nu*b_i>=0`, and `c_i>0 -> beta_i>0`.
Finite-sum linearity proves `w=sum c_i Omega_i+mu*P+nu*Q` at
every coordinate. The individually retained floors and the two foot
inequalities then prove

    x.W >= E+mu*(h-w.P)+nu*(t-w.Q).

This is derived pressure, not a hypothesis. In particular, no unknown
height is subtracted from a known lower bound and no coefficient is
paid twice. Both costs enter one capacity inequality even when their
supports overlap.

The second structural identity is the FULL coordinate square

    W.W = E+2mu*w.r1+2nu*w.r2
             +mu^2*r1.r1+2mu*nu*r1.r2+nu^2*r2.r2.

Neither w.r1 nor w.r2 is assumed zero. Combining this identity with
the derived pressure and `x.x>=2x.W-W.W` proves

    x.x >= E+2mu*(h-w.Z)+2nu*(t-w.U)
              -mu^2*r1.r1-2mu*nu*r1.r2-nu^2*r2.r2.

The reference/residual products cancel algebraically, while the mutual
product remains. No combined pressure, payment identity, norm expansion,
norm conclusion or desired energy inequality is supplied as a premise.

| Declaration | Derived content |
|---|---|
| `two_payment_coefficients` | Original/residual signs and support inclusion from primitive joint capacities |
| `two_weighted_payment_decomposition` | Exact two-cost payment identity from finite-sum linearity |
| `two_shift_square` | Both reference products, both residual squares and the full mutual product |
| `two_foot_paid_pressure` | One paid-vector pressure from separately retained old floors and two actual feet |
| `two_foot_paid_energy` | One completed-square bound for the same vector, with no flatness premise |

The written target is the two nonorthogonal payments in
[the whole-B20 supporting proof](../results/compact_whole_b20_two_end_actual_foot_exclusion_20261021.md),
Sections5 and7; this map does not assert that proof's acceptance or a
concrete Lean instantiation of its geometry.

* Section5 takes one old prefix anchor P_J, J6 or7, with costs
  a_J=1/J, b_J=5/J, mu=2/5, nu=9/100. The joint capacity costs
  17/(20J), rather than two separately spendable copies of beta_J.
  Here U is the actual initial dotted foot Z5. Its residual need not
  be orthogonal to the long-end foot's residual; the written product
  `r1.r2=1-5/J` must be included.
* Section7 takes P=P5/5 and Q=alpha*P17, nu=1/16, with the native,
  finite or dotted long-end mu from Section6. Their anchor indices
  in a zero-based array are4 and16. The helper18 residual has the
  NONZERO reference product `w.r2=b0/sqrt(2)`, where b0 is the
  written first-block height. The theorem keeps the actual w.U,
  rather than replacing it by alpha*w.P17. Its mutual product with
  the head residual is also retained.

Actual farthest-point geometry, finite/parallel/dotted contact inventory,
feasibility and all-wall signs of the auxiliary points, existence of the
concave old reference and its supported floors, primitive capacities,
the actual residual norms/products and their upper bounds, rational or
radical estimates, the strict final comparison above1 and `x.x<1`
remain written application obligations. The module proves neither a B20
exclusion nor a global dimension bound. The fresh [full mathematical
source review](../results/review_two_foot_paid_kernel_before_compile_20261021.md)
binds its preserved predecessor; the delta review binds the explicit
algebraic successor. Actual bounded execution is scoped to the five
generic conditional statements, not to an instantiation of the geometry.
