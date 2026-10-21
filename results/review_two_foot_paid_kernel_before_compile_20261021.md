# TwoFootPaidKernel: independent source review before compilation

21 October 2026. Foundations reviewer. **PASS for the mathematical
statements, proof structure and assumption map; compilation and actual
axiom-output verification remain pending.** No Lean, scientific import
or other proof job ran, and no source/map was changed.

| Reviewed source | SHA256 |
|---|---|
| formal/CoxeterProofs/TwoFootPaidKernel.lean | cb282fd01130927572c447185f14c60f8f6580bad96b871a07c66fc57348d421 |
| docs/two_foot_paid_logic_20261021.md | e87c3c4213b06156c0b8c03123294e103126297f12a80dcd8dcf2effd69ed510 |

I read all196 module lines and the complete map, and inspected the
finite-sum, dot-product, retained-pressure and coordinate-square lemmas
used by the proofs. I did not author this new module. I DID author its
inherited MultiAnchorPaidFootKernel and share the actual-foot background,
including review of the target B20 proof. Thus this is independent review
of the new two-foot derivation, not an independent authorship claim for
all its foundations.

The only direct import is MultiAnchorPaidFootKernel,
`801ccfc1dd19afee9ab65451ddc63640d95efc7a4c0d802d6433dfe4ff897d4f`.
The relevant inspected dependencies are AnchorAndPaidFootKernel
`ac811bf2cd70cc975ccbbcd446d0a5d575311e5506e7141a967d150e7a639457`,
SummedActualFootKernel
`41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28`,
RightAngledProjectionHeightKernel
`2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90`
and RightAngledCubicScalarKernel
`f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee`.

## All five statements

Write w=sum(beta_i Omega_i), P=sum(a_i Omega_i), Q=sum(b_i Omega_i),
r=Z-P, s=U-Q and c_i=beta_i-mu*a_i-nu*b_i.

1. `two_payment_coefficients`39–56 correctly derives beta_i>=0 and
   c_i>=0 from nonnegative costs/multipliers and the JOINT capacity.
   Since beta_i=c_i+mu*a_i+nu*b_i, positive residual support implies
   positive original support. No residual sign or support inclusion
   is independently assumed.
2. `two_weighted_payment_decomposition`58–72 is the exact pointwise
   identity w=sum(c_i Omega_i)+mu*P+nu*Q. It follows from ring
   distributivity, finite-sum congruence, additivity and scaling. It
   charges both payments against the same original coefficient.
3. `two_shift_square`74–83 expands the full square of w+mu*r+nu*s.
   All six terms are correct, including2mu*w.r,2nu*w.s and
   2mu*nu*r.s. Symmetry only reverses dot-product arguments; it does
   not assert that any of these products vanishes.
4. `two_foot_paid_pressure`85–140 applies the original individual
   floor only where c_i>0, using the derived implication beta_i>0.
   The inherited retained-pressure lemma proves w.C<=x.C for
   C=sum(c_i Omega_i). The decomposition gives
   E=w.C+mu*w.P+nu*w.Q and the actual paid vector is
   C+mu*Z+nu*U. Multiplying each supplied foot comparison by its
   nonnegative multiplier yields exactly

       x.W >= E+mu*(h-w.P)+nu*(t-w.Q).

   There is no subtraction of an unknown original height from a known
   lower bound and no use of aggregate reference pressure in place of
   individually available floors.
5. `two_foot_paid_energy`142–188 uses the proved coordinate inequality
   x.x>=2x.W-W.W, obtained from the nonnegative sum of(x-W)^2.
   Substituting the derived pressure and full square, and using
   w.r=w.Z-w.P and w.s=w.U-w.Q, gives exactly

       x.x >= E+2mu*(h-w.Z)+2nu*(t-w.U)
                  -mu²*r.r-2mu*nu*r.s-nu²*s.s.

   The reference products cancel by algebra. Neither residual needs
   to be orthogonal to w, and the mutual product remains with its
   correct sign and factor. The proof completes ONE square on ONE x;
   it never adds two separately established bounds for x.x.

The statements are valid over the declared linearly ordered commutative
ring with all sums explicitly bounded by n,d. They require no division,
real-to-integer replacement, positive reference energy or geometric
realization to prove these conditional algebraic implications. The
proof terms use the stated finite-coordinate identities and ordered-ring
operations; I found no mathematical or obvious source/elaboration blocker.
Only a real compiler run can establish elaboration success.

## Application map and exact scope

The map accurately distinguishes this theorem from the imported one-foot
energy theorem with its flatness premise. In B20 Section5, the shared
anchor costs1/J and5/J with mu=2/5,nu=9/100 consume17/(20J), and
the nonzero product1-5/J must be retained. In Section7, the two paid
anchors are P5/5 and alpha*P17; indices4/16 are the correct zero-based
locations for those prefixes. The helper residual's nonzero w.s term
is respected through the actual w.U. The module supplies the common
algebra; it does not instantiate those geometric pressures or capacities.

The new module has five named theorems and exactly five matching
`#print axioms` commands. It introduces no axiom, sorry, admit or unsafe
proof escape. No actual print result or verified-declaration increase
is inferred from those source commands. Root must separately bind the
compiler, dependency objects, successful compilation and actual outputs.

Actual wall/foot feasibility, cap inventories, the concave reference and
its supported floors, primitive joint capacities, residual estimates,
strict numerical comparison above1 and x.x<1 remain application inputs.
There is no B20 exclusion, dimension theorem or geometric Lean claim in
this module or this source-review verdict.
