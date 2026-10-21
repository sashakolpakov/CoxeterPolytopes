# Anchor and paid-foot kernel: source review

21 October2026. **Mathematical/source-logic PASS, uncompiled in this
review.** No Lean job, import execution, old-kernel rerun or resource
admission was attempted. Elaboration and actual axiom output remain
separate required checks; this verdict does not predict their outcome.

Reviewed source `formal/CoxeterProofs/AnchorAndPaidFootKernel.lean`,
7,260 bytes/159 lines, SHA256
`b62bf37921cb0a2ff42f223f3e5287d3465402a28c2ee258396f088ebac32d99`.
Its direct imported source `SummedActualFootKernel.lean` is unchanged,
SHA256 `41d183ca03e0e405ce3d86ee838bdc3859c03f23ca19bb2cf9cd5239da050c28`.
The import's finite-dot linearity and coordinate completed-square lemmas
were inspected at the exact interfaces used. No new axioms or sorry
appear in the new module. It has six meaningful theorem axiom-print
commands, but their ACTUAL output has not been obtained in this review.

The reviewer authored earlier A18 supporting sources and the original
summed-foot module, with shared root inputs. Those contributions are
disclosed. The new module is root-authored; its full new logical chain
and assumption boundary were inspected here.

`same_k_reflection` is the correct concrete statement on Fin19, including
the harmless index0: equal i(19-i) gives equality or index sum19. It uses
`decide`, rather than a native-computation axiom. Until compiled it is a
source-level intended kernel computation, not actual verified evidence.

`weak_index_unique` derives uniqueness from necessary parameter/k matches,
an actual weak anchor, the reflection identity and a supplied reflected
obstruction. No equality of parameters, census emptiness or actual wall
existence is concluded from a formal model. In an A18 application `weak`
must include the central5..14 restriction (false outside it); the script
does not cover an unrestricted height predicate on all Fin19. Extraction
of actual cap parameters, the complete match implication and arithmetic
Cartan obstruction remain external hypotheses. This application boundary
is essential and should appear in the documentation map.

`remaining_weight_pressure` handles precisely the finite retained sum.
Nonnegative weights and floor comparisons are used only where the weight
is positive; zero weights contribute zero. Thus it does not demand the
unknown omitted height. Dot linearity and finite-sum monotonicity derive
the weighted inequality instead of assuming it.

For `paid_reference_pressure`, the supplied payment is
w=B+mu*a*P with B the retained nonnegative sum. The theorem FIRST derives
E=w.B+mu*a*H from the actual reference norm and height identities. It then
combines the retained-floor inequality with mu>=0 and x.V>=g to obtain
x.(B+mu V)>=E+mu(g-aH). This is the correct paid subtraction; it never
subtracts an unknown height from a known lower bound. No separate sign
on a is necessary for this more general algebraic theorem.

`paid_reference_energy` derives the SAME-vector equality
B+mu V=w+mu(V-aP). Its full square is expanded from coordinate dot
linearity. The supplied reference orthogonality w.(V-aP)=0 and actual
residual square S yield E+mu^2*S; that desired full graft-square identity
is not an input. The imported coordinate nonnegative square then derives
the claimed lower bound E+2mu(g-aH)-mu^2*S on x.x. The final contradiction
uses a strict scalar lower bound above1 and the actual x.x<1. Weak <=
floor/foot pressure is sufficient because that final scalar inequality
is strict. No independent component energy is added.

Concrete coordinates, paid coefficient realization, original reference
norm/height, residual orthogonality and square, actual-foot feasibility,
scalar estimates, field reduction and full cap classification/census are
ALL still application inputs. Their presence is explicit and appropriate
for this internal logic kernel. The module proves no full A18, quadratic
field, polytope or dimension theorem by itself. Source correctness found
no mathematical blocker; compilation must remain deferred until the
unchanged6GiB/40GiB admission gate is met.
