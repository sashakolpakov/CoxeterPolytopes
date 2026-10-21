# TwoFootPaidKernel: explicit pressure repair source review

21 October 2026. **PASS for the exact proof repair; current elaboration and actual axiom verification remain separate obligations.** No Lean command, worker import, proof job or source edit was performed. The reviewer authored the original preparatory module, but did not author the parent's repair. This note independently checks that repair and its unchanged mathematical interface; it is not an independent authorship review of the entire original module.

Current source: `formal/CoxeterProofs/TwoFootPaidKernel.lean`, SHA-256 `ff40c49823d399777d3f93ae869d1ee3beaa33aca616d85487b2fd27ba5f1e65`, 10698 bytes / 209 lines.

## Preserved failure and exact delta

Read `results/two_foot_initial_failure_preservation_20261021.json` and independently checked both gzip files and decompressed bytes:

| Preserved item | Gzip SHA-256 | Decompressed SHA-256 / bytes |
|---|---|---|
| `results/source_archive/TwoFootPaidKernel_before_explicit_pressure_doubling_20261021.lean.gz` | `60ba9d2818cd98ad96b9d2f67cf679dee8ae46f7fa54b5e400b6b4a9c377fc08` | `cb282fd01130927572c447185f14c60f8f6580bad96b871a07c66fc57348d421` / 9997 |
| `results/source_archive/TwoFootPaidKernel_initial_compile_failure_20261021.txt.gz` | `42d148a525db72343ef78f0e10b197d0fd8540a97bca62e70979b1bf1a4091f8` | `d3d8c36392f6577f771e851e7da180c59d9b0c6fd0223ec0825259b3863677e5` / 13750 |

A full byte-derived unified diff contains one hunk only: the final bare `grind` in `two_foot_paid_energy` is replaced by current lines 188–201. Definitions, hypotheses, all five theorem statements, all four earlier proofs, imports and five print commands are unchanged. Direct extraction also confirms equality of all five complete theorem headers. The failed initial compile and its failed energy print containing `sorryAx` earn no successful module or declaration coverage; the repair does not erase that record.

## Independent check of the repaired inference

Write (w=\sum\beta_i\omega_i), (A=\sum a_i\omega_i), (B=\sum b_i\omega_i), (r=Z-A), (s=U-B), (W=w+\mu r+\nu s), (E=w\cdot w). The already derived pressure is

\[
P:=E+\mu(h-w\cdot A)+\nu(t-w\cdot B)\le x\cdot W.
\]

It comes from the primitive nonnegative capacities, supported old floors on the **same** (x), and the two individual foot pressures. No combined-pressure conclusion is supplied as a hypothesis.

The existing full square and the two residual pairings give

\[
W^2=E+2\mu(w\cdot Z-w\cdot A)+2\nu(w\cdot U-w\cdot B)
+\mu^2r^2+2\mu\nu(r\cdot s)+\nu^2s^2.
\]

Consequently the theorem's entire left side is exactly (2P-W^2): both anchor-reference terms cancel, with the mutual residual product retained. This is the equality at lines 191–197, using only the prior square and linear dot identities. It needs neither (w\cdot r=0), (w\cdot s=0), nor (r\cdot s=0).

Line 188 explicitly multiplies (P\le x\cdot W) by the nonnegative scalar two. Subtracting the identical (W^2) therefore gives

\[
2P-W^2\le2x\cdot W-W^2\le x^2.
\]

The last step is precisely imported `SummedActualFootKernel.coordinate_energy_lower`, derived there from (0\le\lVert x-W\rVert^2). All signs and factors of two in the new equality and both inequalities are correct. The repair replaces an overloaded final automation step with the explicit valid monotonicity and transitivity steps; it does not weaken the theorem or add a conclusion-shaped assumption.

**Boundary:** mathematical source PASS only, pending the parent's separately admitted bounded retry and actual five-declaration audit. There is no B20, geometric realization, cardinality or dimension conclusion in this module. The six papers and current selected-kernel/declaration counts are unaffected by this review.
