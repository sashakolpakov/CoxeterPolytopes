# Finite-volume tail scalar kernel: assumptions and scope

`formal/CoxeterProofs/FiniteVolumeTailKernel.lean` now has actual **bounded compilation
and nine-declaration axiom-audit PASS**. Initial compilation took1.18s;
the initial1.68s audit was followed by a separate output-capture audit,
whose exact1023 child bytes match both original execution hashes.
Only the standard axioms occur; no `sorryAx` is present.
It imports only `Std` and the pinned toolchain's `Init.Grind.Ordered`;
it changes no frozen module or manuscript. Its verification catalogue
selection is recorded separately from execution.

Written source: `manuscripts/finite_volume57.tex`, SHA-256
`c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0`.
The relevant inputs occur at lines 1837–1974, 2226–2250, 2498–2519,
2590–2648, and 3132–3188. Existing `ManuscriptFinite57.Earlier58Tail`
still assumes `n<=58`; this successor does not instantiate or alter it.

The source derives three concrete pieces of arithmetic:

* **Common payment subtraction.** A list of separate finite-vertex costs
  is summed from its pointwise cap. A separate list stores each actual
  cusp's cost and curvature and sums their pointwise comparison. A
  positive shape correction contradicts the global payment upper bound
  whenever `U<=beta`. Empty cusp lists are included. These natural-number
  quantities must already have a common positive denominator cleared.
  The displayed factor-eight caps are integral, but individual charges
  can still have larger denominators. The bounded-tail application
  therefore retains an arbitrary additional natural scale `s`, with
  finite cap `s*u8(n)` and shape coefficient `s*beta8(n)`. The geometric
  application chooses a positive `s` clearing every actual charge and
  curvature denominator; the kernel proves the scaled comparison from
  the unscaled one, rather than assuming scale eight clears all costs.
* **Every bounded scalar cap.** The literal `u8` definition includes the
  59–69 exceptional values and the 70–77, 78–116, 117–126, 127–136,
  137–146, 147–156 and 157–163 linear branches. `beta8` uses the displayed
  parity formula. `bounded_cap_table_check` proves all 105 comparisons
  over `Fin 105` with ordinary kernel `decide`, and
  `bounded_cap_comparison` transfers the result to every natural
  `59<=n<=163`. No supplied comparison or external certificate is used
  in that table theorem. Values of `u8` outside that interval have no
  claimed interpretation.
* **Narrow-cusp scalar closure.** Over any linearly ordered commutative
  ring, use `p=-(x,y)`, `H=q+q^{-1}`, `W=sum W_a`, `T=sum T_a`, and
  `N=n-1`. The inputs are `2p=MW-2H`, `H<=2M`, and
  `9T=2N-2`. The near deficit `T-W<=4` gives
  `9p>=M(N-37)`. The far deficit
  `3(T-W)<=4(C0+C1)` together with the separate caps `C0,C1<=2` gives
  `9p>=M(N-43)`. The source derives contradictions with `p<14` when
  `N>=163` and respectively `M>=1` or `20M>=21`. The final theorem
  covers both ratio ranges by an actual order case split.

The scalar kernel does **not** prove existence of a simultaneous payment
rule, face payment, correctness or coverage of the finite forest tables,
the identification of `u8` with actual finite costs, cusp cost/curvature
domination, concentration, or the positive geometric shape correction.
For the unbounded cusp branch it does not prove actual framing-wall
existence, the Lorentz geometric interpretation or common-ratio identity,
the all-factor near loss bound, or the coefficient-four-thirds local
inequalities. Those remain the explicit primitive inputs displayed in
the manuscript. The inequality `M>=1` comes from the geometric definition
`M=max(q,q^{-1})`; it is not derived by this division-free kernel.

In particular this draft does **not** establish the complete `n<=58`
tail, the compact-facet alternative, or an end-to-end `n<=57` theorem.
Its added verification content is the literal 105 scalar comparisons
and the derived payment/Lorentz arithmetic, now compiled and separately source-reviewed. The remaining geometric/table
bridges still require their own mathematical interpretation.
