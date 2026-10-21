# Continuous-loss absorption kernel: independent source review

21 October 2026. **PASS for the stated conditional algebraic interface;
compilation and actual axiom-output verification are not performed here.**

Source: `formal/CoxeterProofs/ContinuousLossAbsorptionKernel.lean`, SHA-256
`f52af796a9cb2f3f77ee3740fc671ef0963d2130ec3c69af34b270c5007e0b22`,
7215 bytes,153 lines. Its sole direct import is the existing
`CoxeterProofs.SummedActualFootKernel`. No Lean, scientific or TeX job
ran; no source was changed. I did not author this module. Earlier shared
actual-point/reference work and the independently reviewed imported
summed-foot kernel remain explicit foundations.

All seven theorem statements and their mathematical proof steps were
checked, including the finite-coordinate integration:

* Lines31–41 derive the circle identity from `s²=v-v²` and the exact
  rotated-square identity. Neither is an assumed maximum.
* Lines44–65 obtain `2(av+s)<=a+r` from the sum of two nonnegative
  squares and `1+a²<=r²`, with r>=0 selecting the upper branch.
  The half-circle estimate follows by the analogous square comparison.
* Lines69–94 cover the ENTIRE closed interval0<=v<=1. On2v<=1,
  nonnegative a and `2s<=1` make the adjusted left side negative
  under `a+1<2ell`. On2v>1, the product
  `(2v-1)(3-2v)>=0` gives `4(2v-v²)>=3`; nonnegative b and beta>=b
  combine this with `4(a+r-2ell)<3b`. The strict conclusion follows
  without division, a field instance, optimization or calculus.
* Lines96–102 multiply the nonnegative absorbed residual by beta.
  This derives exactly the energy-correction sign used subsequently.
* Lines106–143 instantiate the imported derived summed-foot inequality
  with `T_i=w·omega_i+ell-loss_i`. Each positive-weight correction
  is nonnegative by the preceding lemma; a nonpositive alternative
  under beta>=0 is exactly beta=0, where the whole correction vanishes.
  Summing then proves `||w||²<=||x||²`. Empty finite index sets and
  zero weights require no positive-weight loss-absorption claim.

The normal assumptions are primitive unit diagonal and nonpositive
off-diagonal products, nonnegative shift parameters, and the reference
pairings `w·normal_i=beta_i`. The summed shift-square estimate and
completed-square energy conclusion are derived in the imported kernel;
neither energy comparison is supplied as a premise. No nonadjacency,
dimension, complete cap classification, field restriction or concrete
reference-energy claim occurs in the new theorem hypotheses.

One interface detail should remain explicit: lines115–117 require the
foot inequality at ALL indices, including zero weights; only the
absorption hypothesis is restricted to positive weights. This is not
a logical defect. An omitted slot can be completed with v_i=0 and
`loss_i=w·omega_i+ell-x·omega_i`, making its foot equality automatic
without a geometric foot or loss bound. Unit/obtuse/dual hypotheses
still apply to those old normals. A claim that the theorem's literal
foot premise is already positive-index-only would be inaccurate.

The continuous maximum and two-half criterion are generic internal
results. Applications still supply root/radius identities and scalar
comparisons, an actual coordinate/reference realization, primitive
normal signs and the actual available foot pressures. In particular,
the arithmetic A15–A17 inventories, field/time packing, concrete beta
minima, rational references, endpoint constructions and any dimension
consequence are not formalized by this module. Seven print commands
are present; their actual successful outputs and axiom lists await
the separately authorized bounded compilation/audit.
