# The old7 boundary bridge fails at original walls

There is no norm2 cap root K in the original H9 ambient space compatible
with all original walls and having both finite duties

\[
W+7+30+K,\qquad W+30+52+K,
\quad W=\{0,1,2,3,5,6\}.
\]

This is an all-finite-label statement for these two duties. It covers all
norm2 roots, not just the seven profiles used by the numerical searches.
It does not require retaining the old common center, and the three later
walls75,76,77 play no part in the obstruction. Every possible root has a
strictly positive product with original wall4 or47.

The [producer](../scripts/fresh_exact66_old7_boundary_cap_obstruction_20261018.py),
[certificate](fresh_exact66_old7_boundary_cap_obstruction_20261018.json), and
[saved checker](../scripts/check_exact66_old7_boundary_cap_saved_20261018.py)
use the frozen [original physical scaffold](fresh_remote60_exact_original_extension_20261018.md).
They bind the [exact66 seed](fresh_exact65_escape33_cap_20261018.md), whose
66 genuine facets and205 finite vertices are unchanged. No wall is deleted,
no finite-support certificate is debited, and no new chamber or cone is
asserted. Other incidence patterns remain open.

## Why only four roots are possible

The original physical products certify two disjoint H4 factors in W+7+30:

\[
0\mathrel{-5-}1\mathrel{-3-}2\mathrel{-3-}3,
\qquad
5\mathrel{-3-}6\mathrel{-3-}7\mathrel{-5-}30.
\]

In a positive-definite Coxeter diagram, an H4 component has no connected
finite extension. Consequently the first finite duty forces K to be
orthogonal to every one of these eight roots. In particular its white
profile is zero and its contacts with7 and30 are right. This uses the
complete finite Coxeter classification, not a finite label cutoff or a
restriction to the previously sampled profiles.

The original30/52 bond has label4. The three-root principal30/52/K in the
second finite duty has determinant

\[
4-2\theta_m^2,\qquad \theta_m=2\cos(\pi/m).
\]

It must be strictly positive, so m is2 or3. Label4 already gives determinant
zero and all larger finite labels give a nonpositive determinant. The two
remaining possibilities give a disconnected A1+I2(4), or B3, respectively.
Thus the complete remaining contact vectors on frame(7,30,52) are
(0,0,0) and(0,0,−1).

Project orthogonally off the six white roots. Their span is positive
definite, and the remaining physical space has signature(3,1). The
three projected frame roots have signature(2,1); their orthogonal
complement is a positive line. Let G be their Gram matrix, h the projected
old center, f its three frame products, and

\[
e=h-fG^{-1}\,\text{frame},\qquad
q=e^2=h^2-fG^{-1}f^T>0.
\]

For either required contact row b, put
\(\rho=2-bG^{-1}b^T\). Its certified interval is strictly positive.
Every solution of the linear contacts and norm2 equation is exactly

\[
K_\pm=bG^{-1}\,\text{frame}\ \pm\sqrt{\rho/q}\,e.
\]

Both signs use the same explicitly defined residual e. This proves root
completeness without relying on the numerical names of sheets or an
approximate norm solver. The formula concerns the fixed ten-dimensional
Lorentz space of H9; it makes no assertion about additional ambient
orthogonal dimensions.

## Original-wall exclusions

The producer reconstructs the physical roots and products using directed
rational intervals. Exact original Coxeter identities are isolated only
after strict interval comparisons with neighboring possible labels. The
original-native premise is inherited from the frozen scaffold; new-root
products are never rounded to native labels.

For each original blocker j, it evaluates

\[
K_\pm\cdot r_j=
bG^{-1}c_j\ \pm\sqrt{\rho/q}
\bigl(h\cdot r_j-fG^{-1}c_j\bigr),
\]

and also evaluates the same product directly from the physical vector.
The two enclosures may be intersected because they enclose the same
algebraically defined product. All four roots have a strict positive
original contact:

|Label K/52|Residual sign|Original blocker|Certified lower bound|
|---|---:|---:|---:|
|2|−|47|4/5|
|2|+|4|13|
|3|−|47|4/5|
|3|+|4|13|

The displayed bounds are strict and were checked as rational inequalities
against the saved interval endpoints. They imply failure of acuteness
before any test involving walls75–77, common-center orientation or a full
cone. The old numerical observation that the center-facing sheets fail
wall47 is therefore upgraded to an exact obstruction on both sheets.

This rejects the proposed first cap for a changed boundary7–30–52–47–7.
It does not reject all bridges around that boundary, nor any family that
changes one of the two stated finite duties or deletes an original wall.
It is a structural guide for the next incidence choice, not a completion
of the finite-volume construction goal.

## Verification

```sh
python3 scripts/fresh_exact66_old7_boundary_cap_obstruction_20261018.py --verify
python3 scripts/check_exact66_old7_boundary_cap_saved_20261018.py
```

The saved checker imports the frozen stand-alone quadrilateral-cap interval
checker for directed rational arithmetic, Gaussian elimination and adjugates;
it imports no producer modules. It independently reconstructs the white
projection, old native identity isolation, frame signature, positive line,
both complete two-root equations and original positive blockers. The
original scaffold, its native premise and finite Coxeter classification
remain inherited. No cone, historical vertex census or original profile
theorem is rerun.

Generation74626/75022b, separate producer replay23397/df6f24 and the final
saved checker e058c6 terminated exit0/PASS. The strict displayed rational
bounds passed5d87a6. Metadata25f271 passed7 bindings,5 links and2 ASTs.
The root agent completed a full mathematical read of the producer, report
and stand-alone saved checker. The review covered H4 isolation and all-root
profile zero, the all-label restriction at the fixed 30/52 label4 bond,
the physical quotient and positive-line two-sheet completeness, and both
original-wall positive blockers. It inherited the original native scaffold
premise and did not execute the producer or checker, reconstruct a full
cone, or replay the historical finite-vertex census. No correction was
requested. This packet is final and frozen after that review.
