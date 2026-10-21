# Independent review: Bugaenko-facet binary reduction and collar

20 October 2026. Source and paper review by `batch_guard`; no checker
import, scientific execution, root search, Lean or TeX job.

**PASS within the stated structural scope.** I read the complete source
`results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md`,
SHA256 `774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227`,
and independently rederived its algebraic and geometric bridges. This
review does not establish a new dimension bound or a classification of
all binary completions.

## Coverage and independent checks

1. **Field and integral module.** The connected intrinsic 3/5 diagram
   permits common norm normalization over the ambient field. A
   nonphysical embedding fixing `phi` would preserve the indefinite
   rank-nine facet Gram, so the actual field is `Q(phi)`. In dimension
   nine the trace of two norm-two reflections is `6+(v,w)^2`.
   Commensurability gives integral traces, hence integral pairings and an
   even orbit module. Its finite generation follows from bounded
   denominators for a finite coset decomposition. Rank nine would give
   an invariant hyperplane and an unbounded distance function on a
   compact quotient; its rank is therefore ten. E8 unimodularity splits
   the module, and evenness plus the unit `phi` makes the time vector
   primitive. The two residue classes give exactly the displayed even
   binary forms with cross entry zero or one. Here the odd glue has
   determinant ideal `(2)`: its quotient is `O/(2)`, of cardinality four
   as an underlying integral lattice. This clarifies the source's
   O-module phrase “index-two.”

2. **The infinite local family is genuine.** I checked the determinant,
   admissibility threshold `k>=-1`, Lucas norm formula and strict growth.
   The thirteen-wall seed has the required native finite angles or
   ultraparallel pairs, an actual compact base facet and a finite E8
   edge. The displayed escaping point has exactly square
   `D(D-8)/16`; its facet evaluations are the claimed positive multiple
   of those at the known base-facet vertex. For `k>=0` it gives an open
   ideal escape after perturbation. The two displayed integral basis
   changes reproduce the stated closed-module identifications; equality
   of discriminant norms is not used as an isometry test.

3. **The involution statements have the appropriate boundary.** The
   extension of the intrinsic reversal fixes the perpendicular binary
   normal and is integral because `S*u0-u0` is in `2M0`. This alone
   does not preserve a full chamber. The second point inversion is
   integral exactly when its displayed denominator is a unit: the
   relevant vector and covector have unit coordinates. Its Lucas norm
   is a unit only at exponent two in the timelike range. Neither
   involution is silently promoted to a universal period.

4. **The uniform collar is valid for every norm class.** For a primitive
   integral reflector of square `2d`, preservation of the lattice gives
   `(m,R) subset dO`. E8 unimodularity yields `m_E=d*z`, and evenness
   yields `T=d*z^2/2 in O`. At the other real embedding, positivity of
   the complementary block gives `0<T'<=1` when the E8 component is
   nonzero. Equality forces the complementary component to vanish.
   In the mixed case `T>1`, and writing `T=A+B*phi` forces integers
   `A,B>=1`; therefore `T>=phi^2`. Reverse Lorentz Cauchy then gives the
   claimed distance `asinh(sqrt(phi))` from the entire hyperbolic binary
   subspace. Pure E8 reflectors reduce to the ordinary E8 zero
   arrangement: primitivity makes their half norm a totally positive
   unit, and norm-two vectors in E8(O) are ordinary E8 roots after unit
   rescaling. The tube is in an actual full chamber, not merely a
   finite prefix. Its Fermi integral is exactly `phi^4/8`, giving the
   displayed volume lower bound.

5. **The determinant estimate uses the actual binary edge.** For odd
   glue I rederived both formulas for `q` and `q-4`. Reflection
   integrality makes the square factor in equation (5) an element of
   `O`; its nonzero norm has absolute value at least one. The conjugate
   binary Gram gives `0<q'<4`, so `q'(4-q')<=4`. Thus
   `abs(N(D))<=4*q*(q-4)` follows without a restriction on the
   intervening wall's norm. The direct case similarly gives the stated
   weaker constant 64 for the binary determinant. These are actual
   adjacent binary walls and actual axis length, not an unproved
   identification of the first wall with a convenient seed root.

## Inputs and conclusion

The metric-facet rigidity source is
`results/coxeter_metric_facet_angle_rigidity_20261020.md`,
SHA256 `f688e520f3ff7ac43225aaacd89413c03f5dd3d2d008d7e4dcdcc48326ac83fb`.
The inherited native facet, reversal and full-arrangement calibration are
those bound in
`results/compact_h9_discriminant19_d7_edge_period_20261020.md`,
SHA256 `ecc61a2dcff6df1678cc94533aa201c70cd0b262e26deab5ef85f5d64ebaf824`.
Their published arithmetic, finite-diagram and acute-face inputs remain
written inputs; this is not end-to-end geometric formal verification.

The audit validly reduces the metric-Bugaenko-facet class to two binary
glue types and proves a determinant-dependent volume lower bound with a
uniform mixed-mirror collar. It does not force such a facet in general,
prove a volume upper bound for these full chambers, preserve the base
facet under arbitrary extra norm classes, or exclude all large
determinants merely from the lower bound. Those global gaps remain open.
