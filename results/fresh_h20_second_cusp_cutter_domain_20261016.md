# Every residual E5 cutter is dotted to both walls 141 and 144

**PASS.** The [exact verifier](../scripts/fresh_h20_second_cusp_cutter_domain_20261016.py)
and [certificate](fresh_h20_second_cusp_cutter_domain_20261016.json)
prove the following all-real necessary condition for the 151-wall chamber:
every spacelike normal perpendicular to the retained sixteen-root support,
acute to the existing walls, and cutting E5 is strictly dotted to both
wall 141 and wall 144. Generation and fresh saved replay pass. A separate
narrow source review checked the normalization, outer relaxation,
exhaustion and geometric implication; it did not repeat the large physical
geometry reconstruction.

Let C be the verified finite corner perpendicular to caps 141 through 144,
let d=-C^2>0, and let q_i=n_i^2. Acuteness to the four orthogonal caps gives

    r = alpha*C - sum(u_i*n_i),  u_i >= 0.

The strict negative pairings of old wall 133 with C and all four caps
force alpha>0 for any nonzero candidate. After positive scaling take
alpha=1 and put z_i=u_i*sqrt(q_i/d). These are nonnegative coordinates.

Only five further wall inequalities are needed: walls 31,95,140,145,150.
Each has the form A_j.z<=1. The last two rows follow directly from the
exact null-shift recipes:

    A145 = (1/z_parameter, 0, 0, 1),
    A150 = (0, 0, sqrt(q143)/6, 1).

In particular, the second row contains no additional sqrt(d) factor.
The condition r.E5>=0 is B.z>=1. Its coefficients use the independently
reconstructed E5 pairings; the E5 components along the sixteen support
duals are perpendicular to this entire residual space.

The verifier rounds each A coefficient downward and each B coefficient
upward to a rational with denominator 1,000,000. Because z_i>=0, these
changes enlarge the admissible region. Together with nonnegativity this
gives ten rational inequalities in four variables. The strictly positive
row for wall 140 bounds every coordinate, so the outer region is bounded.

All 210 four-row bases are checked by exact Fraction elimination. Of these,
195 are nonsingular and eleven yield feasible vertices. The certificate
lists every vertex and active set. Their maximal projected squared lengths
are approximately 0.984113603 after omitting coordinate 141 and 0.447314009
after omitting coordinate 144. Exact rational comparisons give

    z142^2 + z143^2 + z144^2 < 99/100,
    z141^2 + z142^2 + z143^2 < 9/20.

These convex quadratic functions attain their maxima at vertices of the
bounded outer polytope. The inequalities therefore hold throughout it.
The corresponding orthogonal projections of r have squares bounded above
by -d/100 and -11d/20, respectively. Both are timelike. For a spacelike r
acute to the two walls, this is precisely the strict dotted-contact
condition.

The result has no field, label, norm, height or retained-center cutoff.
The E5-cut and sixteen-root orthogonality hypotheses are essential. It
does not rule out cutters dotted to both walls, walls leaving this residual
space, or finite-volume completion. The numerical discovery polytope is
not a proof input.

```sh
python3 scripts/fresh_h20_second_cusp_cutter_domain_20261016.py --verify
```
