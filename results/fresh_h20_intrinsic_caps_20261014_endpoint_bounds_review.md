# Review of the all-radius endpoint-cap obstruction

The [endpoint bound](fresh_h20_intrinsic_caps_20261014_endpoint_bounds.md) passes source/proof review. A [separate small arithmetic checker](../scripts/fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review.py), using the earlier independent SymPy/rational sign implementation, rechecks all 260 radius bounds, the 90 retained exceptional-right-factor cases, all eight golden-gap tightenings, the five critical half-vertex margins, and the two infinite-tail polynomial identities. Its [certificate](fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review.json) binds the frozen sources. Generation and saved replay pass.

This is not a second implementation of the 16,272 transported point/profile tests or the endpoint census's rejected-profile enumeration. Those portions receive source/proof review, with the primary saved certificate as their arithmetic evidence.

No omitted case or mathematical issue was found:

- For every omitted wall, an acute Coxeter ratio is either zero or at least `1/2`. Above the minimum ordinary bound `A/(B+1/2)`, at least one ratio must be exactly zero. Enumerating every `A/B` with `B>0` and retaining all other ratios that are zero or at least `1/2` is therefore a safe relaxation. Negative ratios are rejected; equality and multiple simultaneous right contacts are retained.
- The special E6-center improvement is valid for every profile satisfying its hypotheses. The untouched E6 ratio `4/R` forces `R<=8`. The untouched triangle ratio is then strictly above `1/2`, so the next permitted finite cosine is `1/sqrt2`, giving `R<=3phi*sqrt2`. Neither omitted factor can have a right exception because its subtractive term is zero.
- The complete profile cost satisfies `p<=1`, including the PSD boundary `p=1`. The ball power `F+2Rt+(p-1)R²` is concave. Checking both `R=0` and the largest relaxed radius therefore supplies the uniform lower bound throughout the entire interval, including isolated exceptional radii. The equivalent nesting argument uses center speed `sqrt(p)<=1`.
- The isolated dihedral tail includes every integer label and the parallel limit. Its omitted E6 right exception would require `R=12/c>=12`, which conflicts with `R<=6phi<12`. The remaining bound `12/(1+c)` is valid over the enlarged interval `1/2<=c<=1`. The half-vertex polynomial becomes `108+210t+58t²`, and the incenter excess is `56+376c+320c²`; both are strictly positive on the stated intervals.
- The profiles include finite and ideal opposite endpoints only as local PSD conditions. The proof does not silently infer that a singular endpoint has a bounded cusp link. Its conclusion is an escape obstruction, for which the broader necessary PSD class is sufficient.
- The common half-vertex point is strict on the 22 seed walls. The uniform ball-power margin permits a positive-square Lorentz lift that is strict on every allowed cap. For any finite selection a sufficiently high point is a common strict timelike interior point. Thus the conclusion is genuinely infinite volume, not just a missed point in a numerical sample.

The scope is precisely the retained intrinsic seed at inradii `(1,phi,1)` and caps that can terminate an original cusp edge. Interior caps, changed seed walls, other inradii, or different incident faces are not excluded. The successful [half-vertex layer](fresh_h20_intrinsic_caps_20261014_golden_extension.md) lies outside this class.

## Independent formula check of the proposed annulus obstruction

The parent's annulus argument is also sound with its radius assumptions explicit. The frozen vertex squared norms are `302+4phi`, `326+4phi` and `446+4phi`, so all vertex radii lie strictly between17 and22. For the first vertex, `x=3z_first/4` has radius strictly between13 and14.

Radius-one vertex balls lie wholly outside radius16; radius-one half-vertex balls lie wholly inside radius12. A quarter-face center is in one quarter of the convex hull of all vertices, so every radius-`3/2` quarter-face ball lies inside radius7, regardless of the subset recipe or compatibility selection. Finally, a cap of radius `R<=2` meeting the central radius-two cap at a finite Coxeter angle or in parallel has center distance

```text
d² = 4+R²+4R*cos(theta) <= (2+R)²,
d+R <= 2+2R <= 6.
```

Consequently every ball in these whole families misses x. This is a source/formula review of the reduction; the parent's separate annulus certificate supplies its frozen arithmetic. The central-contact statement is conditional on `R<=2`; it is not a radius bound for every imaginable central-contact cap. Larger-radius caps, other face placements and replacement remain outside this obstruction.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_endpoint_bounds_review.py --verify
```
