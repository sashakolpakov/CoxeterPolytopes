# Nonright S61 contacts: an exact finite-bank obstruction

The retained 133-wall H20 chamber admits many numerical nonright-S61 escape cutters, but the complete saved **847-recipe bank still has a common exact strict spacelike escape**. Every subset of these recipes therefore leaves infinite volume. No 135-wall construction is promoted merely for its compatible right-angle pair: the three such numerical pairs create no finite corner in the saved selection screen.

There is also an all-real contact constraint. For every nonzero coherent acute added wall, the component perpendicular to S61 must be **future timelike**. Thus its support projection has square strictly greater than its ambient square. This constraint covers arbitrary real normals and all contact labels; it does not exclude the remaining timelike-residual families.

The prior [133-wall compact-face construction](fresh_h20_s61_multicap_20261014.md) and its frozen sources are unchanged. The present result does not establish nonexistence of a finite-volume H20 completion.

## Exact certificate and scope

The [checker](../scripts/fresh_h20_nonright_s61_20261014_bank_gate.py) and [certificate](fresh_h20_nonright_s61_20261014_bank_gate.json) reconstruct every archived proposal from:

1. its named one- or two-node S61 contact profile;
2. its two named anchor contacts, with labels 2–6;
3. the norm-two equation;
4. a rational box isolating one of the two quadratic branches.

Only anchors 0, 65, 84, 92 and 95 occur. Their physical pairings are reconstructed using the frozen radical-tower data. Named cosines are exact radicals. The rounded scout coordinates only choose a branch: the checker proves that exactly one algebraic root lies in the displayed rational box, and the other lies outside it. Nonsingular linear systems, nonzero quadratic coefficients and strictly positive discriminants are certified by outward rational intervals.

The certificate proves all 847 roots have square two and future-timelike S61 residuals. It checks their strict signs at a common timelike interior point and at the new rational escape. It also verifies that every one cuts the previous escape. **It does not certify their remaining numerical Coxeter-angle fits or mutual compatibility.** Those extra conditions are unnecessary for the common-escape conclusion: even admitting every recipe as a halfspace leaves the escape.

All 133 retained-wall escape signs are separately re-evaluated, including the unexpanded six-cap hexagon and 26-cap algebraic chain. The new rational point has exact square

`59617546381486495865787615871 / 10000000000000000000000000000`.

The largest retained normalized-wall pairing is below `−6.03975×10⁻⁷`; the largest bank pairing is below `−1.5677961×10⁻⁶`. These decimals summarize rigorous rational bounds. All 847 common-center pairings are below `−0.8643`. The point has negative pairing with that future timelike center.

The segment from the common timelike center toward the strict spacelike point crosses the future light cone while every inequality remains strict. An open ideal-boundary patch is consequently retained, proving infinite volume for the full bank intersection and every subset. This is a **finite archived recipe-bank theorem**, not an exclusion of all continuous two-anchor solutions, all S61 contact profiles or all possible added walls.

## Why every valid residual is timelike

Let the 18 retained S61 roots be `s_i`, with positive-definite Gram `G` and entrywise nonnegative inverse `J`. Let `ω_i` be their dual weights. For a proposed acute wall `r`, write

`c_i = −(r,s_i) ≥ 0`, `w = Σ c_i ω_i`, and `r = −w+y`, with `y ⟂ S61`.

Since `w = Σ (Jc)_i s_i` has nonnegative coefficients, every nonsupport old wall `n` satisfies `(w,n)≤0`. In particular this holds for all 29 sides of the compact S61 polygon. Acuteness gives

`(y,n) = (r,n)+(w,n) ≤ 0`.

Thus `y` lies in the cone over that compact polygon. Its time orientation follows explicitly from the positive normal relation

`n119+n120+√2√(a²+1)n95 = (2a,0,0)`.

The three spatial normals in this relation span the plane. A vector satisfying all side inequalities has nonnegative time coordinate; time zero forces the vector to vanish. At positive time, dividing by that coordinate puts its spatial part in the compact 29-gon, all of whose 29 vertices lie strictly inside the unit disk. Every nonzero such vector is therefore future timelike.

The alternative `y=0` is incompatible with a common strict interior point `h`: then `r=−Σ(Jc)_i s_i`, with at least one positive coefficient, so `(r,h)>0`. Hence every coherent nonzero acute addition has `y²<0`, and

`w² = r²−y² > r²`.

For a norm-two cap, the necessary budget is `cᵀJc>2`. This also means that **no new coherent acute wall perpendicular to all S61 roots can be added to this retained compact-face chamber**. The argument uses the compact polygon and acuteness; it assumes neither arithmeticity nor an intrinsic Coxeter structure on the polygon.

The exact budget census of the 225 profiles used in the scout is:

| Support-projection square | Profile count | Consequence |
|---|---:|---|
| below 2 | 59 | impossible for a coherent acute addition |
| equal to 2 | 5 | also impossible; a null residual cannot survive the compact-face cone |
| above 2 | 161 | not excluded by this necessary constraint |

The numerical scout originally omitted the 59 positive-residual profiles and allowed null profiles in its convex relaxation. The new geometric proof closes both omitted residual types. It does not make the remaining numerical relaxation an exact global optimizer.

## Bounded constructive search and selection

The [scout](../scripts/fresh_h20_nonright_s61_20261014_scout.py) and [saved bank](fresh_h20_nonright_s61_20261014_scout.json) test 72 single-node profiles with labels 3, 4, 5 and 6, plus all 153 two-node profiles with both labels 3. Physical inputs are evaluated at 100 digits before the numerical optimization. For each positive relaxation, only its eight nearest nonsupport anchors and anchor labels 2–6 are used. The run has 150 positive numerical relaxations, 21,908 quadratic branches and 847 saved proposals.

The separate [selection diagnostic](fresh_h20_nonright_s61_20261014_selection.md) finds three robust finite-contact pairs, all right: 103/811, 104/812 and 105/813. Each uses support contacts on walls `{1,7}` versus `{9,10}`, all at π/3. Their support costs are `46/11`, their cross product is `24/11`, and the common residual square is `−24/11`, explaining the right-angle cancellation if the full anchor recipes are accepted. The shared anchor labels are respectively `(5,6)`, `(6,5)` and `(6,6)` at walls 65/95.

None of these pairs has a compatible third member in the numerical bank graph. The pairwise finite-corner screens find no positive-definite rank-20 support containing both proposed roots. The only three-member compatible clique is all dotted. These graph and vertex-screen claims remain numerical; the exact common-escape theorem does not depend on them.

The practical conclusion is to change the contact/anchor family or retained core before doing more subset selection in this bank. More general multi-node profiles and different anchor equations remain open. No larger cap bank or further construction was started for this checkpoint.

## Verification

Exact generation and saved replay pass. The original numerical scout's deterministic saved replay also passes. The checker uses inherited exact field arithmetic and a separate outward interval reconstruction of all 847 quadratic recipes; this is not an independent reimplementation of the frozen 133-wall chamber.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_nonright_s61_20261014_bank_gate.py --verify
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/fresh_h20_nonright_s61_20261014_scout.py --verify
```
