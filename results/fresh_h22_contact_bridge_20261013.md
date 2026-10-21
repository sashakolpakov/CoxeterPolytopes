# An all-real restriction on bridges from the four H22 corners

No bridge exists in the following defined family. Retain the [860-wall π/6 construction](fresh_h22_multiquadratic_20261013.md), choose any one of the 22 edge supports at any of its four verified positive-side quadrilateral vertices, and require the new cap to be perpendicular to those 21 walls. A positive-height cap in this family cannot meet any of the 301 positive-height rational parent caps at a finite Coxeter angle while remaining acutely compatible with all 860 walls. Even parallel contact with one of these parents is excluded.

A separate exact extension checks the other four caps of the existing 305-cap finite-contact component. All 312 additional real contact intervals are empty. Thus **direct finite or parallel contact with any cap in that entire 305-cap component is impossible in these 88 edge families**. This does not prohibit an indirect connection through further new caps.

The result covers all real positions of each such cap and every integer angle label. It does **not** cover nonright edge profiles, other vertices, lower-height caps at these same vertices, different anchor families, or replacement of retained walls. It supplies a reason to change those constraints; it is not a general completion obstruction or a new finite-volume polytope.

## Continuous reduction

At a specified finite vertex let its22 incident normals have positive-definite Gram matrix A. Normalize the future vertex vector to `v²=−1`, and let `d_j` be the dual of the released wall in the incident span. Put `e=d_j/sqrt((A⁻¹)jj)`. The exact certificate reconstructs these vectors and checks

`e²=1`, `v·e=0`,

with v and e both orthogonal to the21 retained edge walls. The inverse A⁻¹ is rational and entrywise nonnegative.

Every norm-two normal perpendicular to the edge lies in this Lorentzian plane. Acuteness with the released wall forces its e coefficient to be negative. Its v coefficient must be positive for any completion with a strict timelike interior, without requiring the original center to remain inside. Indeed, every future strict point in the old vertex sector decomposes as

`x=αv−Σγ_i d_i`, with `α>0` and `γ_i>0`.

Because `(e,d_i)≥0`, a cap `r=a v−sqrt(a²+2)e` with `a≤0` has `(r,x)≥0`. It cannot contain such a strict point. Thus the complete admissible branch has the parameterization

`r=(t(v−e)−t⁻¹(v+e))/sqrt(2)`, `T=t²>1`.

For an old wall w of norm n write `a_w=(v−e,w)` and `b_w=(v+e,w)`. Its acuteness condition is the linear inequality

`a_w T−b_w≤0`.

The positive-height requirement supplies one more inequality of the same form. For a proposed finite or parallel anchor contact, the correctly normalized cosine and necessary quadratic inequality are

`cos θ=(−a_w sqrt(T)+b_w/sqrt(T))/(2sqrt(n))`,

`(a_w T−b_w)²≤4nT`.

The coefficient4 here includes the norm 2 of the new cap. The proof does not use a finite list of anchor angle labels.

## Exact interval certificate

The [verifier](../scripts/fresh_h22_contact_bridge_20261013_verify.py) reconstructs all 860 normal recipes using rational arithmetic and sparse squarefree radicals. It checks their norms, reconstructs each incident rational Gram matrix, proves positive definiteness and nonnegative inverse, and verifies every orthogonality and normalization in the parameterization. The actual vertices and the old Coxeter pair checks are inherited from the frozen, source-bound 860 certificate and its [independent review](fresh_h22_multiquadratic_20261013_review.md).

For each family a numerical scout proposes a small set of useful inequalities. The verifier computes their exact radical coefficients and outward-rounded rational bounds; it never trusts the scout's floating bounds or classifications. Its [saved certificate](fresh_h22_contact_bridge_20261013_verify.json) records:

| Necessary gate | Exact count |
| --- | ---: |
| Edge families considered | 88 |
| Families excluded by incompatible acute/height bounds | 10 |
| Remaining families | 78 |
| Parent-anchor tests on these families | 23,478 |
| Empty finite/parallel anchor intervals | 23,003 |
| Remaining intervals with a uniformly forbidden wall angle | 475 |

In the 475 remaining cases, another retained wall has normalized cosine strictly between 0 and 1/2 throughout the entire allowable anchor interval. No angle π/m for an integer `m≥2` lies there. The range is certified by the sign of the derivative and rigorous endpoint bounds. This excludes every finite anchor label simultaneously.

All interval exclusions are strict. Closed acuteness and parallel endpoints, including possible singleton intervals, are preserved in the overestimate unless a strict inequality proves them impossible. Positive height is relaxed to nonnegative height in the necessary bounds. Thus boundary cases are not discarded by numerical tolerance.

The initial bounded scout used the same 88 edge supports, the 301 parent anchors and prescribed labels 2,3,4,5,6,8. It found no compatible candidate among 314,244 numerical quadratic branches. That numerical failure motivated the continuous proof but is not its justification. The certificate imports no numerical scout module, optimizer or floating geometric calculation.

The [four-anchor extension](fresh_h22_contact_bridge_20261013_big_verify.json) preserves the preceding certificate unchanged. The [certified contact graph](fresh_h22_covering_cells_20261013.md) identifies the extra types as 0,2,6,22, with positive-height physical indices 812,816,824,854. Its [separate verifier](../scripts/fresh_h22_contact_bridge_20261013_big_verify.py) reconstructs the 88 frames and necessary bounds again, checks those source/index mappings, and proves that all `78·4=312` further finite/parallel upper bounds lie strictly below the necessary lower bound. No additional angle classification is needed for these four types.

## Consequence and replay

These all-right edge insertions cannot supply a direct finite-contact bridge from the separated new-cap component to the existing 305-cap component. A constructive continuation must change at least one of the theorem's hypotheses, for example allowing nonright contacts on the edge, choosing another face or replacing retained walls. An indirect chain of additional caps also remains outside this theorem. Adding distant caps dotted to the entire component does not itself provide the requested bridge.

```sh
python scripts/fresh_h22_contact_bridge_20261013_verify.py --verify
python scripts/fresh_h22_contact_bridge_20261013_big_verify.py --verify
```

The corrected exact producer and saved-result replay pass. The 860 construction and all its earlier files remain unchanged. No new normal or finite-volume claim is made in this checkpoint.

The [independent source/theorem review](fresh_h22_contact_bridge_20261013_review.md) confirms the normalization, the strict-interior argument for `t>1`, closed endpoint handling, full anchor coverage and derivative bounds, and independently invokes the frozen parent-bank replay. It uses the same exact implementation for that replay rather than claiming a second arithmetic reconstruction.

The separate [four-anchor review](fresh_h22_contact_bridge_20261013_big_review.md) also passes, including a fresh saved-result replay, the exact type/index/norm mapping and all312 strict additional interval exclusions.

The separate [four-anchor review and fresh replay](fresh_h22_contact_bridge_20261013_big_review.md) also pass. Together these establish the stated direct-contact obstruction for the entire 305-cap component within the 88 specified families.
