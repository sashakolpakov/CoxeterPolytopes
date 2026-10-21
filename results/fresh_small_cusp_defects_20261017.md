# Small-cusp defects and a one-cusp consequence

This package proves a local refinement of the fine-payment theorem and an exact cusp-incidence balance. A deficient triangle or quadrilateral has a small affine factor in its original wall star. After paying small ideal corners sufficiently, every remaining carrier cusp lies **outside** the defective face. Separately, dimension is at most155 if every cusp has some irreducible factor of rank at most8; in particular, every one-cusp chamber has dimension at most155.

This package does not improve the unrestricted156 ceiling. Its155 consequence is retained as a historical refinement of that bound; the separate selective-distance13 proof, once finalized, supersedes it with a stronger unrestricted result. The local nonincident-defect statement remains a distinct ingredient. No bound on carrier multiplicity, new example, or novelty claim is made.

## Localizing the fine-payment failure

Let P be a finite-sided finite-volume hyperbolic Coxeter n-polytope. For an actual triangular or quadrilateral two-face F let H(F) consist of its n−2 supporting walls and its side walls. This is the original wall star, without artificially adding all walls at nonsimple ideal vertices.

**If H(F) contains no connected affine subdiagram of order at most9, then F receives fine payment at least2 for a triangle or at least1 for a quadrilateral.** Use the finite Vinberg weights and the distance-three/end-segment rule at ideal corners.

This is the face-local form of the frozen [rank9 theorem](fresh_cusp_minimum_factor_rank9_dimension29_20261017.md). The following audits the hypothesis change. The triangle's minimal indefinite core and its proper affine subdiagrams belong to H(F), so its Lanner reduction and finite-corner alternatives are unchanged. At an ideal corner the varying affine factor belongs to H(F), and thus has order at least10. The [ideal-corner proof](fresh_affine_order10_quad_payment_20261017.md) uses only this factor, the other two sides and the supporting spectators. Every excluded low-contact affine, including the new dual-E8 witnesses, is a principal subdiagram of H(F).

For four finite corners, both complements in the private/shared decomposition lie in H(F). Every witness in the [carrier proof](fresh_affine_carriers_order10_20261017.md) uses only its carrier, the actual white connecting path and the four sides. The later published cases use the order-at-most9 tests in the [refined source audit](fresh_vinberg_finite_quad_affine_cutoff9_20261017.md); their large orthogonal pieces are positive definite by the actual Lorentz realization. No affine prohibition elsewhere is used. The published finite case tables remain inherited inputs, not a new classification by this checker.

A connected affine diagram in H(F) determines an actual cusp and is one of its full irreducible factors. Its positive kernel combination is a nonzero feasible null vector: nonzero follows by pairing with a strict center, and feasibility follows from acuteness. Strict convexity of the Klein ball makes its ideal point an actual vertex. The positive-semidefinite incident Gram then forces every other incident root right to this affine factor. This is the [actual-factor lemma](fresh_cusp_factor_payment_framework_20261017.md), consistent with the primary unique-ideal-vertex statement in [Felikson–Tumarkin, §1, PDF p.4](https://arxiv.org/pdf/math/0502413), reopened for this package.

## A defect demand that charges only nonincident carriers

Use this nonnegative rule:

* At finite corners keep Vinberg's weights, whose vertex total is at most3n.
* At ideal corners whose varying factor has rank r<=8, use weight2.
* At other ideal corners retain the fine distance-three/end-segment rule.

Let s_large(v) count the factors of rank at least9, V_v=product_i(r_vi+1), and

` b(v)=3 s_large(v)+sum_(r_vi<=8) r_vi. `

The total cusp weight is at most V_v b(v). Each marked pair in factor i occurs in exactly V_v/(r_vi+1) link fibers, so a small factor costs
`2 binom(r+1,2)V_v/(r+1)=rV_v`.
A large factor costs at most3V_v by the [all-rank affine cost proof](fresh_cusp_spectator_role_payments_20261017.md).

Call a triangle or quadrilateral deficient when its weight is below2 or1 respectively. Such a face has no small varying-factor corner, since that corner alone pays2. The local theorem nevertheless forces a small affine factor C inside H(F).

**The cusp determined by C is not on F.** Otherwise all its walls contain some ideal vertex v of F. The walls of H(F) containing v are exactly its supporting whites and the two incident sides. The factor-fiber description of this corner star has one complete affine varying factor and one-node-deleted finite pieces from every other factor. Thus C must be the complete varying factor, contradicting the absence of a small varying-factor corner. This allows nonsimple ideal vertices.

Writing T_bad,Q_bad for the deficient counts, the weight upper bound and [exact link-shape curvature](fresh_cusp_link_shape_curvature_20261017.md) imply

` 2T_bad+Q_bad >= (beta_n−3n)v_f
                  +sum_v [kappa(L_v)−V_v b(v)]+5(n−1). `

Here beta_n=n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n, and v_f counts finite vertices. A deficient triangle is short by at most2 and a deficient quadrilateral by at most1; polygons with at least five sides need no positive payment. This proves the inequality directly, without assuming all deficient payments are equal to their maximum shortfall.

A positive right side forces face-star charges to **nonincident** small cusps. No multiplicity bound for these charges follows from the argument. In particular, counting only corners of faces actually passing through the carrier cusp does not finish a dimension proof.

## Exact example distinguishing carriers from incident cusps

In upper half-space H3 consider

` −1<=X<=1, −1<=Y<=1, X²+Y²+Z²>=4, Z>0. `

This Coxeter square pyramid has four finite base vertices (±1,±1,sqrt2) and one ideal apex at infinity. Adjacent vertical walls are right, opposite ones parallel, and the base hemisphere has angle pi/3 to each. In metric diag(−1,1,1,1), outward normals of square one are

```
R=(1, 1, 0,1), L=(1,−1, 0,1),
U=(1, 0, 1,1), D=(1, 0,−1,1),
B=(−3/4,0,0,−5/4).
```

The strict center (5,0,0,4) has square−9. The full cone has exactly the rays (5,±2,±2,3) and (1,0,0,1). The checker exhausts all ten triple frames in the section defined by the timelike sum of the normals. That section is bounded: its recession direction would pair nonpositively with every normal and zero with their sum, hence zero with all five; they span the Lorentz space. Every facet has an explicit relative-interior point from a positive sum of its vertex rays. There are no exterior rays.

The relations R+L=U+D=2(1,0,0,1) certify the two affine A1 factors and the proper rectangular cusp link. The bounded link and complete finite/ideal cone give finite volume.

The compact base's full star contains both affine factors although their cusp is not on that base. This disproves identifying a carrier in a face star with a cusp on the face. **It is not an unpaid-face example:** each base corner is A3 marked at its two end roots, and the published rank3 rule pays that pair weight2, total8. No adjacency-only finite rule is being substituted for Vinberg's actual weights.

## Exact cusp-type balance at dimension156

A separate consequence uses the [unrestricted distance14 theorem](fresh_distance14_unrestricted_dimension156_20261017.md), rather than the unbounded nonincident charge. Its cusp cost per link vertex is

` c(r)=sum_i min(14,r_i/2). `

At n=156 the total link rank is155. The frozen convex-order transfer (a,b)->(a−1,b+1), applied when both ranks are strictly below28 or both strictly above28, preserves c and cannot increase normalized curvature. Thus a reduced profile has copies of28, at most one positive rank below28 and at most one above28. All128 reduced shapes of total155 give

` kappa(L)/V−c(r) >= −4172075/82044596. `

If a cusp has a rank at most8, hold one such rank a fixed and transfer only the other factors. For a=1,...,8 this gives988 marked reduced shapes, proving

` kappa(L)/V−c(r) >= 1395158/20511149 >0. `

This does not assume that unrestricted transfers preserve a small factor. The exact minima occur at (15,28,28,28,28,28) and (1,14,28,28,28,28,28). Standalone rational convolution and independent inclusion-exclusion for the centered absolute moment agree on all1116 profiles. Completeness comes from the proved concentration reduction, not a sampled partition cutoff.

Let V_small sum V_v over cusps having a rank at most8, and V_large sum V_v over all other cusps. The exact curvature lower bound minus the distance14 payment upper bound gives

**`4172075 V_large >= 60795045636 v_f +5580632 V_small +63584561900.`**

The coefficients are exact: beta_156−14*156=741 and the global slack is5*155=775, multiplied by82044596. V_v counts the actual edges incident to a cusp; it is not merely the cusp count or the factor count. Thus the multiplicity in this incidence inequality is explicit.

The [rank9 theorem](fresh_cusp_minimum_factor_rank9_dimension29_20261017.md) forces V_small>0 in dimension156, and the displayed balance forces V_large>0. Every hypothetical156-dimensional chamber must therefore have at least two cusps of these different types. Consequently:

* Every chamber whose each cusp has some factor of rank at most8 has dimension at most155.
* Every one-cusp chamber has dimension at most155: a large-only cusp would instead put it in the rank9 class, bounded by29.

The unrestricted156 theorem excludes all higher dimensions; compact chambers already satisfy29. No unrestricted155 conclusion follows, since mixed cusp types can satisfy this balance.

## Verification scope

The [checker](../scripts/fresh_small_cusp_defects_20261017.py) and [certificate](fresh_small_cusp_defects_20261017.json) use standalone Fraction arithmetic. They check128 general plus988 marked reduced profiles, two moment formulas, exact balance coefficients, eight small-factor cost identities and the complete rational five-wall example. Frozen proofs supply the localization cases, curvature, concentration and the two earlier dimension bounds. Their published tables and prior arithmetic are not reconstructed.

Generation and a fresh saved replay both pass. Syntax, all13 direct source hashes and11 local report links pass. The parent independently read the localization, nonincident-cusp argument, modified weight budget, deficit inequality, square-pyramid caveat and exact156 balance, finding no gap. That conceptual/source review does not duplicate the arithmetic or published finite tables. These are scoped incidence and historical restricted-dimension results, not an unrestricted improvement or external refereeing. Source, data and report are final/frozen.
