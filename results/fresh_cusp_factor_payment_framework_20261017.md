# Cusp-factor payments and a positive-semidefinite spectator lemma

This package proves a restricted factor-weighted dimension inequality and an exact Schur-complement monotonicity lemma. The restriction is essential: **every irreducible affine factor at every cusp has rank at least9**. Under that restriction, at most two factors per cusp imply dimension at most117. This does not improve the unrestricted bound163. The determinant lemma is intended to support a finer argument without discarding nonsimple cusps.

The [checker](../scripts/fresh_cusp_factor_payment_framework_20261017.py) and [saved data](fresh_cusp_factor_payment_framework_20261017.json) certify the exact algebra below. Finite graph/model checks corroborate the all-rank proofs; they do not supply a cutoff argument.

## Cusp factors and actual two-face corners

Let P be a finite-sided finite-volume hyperbolic Coxeter polytope of dimension n. Write an ideal vertex diagram as a disjoint union of irreducible affine diagrams with ranks r_vi, where sum_i r_vi=n−1. Its Euclidean link is the product L_v=product_i Delta^(r_vi). Put

` s_v = number of factors,       V_v = product_i (r_vi+1). `

Thus V_v is the number of edges at v and the link has (n−1)V_v/2 edges. These are the standard finite-volume cusp-link and simplicity-at-edges facts used by [Khovanskii, §6 and Theorem10](https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf), also recorded in the [frozen truncation note](fresh_finite_volume_simplex_truncation_20261017.md).

Every connected affine subdiagram is an entire irreducible factor of its actual cusp. Indeed its positive kernel gives a nonzero null vector: pairing its positive root combination with a strict chamber center is negative. Acuteness makes that vector feasible against every other wall. A feasible ideal point is extreme in the projective closure, by strict convexity of the closed Klein ball, hence is an ideal vertex. Its affine subdiagram lies in one irreducible cusp component, whose proper principal subdiagrams are positive definite. The published containment of every parabolic subdiagram in its unique ideal vertex diagram appears in [Felikson–Tumarkin, §1](https://arxiv.org/pdf/math/0502413).

An actual two-face has an elliptic supporting set W of n−2 walls. Each of its edges adds one side normal, giving a positive-definite set of n−1 walls because the edge has an interior hyperbolic point. At an ideal corner, the corresponding link edge varies in exactly one simplex factor i. W retains r_vi−1 nodes in that factor and r_vj nodes in every other factor. Adding its two marked sides restores the full affine factor i; all other factors have one node deleted and are elliptic. In particular the marked pair lies in the varying factor. Every pair of nodes of factor i occurs in exactly V_v/(r_vi+1) such link-edge fibers.

This description does not assume that the cusp itself is simple. An independent internal conceptual read checked this corner/fiber multiplicity and the closed-edge argument below; it did not repeat the scalar or graph checks.

## A restricted distance-seven payment theorem

Assume from now on that every r_vi>=9. There is no affine subdiagram of order at most9. A minimal indefinite diagram is Lannér or quasi-Lannér of order at most10; the quasi-Lannér case contains an affine proper subdiagram of order at most9 and is therefore excluded. Every minimal indefinite diagram is consequently Lannér of order at most5. The finite classification and its use here are inherited explicitly from the [affine-exclusion proof](fresh_simple_cusp_affine_exclusion_20261017.md).

Assign a two-face corner weight1 when its two side nodes are at graph distance at most7, and weight0 otherwise. For a triangle, a minimal indefinite subdiagram of its star contains all three sides. Its order is at most5. A tree joining the three marked nodes has at least two marked-to-marked paths avoiding the third node, each of length at most3; these give two paid corners. Cyclic minimal diagrams only supply additional short paths. Thus every triangle receives at least2.

For a quadrilateral, W plus either opposite pair is indefinite. If it were positive definite, the usual acute projection of the strict center would give a point on both opposite closed sides. If it were positive semidefinite and singular, its positive affine kernel would give a nonzero feasible ideal point on both opposite closed sides. Both contradict their disjointness, including when the quadrilateral has ideal vertices. Let L and M be minimal indefinite subdiagrams for the two opposite pairs. Each contains both of its marked nodes, has order at most5, and has distance at most3 from a white node to a nearest marked node. The two subdiagrams cannot be orthogonal, because that would give two negative directions. Their overlap or a joining bond therefore supplies a path of length at most3+1+3=7 between adjacent sides, avoiding the other two marked nodes. Every quadrilateral receives at least1. This is the small-core proof of the frozen distance theorem, using only actual edge/corner supports rather than simple ideal vertices.

At a finite vertex the total weight is at most7n. At an ideal vertex, a factor of order m=r+1>=10 is classical, and its number of pairs at distance at most7 is at most7m. The all-rank path, cycle, single-fork and double-fork formulas are proved in the [distance package](fresh_simple_cusp_distance_bound_20261017.md). The current checker verifies the245 inherited literal graph rows and their source hashes; these checks are not the justification for arbitrary rank. Multiplying by the fiber count yields the exact budget

` C(v) <= sum_i [V_v/(r_vi+1)] * 7(r_vi+1) = 7 s_v V_v. `

## The exact weighted inequality

Let v_f count finite vertices and define

` beta_n = n(n−6)/8                 (n even), `

` beta_n = (n−1)(n−5)/8             (n odd). `

For the original polytope, not its truncation, Khovanskii's reduction yields a strict average-face bound for a bounded polytope simple at edges, as justified immediately below, and gives

` D5 = sum_(two-faces F) (5−number_of_sides(F)) `

`      > beta_n v_f + (beta_n/n) sum_v V_v. `

Here the total two-face/vertex incidence is binom(n,2)v_f+(n−1)sum_v V_v/2. Here is the strictness step rather than a tacit upgrade of the published weak bound. Choose a generic linear function on the bounded projective polytope. For each vertex w, let E_w be the number of edges in its simple vertex figure and E_w^0 the number of those edges not cut by the hyperplane through w parallel to a level hyperplane. Theorem9 gives the exact identities I=sum_w E_w and2f2=sum_w E_w^0. Theorem8, as applied in Theorem10, gives2E_w/E_w^0<=eta_n, where eta_n=4(n−1)/(n−2) for even n and4n/(n−1) for odd n. All denominators are positive. At a global minimum vertex the cutting hyperplane misses the vertex figure, so E_w^0=E_w and this ratio is2<eta_n. Thus I/f2, the positively weighted average of these local ratios, is strictly less than eta_n. This proves the displayed strict lower bound for D5=5f2−I. The argument applies to the bounded projective closure, including ideal vertices. The checker verifies the resulting coefficient identities and parity expressions. The corner payments give

` D5 <= 7n v_f + 7 sum_v s_v V_v. `

Consequently

` beta_n/(7n) < [n v_f + sum_v s_v V_v] / [n v_f + sum_v V_v]. `

This is the proved factor-weighted inequality. If every s_v<=s, the right side is at most s, and hence n<=56s+5. At N=56s+6, beta_N=7sN, which the strict inequality excludes without any additional finite-vertex assumption. The all-rank threshold check is

` 8(beta_(N+2t)−7s(N+2t)) = 4t^2+2Nt, `

` 8(beta_(N+1+2t)−7s(N+1+2t)) = 4t^2+(2N+4)t+N+6. `

Thus the rank-at-least9, at-most-two-factor class has n<=117. These are restricted consequences, not claims about arbitrary finite-volume polytopes. The frozen all-simple bound29 is already stronger for s=1.

## Orthogonal positive-definite spectators improve a scalar contradiction

The following identity has no factor-rank restriction. After eliminating the supporting whites of an affine marked corner, suppose alpha>0 and h,k are nonzero, and its four-root Schur matrix is

```
[ alpha  −alpha   B−h   −A  ]
[−alpha   alpha  −B     A−k ]
[ B−h     −B      C     −E  ]
[−A        A−k   −E      D  ].
```

Its determinant is

` (hk−kB−hA)^2 − alpha(Ck^2+Dh^2+2Ehk). `

Additional elliptic white spectators orthogonal to the old affine factor may couple to the two new roots c,d. Eliminating them subtracts a positive-semidefinite matrix [[p,q],[q,r]] from the c,d block, so C becomes C−p, D becomes D−r and E becomes E+q. With u=1/h and v=1/k, the normalized determinant changes by

` alpha(pu^2−2quv+rv^2) = alpha ||u z_c−v z_d||^2 >= 0. `

The equality uses any Gram realization of that positive-semidefinite correction; h,k are unchanged. Therefore every strict positive determinant contradiction of this form survives such spectator attachments. For an actual quadrilateral the complete supporting set has n−2 independent roots and the four residual normals lie in a Lorentz three-space, so the determinant must vanish. This explains why a strict scalar contradiction still excludes the configuration after the spectators are restored.

This identity transfers scalar portions of the [affine-cycle proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) and [affine-tree proof](fresh_simple_cusp_affine_tree_completion_20261017.md). It does not itself transfer their entire payment theorem. A generic finite-component payment can change upon extension (for example, an abstract E7 diagram embeds in E8), while role-compatible spectator extensions in the actual proof require a separate literal check. This sentence does not assert that the specific D6-to-E7 configuration admits such an extension through its permitted new-root position. Also, a small affine diagram that was forbidden under simple cusps can now be a genuine cusp factor. Those separate cases require a new argument before a3s_v V_v budget can be asserted.

## Two factors do not force absence of small affine diagrams

There is an exact local cusp model in every dimension m+1, m>=3, with link Delta^1 times Delta^(m−1). It shows why a factor-count hypothesis alone cannot replace the minimum-rank hypothesis above.

Use a Lorentz space with null nu,eta satisfying nu.eta=−1, an orthogonal unit e, and an orthogonal A_(m−1) root system alpha_1,...,alpha_(m−1). Put theta=sum alpha_i. The normals are

` e, nu−e; alpha_1,...,alpha_(m−1), nu−theta. `

The two components are affine A1 and affine A_(m−1), mutually orthogonal, and both positive null relations sum to nu. Their Gram rank is m and their normal-span rank is m+1. A strict timelike center is

` H=m nu+eta−e/2−sum_(i=1)^(m−1) [i(m−i)/(2m)] alpha_i. `

Its products with the first pair are−1/2 and with each remaining normal−1/m, while

` H^2=−2m+1/4+(m^2−1)/(12m)<0. `

The horospherical link is the bounded product of the two standard alcoves. These formulas prove the all-rank local realization;22 exact matrix checks corroborate them. No finite-volume completion in arbitrary dimension is claimed.

## Verification and scope

Generation passes:245 inherited graph rows with their source manifest,171 scalar incidence checks, eight exact factor-count thresholds, the symbolic four-root determinant and positive-semidefinite correction, and22 rational cusp models. A fresh saved replay also terminates with PASS; syntax, all eight direct source hashes and all seven local report links pass. The separate geometric review described above is limited to the factor-corner and closed-edge steps. A separate root mathematical/source review passes the strict Khovanskii reduction, factor geometry, distance-seven scope and positive-semidefinite identity; see the [root review](fresh_cusp_factor_payment_framework_20261017_root_review.md). That review does not duplicate the producer arithmetic. The source, data and this report are final/frozen. This package does not assert a full fine-weight extension, an unrestricted dimension improvement, or novelty of the restricted consequence.
