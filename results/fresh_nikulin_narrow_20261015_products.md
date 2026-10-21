# A product-cusp separation bound for a narrow facet

**Conditional theorem.** If a Nikulin narrow facet of a finite-sided finite-volume real hyperbolic Coxeter polytope (P\subset\mathbb H^n) contains an ideal vertex, then (n\le225).

This treats products of affine cusp components, including arbitrarily many factors of small rank. It is a further explicit step in reconstructing the historical narrow-facet argument. It does **not** yet give an unrestricted bound of 225: a chosen narrow facet may have only finite vertices. The latter case, which Nikulin's outline proposes treating combinatorially, is still open in this reconstruction. No novelty or current-best claim is made.

The [checker](../scripts/fresh_nikulin_narrow_20261015_products.py) and [certificate](fresh_nikulin_narrow_20261015_products.json) verify the finite exceptional constants and the all-rank scalar identities. The argument uses the inverse and coordinate formulas proved in the separate [classical-cusp package](fresh_nikulin_narrow_20261015.md), whose files remain frozen. The published narrow-facet input is [Nikulin, ICM 1986, printed p.655](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

## Setup and the common height ratio

Use outward normals of square 2 in Lorentz signature ((n,1)). At an ideal vertex (v), a horospherical link is a compact Euclidean Coxeter polytope, hence a product of affine Coxeter simplices. Let their ranks be (r_a\ge1), with

\[
N:=n-1=\sum_a r_a.
\]

Every factor has (r_a+1) nodes. Deleting one node in each factor specifies a vertex of the horospherical link, hence an outgoing edge of (P). Its retained support (S) has (N) roots and positive-definite block-diagonal Gram matrix. The standard Euclidean reflection classification, including the parallel two-wall affine (A_1) case, is the only diagram classification imported here.

Choose two such deletion tuples as described below, and let (x,y) be terminating walls at the other endpoints of their edges. The next section ensures these can be chosen to meet the narrow facet in actual ridges. Their contact rows against the initial cusp roots are nonnegative after negating pairings. In a factor (a), write these rows as (p_a,q_a), its positive null coefficients as (\delta_a), and

\[
c_a=\sum_j\delta_{a,j}\alpha_{a,j},\quad
\sigma_a=\delta_a\cdot p_a=-(x,c_a),\quad
\tau_a=\delta_a\cdot q_a=-(y,c_a).
\]

All (c_a) are positive multiples of the same future null vector. Therefore the ratio

\[
q:=\sigma_a/\tau_a>0
\]

is **the same in every factor**, although their individual scales need not agree. Put (M=\max(q,1/q)\ge1).

The cusp quotient is an orthogonal direct sum of the Euclidean factor spaces. Let (E_a) be the squared norm of the factor component of (x/\sigma_a-y/\tau_a). The Lorentz identity from the classical package becomes

\[
\boxed{-(x,y)=\frac12\sum_a\sigma_a\tau_aE_a-q-1/q.}
\]

The relative null scales cancel in each summand. This is why the proof does not require commensurable cusp scales, arithmeticity, or a prescribed horosphere.

## Choosing terminating walls that really meet the narrow facet

Let the narrow facet have normal (\alpha_s). In its affine factor, both deletion tuples avoid (s), so both outgoing edges lie in the narrow facet. At a finite second endpoint, every additional wall has finite angle with (\alpha_s).

At an ideal second endpoint, let (T_b) be its affine components. Its full rank is again (N). The retained edge support (S) is positive definite with exactly (N) nodes, so it contains at most (\operatorname{rank}T_b) nodes from each component. Equality of the total counts forces equality in every component. Thus exactly one node in each (T_b) is missing from (S), and these are additional terminating walls.

If the component containing (\alpha_s) has rank at least 2, choose its missing node. Every two nodes in such an irreducible affine diagram have finite angle, so this wall meets (\alpha_s). If that component has rank 1, its missing node is parallel to (\alpha_s); instead choose a missing node of another component. That wall is right to (\alpha_s). Another component exists for (n\ge3), the only range relevant to the claimed contradiction.

In either case the chosen terminating wall and (\alpha_s) have positive-definite rank-two Gram. Projecting a strict interior point of (P) onto their common orthogonal complement gives a timelike point still strict against every other wall, because the inverse Gram is nonnegative and all distinct wall pairings are nonpositive. Hence they share an **actual ridge**, not merely a point at infinity. This proves that (x,y) belong to the exact framing-wall set controlled by the narrow-facet theorem.

## At most eight factors can carry retained nonright contacts

At either second endpoint, the Gram matrix of (S\cup\{x\}) is positive semidefinite. Thus

\[
p_S^T G(S)^{-1}p_S\le2.
\]

A positive-definite matrix with diagonal 2 has inverse diagonal at least (1/2), by its one-node Schur complement. A positive-definite Coxeter Gram is a symmetric (Z)-matrix, so its inverse is entrywise nonnegative. Every nonzero retained contact has magnitude at least 1. Consequently at most four entries of (p_S) can be nonzero. The same holds for (y).

Call a factor **spoiled** if either terminating wall has a nonright contact with any retained node in that factor. There are at most eight spoiled factors. Also mark the factor containing (s), regardless of its contacts. Thus at most nine factors are marked. In all unmarked factors each contact row is supported only on the deleted node. A terminating wall cannot have zero contact with an entire factor, since its null sum is a positive multiple of the cusp vector and its height is positive.

## Point-contact energy in every affine type

For a factor with Gram (G), null vector (\delta), and two chosen omitted nodes (i,j), define

\[
R_{ij}=\left(e_i/\delta_i-e_j/\delta_j\right)^T
G^+\left(e_i/\delta_i-e_j/\delta_j\right),
\quad L_{ij}=\min(\delta_i^2,\delta_j^2)R_{ij}.
\]

The definition is unchanged by rescaling (\delta). In an unspoiled factor the two contact magnitudes are at least 1. Since their null-height ratio is (q),

\[
\sigma_a\tau_a\ge
\max(q\delta_j^2,\delta_i^2/q)
\ge M\min(\delta_i^2,\delta_j^2).
\]

The omitted nodes can be chosen so that (L_{ij}\ge2r/9) for **every** irreducible affine type of rank (r):

| Type | Chosen nodes | Exact or lower point-contact energy (L_{ij}) |
|---|---|---:|
| (\widetilde A_r, r\ge2) | opposite cycle nodes | (\lfloor(r+1)^2/4\rfloor/(r+1)\ge r/4) |
| (\widetilde B_r) | a fork leaf and the short end | (r/4) |
| (\widetilde C_r) | the two short ends | (r/2) |
| (\widetilde D_r) | fork leaves at opposite ends | (r/4) |
| (\widetilde A_1) | its two nodes | (1/2) |
| (\widetilde E_6) | two outer leaves | (4/3=2r/9) |
| (\widetilde E_7) | equal-weight nodes specified in JSON | (2) |
| (\widetilde E_8) | equal-weight nodes specified in JSON | (2) |
| (\widetilde F_4) | its two path ends | (1) |
| (\widetilde G_2) | its two path ends | (2/3) |

The classical values follow from the explicit coordinates in the first package. For (B,D), the normalized point difference has (r) coordinates of magnitude (1/2); for (C), it has (r) coordinates of magnitude 1 and the minimum squared null coefficient is (1/2). The finite exceptional calculation is exact: the checker uses ((G+\delta\delta^T)^{-1}) on the orthogonal complement of (\delta), where it agrees with (G^+). The certificate records each full matrix and node choice, rather than referring only to diagram names.

For every unmarked factor we use these choices. They do not need to avoid a narrow-facet node, since that factor was separately marked.

## Every marked factor loses at most 2

We claim that in each marked factor of rank (r), the omission choices can be made so that

\[
\sigma_a\tau_aE_a\ge M(2r/9-2).
\]

For (r\le9), simply use nonnegativity of energy. If it is the narrow-facet factor, choose any node other than (s) for both omissions; they may coincide. Since (2r/9-2\le0), this proves the claim without a low-rank search.

For (r\ge10), only the classical families occur. Use the following choices and Schur-localization estimates. Importantly, the choices are made **before** seeing the terminating contact rows.

* For (A_r), take opposite cycle nodes, rotating both to avoid (s) when necessary. Deleting a node leaves (A_r); the inverse-diagonal test localizes either row within two cyclic nodes of its deletion. The two long-arc energy estimate gives (E_a\ge(r-8)/4). All null coefficients equal 1, so (\sigma_a\tau_a\ge M).
* For (B_r), take a left fork leaf and the right short end. Deletion leaves respectively (B_r) and (D_r). In the first row only nodes (0,1,2) can be contacted; the other row is supported within the last three nodes. If the short end is (s), delete the adjacent interior node instead; the retained diagram is (D_{r-1}+A_1), and the second row uses at most the last four nodes. If the chosen left leaf is (s), use the other left leaf. Thus at least (r-5) intervening spatial coordinates have normalized difference (1/2). Every allowed null coefficient is at least 1, giving (\sigma_a\tau_aE_a\ge M(r-5)/4).
* For (D_r), take fork leaves at opposite ends, switching to the other leaf if one is (s). Each deletion leaves (D_r); inverse-diagonal localization allows only the adjacent end cluster. At least (r-4) coordinates have normalized difference (1/2). The weaker bound (M(r-5)/4) suffices.
* For (C_r), take its two ends. If one is (s), move that deletion inward by one node; the retained diagram there is (A_1+B_{r-1}). Localization leaves at least (r-5) coordinates with normalized difference 1. The allowed null coefficients are at least (1/\sqrt2), giving (\sigma_a\tau_aE_a\ge M(r-5)/2).

These are applications of the same closed (A,B,D) inverse formulas and cumulative-coordinate calculations already proved in the first package. In all cases (r\ge10) meets the needed rank thresholds, including (D_{r-1}) in the modified (B) choice. The elementary identities

\[
\frac{r-8}{4}-(2r/9-2)=r/36,\qquad
\frac{r-5}{4}-(2r/9-2)=r/36+3/4,
\]

and

\[
\frac{r-5}{2}-(2r/9-2)=5r/18-1/2>0
\]

prove the claimed uniform loss bound. No assumption about the contacts' maximum labels or magnitudes is used.

There is no circular choice here: in every non-narrow factor use the point-contact omissions from the preceding table. If that factor subsequently turns out spoiled, the same classical end choices satisfy these localization estimates; the low-rank estimate needs no change. Only the narrow factor uses the stated modified choices, and it is already included in the loss allowance whether or not it is spoiled.

## Summing the energies

There are at most nine marked factors, so

\[
\sum_a\sigma_a\tau_aE_a\ge M\left(\frac{2N}{9}-18\right).
\]

Since (q+1/q\le2M), the Lorentz pairing identity yields

\[
-(x,y)\ge M\left(\frac N9-11\right).
\]

If (n\ge226), then (N\ge225) and the right side is at least 14. But both (x,y) are actual framing walls of the narrow facet, whose pairwise hyperplane distance is strictly below (\operatorname{arccosh}7). For square-two normals this contradicts (-(x,y)\ge14). Therefore (n\le225) under the stated ideal-vertex hypothesis.

The bound remains conditional. It proves that in dimension at least 226 a Nikulin narrow facet, if one existed in a finite-volume Coxeter polytope, would have only finite vertices. Turning that remaining case into an unconditional dimension bound requires a separate argument; it cannot be obtained by silently declaring the facet itself a Coxeter polytope.

```sh
python scripts/fresh_nikulin_narrow_20261015_products.py --verify
```
