# A full-host obstruction from a repeated three-dimensional vertex configuration

The conductor-nine host

\[
L=\langle-(T-1)\rangle\oplus
\begin{pmatrix}2&-T\\-T&2\end{pmatrix}\oplus\langle1\rangle,
\qquad T^3-3T-1=0,\quad T>1,
\]

is **nonreflective**: its full integral reflection chamber in H3 has infinite volume. This is a full-host result, rather than another finite-prefix escape. The necessary-face argument below consequently rules out a compact chamber in the corresponding **H9 host `⟨−(T−1)⟩+B⁴+⟨1⟩`**, the previously searched coefficient `α=Tβ`, where `β=4−T²` and `Tβ=T−1`.

This does not rule out H9 Coxeter polytopes generally, other conductor-nine coefficients, or other lattices over the field. No worldwide novelty claim is made for this nonreflectivity result.

## Exact infinite-order isometry

Let G be the displayed rank-four Gram matrix. The following matrix acts on column vectors in its coordinate basis:

```
g = [ 43+152T+81T²,  -6−18T−9T²,   -17−57T−30T²,  -16−55T−29T²;
      71+254T+136T², -8−29T−16T²,  -27−95T−51T²,  -26−92T−49T²;
      73+257T+136T², -8−30T−16T²,  -27−96T−51T²,  -27−93T−49T²;
      19+68T+36T²,   -3−8T−4T²,   -6−25T−14T²,   -7−24T−13T² ].
```

Independent Fraction arithmetic verifies `gᵀGg=G`, and both g and `G⁻¹gᵀG` are integral. Thus g is a lattice automorphism. Its trace is `1+3T+T²>4`, so it has infinite order: a finite-order real four-dimensional matrix has trace at most four.

Use the rational point

\[
p=(1,55/32,29/16,5/16).
\]

Its image is

\[
gp=\left(
\frac{-100+18T+67T^2}{32},
\frac{6+103T+24T^2}{32},
\frac{60+76T+24T^2}{32},
\frac{25+46T-10T^2}{32}
\right).
\]

Both points are future timelike, with identical square `(6951−3702T)/512<0`. The full-chamber membership check below proves that they lie in the same open reflection chamber C. Since g permutes the full integral mirror arrangement, it follows that `gC=C`.

## Why the bounded check proves membership in the full chamber

The [earlier field and norm audit](fresh_cubic9_dihedral_report_20261013.md) proves `O=Z[T]`, class number one, and that totally positive units are squares. Here the dual exponent is β. Consequently every primitive reflection can be rescaled by a unit so that its positive norm belongs to `{1,2,β,2β}`. The complete zero-height chamber is `I₂(18)×A1`, with 38 oriented roots. The verifier checks that p and gp have the same nonzero sign on every zero-height root.

Orient a positive-height root `r=(x,y)` with `x>0`, and write its norm as q and its height as `h=(T−1)x²/q`. For any future timelike point `u=(u₀,u_s)`, if `(u,r)≥0`, spatial Cauchy gives

\[
(T-1)^2x^2u_0^2
\leq(q+(T-1)x^2)\,u_s^2,
\qquad
h\leq\frac{(T-1)u_0^2}{-u^2}-1.
\]

The last bounds for p and gp are approximately **68.438** and **1969.772**. Their strict inequalities below 3000 are checked exactly; the decimal values are explanatory only.

It remains to justify restricting enumeration to dominant roots. If r pairs positively with a simple zero-height wall a, reflecting r in a preserves its norm, time coordinate and height. Because u lies strictly inside the zero chamber,

\[
(u,s_a r)=(u,r)-2(r,a)(u,a)/a^2>(u,r).
\]

The finite zero-height reflection group folds r into the closed dominant chamber. Therefore any mirror with `(u,r)≥0` would yield a dominant root of the same height with nonnegative pairing with u. This applies separately to both endpoints. It also excludes a mirror through either endpoint.

A complete independent enumeration through height 3000 has **142 time/norm controls and 136 dominant positive-height candidate roots**, before wall selection. The verifier checks **every one of those 136 roots** strictly negatively at **both** p and gp. It does not merely test the fifteen accepted walls. Thus no positive-height mirror separates either endpoint from the chosen chamber adjacent to `(1,0,0,0)`, and the zero-height sign check places both endpoints in the same sector. This proves `p,gp∈C` for the full mirror arrangement.

## Independent enumeration and the global conclusion

The producer's higher prefix has fifteen walls and a spacelike escape. Its role here is to propose vertex configurations and provide a candidate census; the proof uses a different complete enumerator.

For any totally positive spatial cap c, the embedding bounds `x_e²≤c_e` imply `Tr(x²/c)≤3`. The [independent enumerator](../scripts/fresh_cubic9_ellipsoid_audit_20261013.py) forms the rational coefficient matrix of this trace form, checks its exact LDL decomposition, and recursively enumerates every integral vector in the ellipsoid. Integer square-root bounds include both endpoints. Final exact embedding tests remove the excess. Time coordinates use unpruned certified boxes, and spatial blocks use a Cartesian join. No producer field arithmetic, numerical geometry, or interval alphabet is imported. The complete enumeration is repeated inside the [global verifier](../scripts/fresh_cubic9_global_nonreflectivity_verify_20261013.py), followed by the all-candidate endpoint checks. A corrupted isometry and a nonregular endpoint are rejected.

The other two real embeddings of G are positive definite, so its integral Lorentz isometry group is discrete. Its infinite-order element g stabilizes C. Proper discontinuity supplies infinitely many disjoint g-translates of a sufficiently small open ball inside C. These have equal positive volume, proving that C has infinite volume and the integral reflection subgroup has infinite index. This establishes nonreflectivity of the full H3 host.

## Consequence for the H9 search

The complete zero chamber of `⟨−(T−1)⟩+B⁴+⟨1⟩` gives the finite vertex diagram `I₂(18)^4×A1`. If its full chamber were compact, removing three `I₂(18)` components would give a compact Coxeter H3 face. The face criterion applies because these spherical components are neither `A_n` nor `D5`. No connected spherical extension of `I₂(18)` has rank greater than two, so its surviving face walls are orthogonal to the removed blocks. This is an application of [Allcock's Theorems 2.1–2.2, attributed there to Borcherds](https://arxiv.org/pdf/0903.0138).

Orthogonality to `−e₁` and `D(e₁−e₂)` forces both coordinates in each removed B block to vanish. Their orthogonal complement is exactly the H3 lattice above, and each face reflection restricts to an integral reflection of that lattice. A compact reflection group generated by these face walls would be a finite-covolume subgroup of the full H3 reflection group. That is impossible for the infinite-volume full group just certified. Hence this **particular H9 host has no compact full integral reflection chamber**, and no compact integral reflection subgroup.

Moreover the admissible form is anisotropic over K, since it is positive definite at another real place. Its ambient arithmetic quotient is cocompact; a finite-covolume reflection subgroup would have finite index and would also be cocompact. Thus there is no finite-covolume integral reflection subgroup in this H9 host either. The same face argument applies to `⟨−(T−1)⟩+B^k+⟨1⟩` for every `k≥1`; it does not apply automatically to other coefficients or positive blocks.

## Artifacts and replay

- [Global certificate, matrix, endpoints and all root pairings](fresh_cubic9_global_nonreflectivity_20261013.json).
- [Full independent height-3000 enumeration](fresh_cubic9_h3_h3000_full_enumeration_audit_20261013.json).
- [Saved prefix and vertex census](fresh_cubic9_dihedral_h3_theta_beta_h3000_20261013.json).
- [Vertex-isometry proposal search](../scripts/fresh_cubic9_vertex_symmetry_20261013.py).
- [Independent source review and final replay](fresh_cubic9_stabilizer_source_review_20261013.md).
- [Trace-ellipsoid regression comparisons](../tests/test_fresh_cubic9_ellipsoid_20261013.py).

```sh
python scripts/fresh_cubic9_global_nonreflectivity_verify_20261013.py --verify
```

PARI was used only to propose square roots in the vertex-isometry scout. The global verifier uses standard-library rational arithmetic and imports no PARI, SymPy, NumPy or project arithmetic implementation. The successful result closes one family of arithmetic hosts; the search for a new-dimensional Coxeter example remains open.
