# Full outside-author review: cubic B16 cap and the actual-I strip

21 October2026. **FULL MATHEMATICAL PASS at the stated cubic-branch scope.** I read the complete175-line new source and both complete direct written inputs, and independently reconstructed the root bounds, integral-basis strip and same-vector contradiction. No mathematical process, enumeration, worker import, symbolic computation, test, Lean or TeX job ran. No reviewed source was changed.

Reviewed source: `results/arithmetic17_b16_isolated_cubic_cap_actual_i_strip_exclusion_candidate_20261021.md`, SHA-256 `f89214ff279908fa2d1ebb1621102adab8a4c17785418c2519d00b2d3811bc45`,6975B. Its two literal direct pins match current bytes:

| Input under `results/` | SHA-256 |
| --- | --- |
| `arithmetic17_b16_actual_i_short_opposite_isolation_and_degree_restriction_20261021.md` | `51910bc233e16870d2322e66c9f2c2766076d757749088e7743277b8ddbed248` |
| `arithmetic17_b16_isolated_cubic_cap_unique_polynomial_candidate_20261021.md` | `2fe5f3c6129bb1cf264d00de2f062e068ed94e66dca290be4490d79a2d328b22` |

The latter's separate full outside-author review is `review_arithmetic17_b16_cubic_unique_polynomial_full_scope_foundations_20261021.md`, SHA-256 `3158e84903f36e56bd29a8326454315de99c7e1ee26cc8090d93414707ae76b4`. Its PASS label was not substituted for checking the interfaces or the new proof. I contributed older B17 normalization/payment foundations used transitively and previously reviewed51910bc2; I did not author the new cubic polynomial theorem or this actual-I strip theorem.

**Actual scope and field, source lines10–39.** The prior reduction supplies an actual isolated dotted short cap at B16+I, independent strength at6..14, and ACTUAL K=Q(B) of degree2 or3. The second input, restricted to degree3, forces beta0 and `F(B)=B^3-25B^2+10B-1=0`, with physical root a in(24,25), two other roots in(0,1/4), norm1 and `O_K=Z[B]`. Its coefficient argument covers the entire cubic branch: the norm forces constant coefficient−1; the trace lies17..28; `4T-S<=15` and `T^2>=4S` force S25,T10. The displayed discriminant473=11*43 is squarefree, so the integral-basis conclusion follows from the discriminant/index identity. These are properties of the actual field, not a selected literal lattice.

For the ACTUAL isolated-I opposite, the invariant two-cycle with the old I normal is `J=4(1+h)`. Thus J is integral over that same K, physically4<J<6. At each nonphysical embedding, positivity of the old-I/c_I Gram and its positive norm factors gives0<J_sigma<4. The already locked K=Q(B) leaves no extra field embedding to omit. The published automatic admissible-form input is retained from51910bc2 (Vinberg's structural reflection lemma, recorded in Agol–Belolipetsky–Storm–Whyte); the source does not introduce a restrictive first-type assumption or claim to prove that external result.

**Root localization and basis, lines43–66.** Direct integer expansion verifies the four displayed numerators−333,2193409,1716121,−144869 with their stated positive denominators. Together with the known three-root ordering they give `.184<u<.1845`, `.2205<v<.222`, `a=25-u-v>24.5935` and v−u>.036. Division of F(B)=0 by the norm-one unit B gives `eta=1/B=B^2-25B+10`. The basis change from(1,B,B²) to(1,B,eta) is triangular with determinant1, so every integral J has unique integer coefficients l,m,n. Since auv=1, subtracting the embeddings and dividing by a−u or a−v gives exactly m−nv and m−nu, with the inverse-root term's minus sign as written.

**All integer strips, lines68–102.** Physical J lies strictly above BOTH nonphysical values, and both differences are<6. Each root denominator is>24.368; `24.368*.247=6.018896>6` proves the stated upper. Hence both m−nv and m−nu lie in(0,.247). Their difference has magnitude |n|(v−u), not twice the interval length, and7*.036>.247 proves |n|<=6.

The source exhausts every integer in that range. For n1..4 it requires simultaneously m>=1 and m<.985; for n5..6, m>=2 and m<1.354; n0 is impossible. For n=−k, k1..4, the second strip makes m>=0, excluding k2..4 by the first strip and leaving only m0,k1. For k5..6 the two strips force m=−1; k5 then has−1+5u<0 and k6 has−1+6v>.323>.247. Thus m0,n−1 is the sole possibility. Physical4<l−1/a<6 leaves l5 or6. The u-embedding rejects l5, so J=6−1/B is necessary and unique. Its nonphysical values are indeed in(0,2), since1/6<u,v<1/4; the proof correctly does not claim this coarse strip is itself inconsistent.

**Actual-I energy, lines106–117.** Substitution gives h=(2−1/a)/4>47/96. If r is the original I half-tangent square, then h=4r/(1−r)², strictly increasing for0<r<1. At r=99/1000 the value is396000/811801, and `47*811801-96*396000=138647>0`. Hence r>.099 and the SAME original-vector component energy e_I>=r is>.099. No auxiliary point's energy is used here.

**Paid reference, lines119–160.** The first-six-linear reference is the UNLIFTED one of51910bc2: all independently strong knots6..14, energy E0>.9084−.01814=.89026, endpoint f14<3.311 and available P14 coefficient>.13. It assumes neither strong5 nor any conclusion about the short cap. Beta0 gives actual short-edge cosine C=4gamma=2sqrt(a)>9. Its direction is P16/4; actual farthestness therefore implies S16>=4sqrt((C−1)/(C+1))>7/2. The weaker7/2 is valid because64/5>49/4.

For mu=.13, the added residual P16−P14 is supported only on coordinates15,16, so it is orthogonal to the complete reference and has squared norm2. The lower pressure on W is justified by first removing mu from the positive P14 coefficient and then using the ACTUAL S16 bound. This gives `x.W>=E0+mu(S16-f14)` with all other multipliers nonnegative; it is not inferred by subtracting two lower bounds. Consequently

```text
e_B >= 2x.W-||W||^2
    > .89026+.26*(3.5-3.311)-.0338
     =.90560.
```

The components are orthogonal in the original vector, so adding e_I>.099 yields1.00460>1. The I energy was used earlier to derive a universal upper edge length and hence a necessary polynomial range; that logical consequence is not a second numerical debit. The final calculation uses the raw E0 reference, with no earlier I lift or joint-I gain added, and counts the actual I energy once.

No omitted embedding, integer coefficient, sign, payment or own/contact case was found. The result excludes the ENTIRE cubic actual short-cap branch only. Combined with51910bc2, the isolated dotted quadratic actual-field branch remains; the whole B16+I family and all global ceilings remain open/unchanged at this source's boundary. There is no actual example, finite-field census, new runtime or Lean credit. No correction is requested.
