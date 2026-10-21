# Full review: the unique cubic B16 short-cap polynomial

21 October2026. **PASS at the stated necessary-condition scope.** I read
the complete145-line source and independently reconstructed its norm,
coefficient, interval and integral-basis arguments. No enumeration
program, scientific import, worker, test, Lean or TeX process ran. No
source, receipt, guide or frozen paper was edited.

Reviewed source:
`results/arithmetic17_b16_isolated_cubic_cap_unique_polynomial_candidate_20261021.md`,
SHA-256 `2fe5f3c6129bb1cf264d00de2f062e068ed94e66dca290be4490d79a2d328b22`,
5252B. Its sole direct source pin is unchanged:
`arithmetic17_b16_actual_i_short_opposite_isolation_and_degree_restriction_20261021.md`,
SHA-256 `51910bc233e16870d2322e66c9f2c2766076d757749088e7743277b8ddbed248`,
15626B. I authored that inherited reduction, but did not author this
new root-authored cubic theorem. This is an outside-author review of
the new argument, not an independent reproof of my own foundation.

The source assumes the actual degree-three branch only after the
accepted B16+actual-I short-opposite reduction. Thus B is integral,
ACTUAL K=Q(B), the physical root a lies in(4,441/16], and its two
nonphysical conjugates u,v lie strictly in(0,1/4). It neither imports
D19 quadraticity nor assumes a finite field list.

The positive integer norm is strictly below441/256<2, so auv=1
and a>16. The complete integral trace range is17..28 because
`16<S=a+u+v<449/16`. The remaining elementary symmetric coefficient
T is also a positive integer. With constant coefficient-1 the
minimal polynomial is exactly `X³-SX²+TX-1`.

I checked both necessary coefficient inequalities. At X1/4 the signs
of its three factors give `16T-4S-63<0`, hence the integral bound
`4T-S<=15`. At the small root u, u³<(1/4)u² gives

```text
0<-(S-1/4)u²+Tu-1.
```

To see explicitly why its discriminant is positive, multiply this
inequality by the positive4(S-1/4) and complete the square:
`T²-4S+1>(2(S-1/4)u-T)²>=0`. Integrality then yields T²>=4S.
The source's three trace ranges are exhaustive:17..20 gives T<=8
but T²>=68;21..24 gives T<=9 but T²>=84;25..28 gives T<=10
and T²>=4S>=100. Thus only S25,T10 is possible. This is a full
integer argument, not a report of an unexecuted coefficient search.

The five displayed polynomial signs recompute correctly. They locate
three distinct roots in(0,1/5),(1/5,1/4),(24,25), exhausting the
cubic's roots. The physical a>16 is the last. A positive finite old-I
contact has beta>=1/2 and would force a<=1323/64<24, a contradiction.
Hence beta0: the cap is perpendicular to the **old I wall**, with no
claim about its opposite across I or any other actual wall.

The rational-root test excludes both possible rational roots+/-1.
The cubic discriminant re-expands to473=11*43, which is squarefree.
The standard discriminant/index identity therefore proves
`O_K=Z[B]` and field discriminant473. Finally N(B)=1 and
`N(4B-1)=-64F(1/4)=3` have the correct signs and scale.
No denominator from a raw-root time coefficient is used in these
field computations.

No correctness blocker or missing coefficient case was found. The new
result singles out one cubic field and its physical generator, with
zero old-I contact. It supplies no actual-I parameter classification,
time compatibility or geometric contradiction for that field; the
quadratic branch is untouched. The B16+I forest, seven arithmetic
groups and global ceiling17 remain unchanged. A polynomial with the
right coarse embeddings is not an actual cap or a polytope. These
limits and the absence of runtime/Lean credit are stated honestly.
