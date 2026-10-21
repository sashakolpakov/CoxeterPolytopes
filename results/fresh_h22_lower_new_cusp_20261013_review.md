# Source review and fresh replay of the two new H22 cusp defects

PASS. The [incidence checker](../scripts/fresh_h22_lower_new_cusp_20261013.py)
and its [saved result](fresh_h22_lower_new_cusp_20261013.json) were reviewed and
freshly replayed from the root task after the construction agent's run.

For each of the two norm-two singleton caps, its pairing with old wall 26 is
exactly -2. Their sum is null and future-directed. This is an additional cusp;
the cap's original parallel partner is old wall 32 or 63 respectively.

The complete incidence calculation uses every signed family in the
5,397,548-facet lower-placement chamber. I checked the normalization constants
separately: the timelike family inequality is `u/12+sqrt110<0`, the old null
family is `d/6+f*sqrt11/2<0`, the k=6 family uses the numerator `r+792` over
36, and the new shifted families use `c+66*t*DEN/5` over `6*DEN`. Their negative
signs and exact squared comparisons cover both height signs. The center bound
uses height at most four and thus covers every old center label `m>=3`.
The code bounds native integer operations and verifies its frozen construction
dependencies by hashes.

Exactly seventeen original walls and the single new cap are incident at each
ray. Every other inherited or new wall is strictly negative. In the rational
coordinate convention with form `diag(Q,11)`, separate Fraction/SymPy local
calculations give normal rank 18 and Gram rank 17. The components are finite
`A5+A4+2A3+A1` and one affine A1. An explicit four-dimensional positive space
is orthogonal to all incident normals and to the center. Since the null ray is
the sum of the affine pair, these vectors are also tangent to the horosphere.
They give four independent lineality directions in its incident-wall link.
The positive leading minors of their Gram matrix are checked exactly.

These are two fully diagnosed unbounded links, not a classification of all
cusps, nor a claim that no further walls could repair them. They supplement
the already established strict global escape. The fresh replay uses the same
large-array implementation; this review is not a second independent enumeration
of all distant pairings.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h22_lower_new_cusp_20261013.py --verify
```

Output: `PASS_EXACT_TWO_NEW_UNBOUNDED_H22_CUSP_LINKS`, with incident counts
18, quotient ranks 17 and lineality dimensions 4 for both representatives.
