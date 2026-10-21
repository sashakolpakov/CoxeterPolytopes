# Independent review of the core shared-single family

**PASS.** No mathematical issue was found in the
[primary proof](fresh_h9_star_mixed_singles_20261014.md). The conclusion is
limited to star case 0 with first three finite internal labels at least 4,
last finite internal label at least 7, leaf 14 having own contacts 3/3, and
leaves 11 and 17 sharing one own contact 3 while both are right to the other
last-block root. All noncore contacts remain unrestricted. Both last-root
orientations are included. This is not an exclusion of the whole incidence.

The [independent checker](../scripts/fresh_h9_star_mixed_singles_20261014_review.py)
and [certificate](fresh_h9_star_mixed_singles_20261014_review.json) reuse only
the previously independent Fraction field with basis
\(1,\sqrt2,\sqrt5,\sqrt{10}\). They import no producer arithmetic or
symbolic algebra. The new fallback calculation uses Laplace determinant
expansion and Cramer's rule, rather than the producer's LDL recurrence.

The review independently checks the literal eight-root support, all core
meeting pairs, the nine maximal induced core faces, and every selected
containing vertex among the 378 saved vertices. Each of the 27 domain
guards is reconstructed from the stated own equalities, internal lower
bounds and one trial contact. No nonright matrix entry is assumed from an
unrestricted leaf. Congruence by its positive kernel produces a rational
Laplacian with zero row sums, directly certifying failure of positive
definiteness. This covers the 24 cross-block guards, the mutual-anchor
triangle and the two fresh affine-F4 upper guards at vertex 295.

All three eight-root budget supports are checked literally. Retaining one
own contact contributes \(1/4\). A foreign nonright contact contributes at
least \(1/2\), so at most one can occur and its label must be 3. A foreign
block of internal label at least 6 has inverse diagonal at least 4 and
cannot receive even one such contact. The actual triangle
\(\{12,13,14\}\) forces leaf 14 right to the last block. These arguments
use positive definiteness and nonnegative inverse entries; they do not
require noncore labels.

The high-first proof is sound. After projecting off the positive
eight-root support, the vectors \(r_{14}\) and \(r_{11}+r_{17}\) have norms
\(-(K+H)\), \(-(L+F)\), where \(K\ge1\), \(L>2\) and \(F,H\ge0\).
Their own-block cross projections vanish by the fresh budgets. Their
ambient product has absolute value at most 1, and Cauchy's inequality
bounds the remaining projection product by \(\sqrt{FH}\). The independent
formal coefficient check verifies

\[
 (1+H)(2+F)-(1+\sqrt{FH})^2
 =1+(\sqrt F-\sqrt H)^2+H>0.
\]

The residual Gram would be negative definite, contradicting Lorentz index
one. The proof covers arbitrarily large finite internal labels.

When the first internal label is 4 or 5, any high middle block has zero
core contacts, so changing its internal label to 4 leaves the residual
matrix unchanged. The review checks the complete Cartesian set of 7,840
inherited affine-case keys, their signs and boundary identities, and the
hash binding to the prior
[independent shared-high arithmetic review](fresh_h9_star_shared_high_20261014_review.md).
It does **not** recompute those determinant coefficients. Their numerical
formulas involve only the core; the earlier noncore domain guards are
replaced by the new guards above.

There are 7,392 positive affine cases. All remaining 448 cases now have
fresh, independently reconstructed witnesses supported entirely in the
core. They use 39 distinct selected full matrices; the producer searched
60 matrix types including unsuccessful trials. For each selected matrix,
the review checks all leading principal determinants through the first
failure and reconstructs a vector \(v\ge0\), \(v\ne0\), with
\(v^TGv\le0\). Every placement is contained in its saved actual vertex.
These explicit vectors also justify monotonicity: increasing the relaxed
internal labels or adding other acute contacts cannot remove the failure.

Finally, the swap of roots 12 and 13 is checked on the complete vertex set
and minimal-nonface collection, and both theorem predicates are verified
entry by entry. No global symmetry beyond this literal automorphism is
assumed.

Generation and saved replay pass. All recursive source manifests and local
report links pass. Frozen predecessor files were not changed.

```sh
python scripts/fresh_h9_star_mixed_singles_20261014_review.py --verify
```
