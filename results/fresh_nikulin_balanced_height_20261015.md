# Balanced heights do not restore a local unit charge

The proposed local inequality

\[
W+C\ge 2r/9,\qquad
W=\min(\sigma,\tau)^2 E,
\]

fails even at **exactly balanced heights**, \(M=\max(\sigma/\tau,\tau/\sigma)=1\). It fails for every preset pair of omitted nodes in affine E8, not only the pair used in the frozen [4/3-charge proof](fresh_nikulin_unit_charge_20261015.md). Thus restricting the height ratio to \(1\le M\le21/20\), or to any smaller interval containing one, does not establish that local unit-charge lemma.

The [standalone Fraction checker](../scripts/fresh_nikulin_balanced_height_20261015.py) and [exact certificate](fresh_nikulin_balanced_height_20261015.json) have passed generation and saved replay. They import no project arithmetic. No frozen earlier files were changed. Independent review is separate from these producer checks.

## Explicit counterexamples

Use the affine E8 star with arm lengths `(1,2,5)`, numbered from its center along the arms as in the frozen proof. Its Gram matrix has diagonal two and edge entries minus one. Its positive kernel, normalized to minimum one, is

\[
\delta=(6,3,4,2,5,4,3,2,1).
\]

For the current omissions \((i,j)=(1,6)\), choose identical contact rows

\[
p=q=e_8.
\]

Both heights equal one. Their normalized difference is zero, so \(E=W=0\). Exact inverses of the two retained finite Grams give

\[
C_p=(G_{\widehat1}^{-1})_{88}=8/9,
\qquad
C_q=(G_{\widehat6}^{-1})_{88}=2/3.
\]

Consequently

\[
W+C=14/9<16/9.
\]

Every nonzero contact is exactly one, hence is an allowed norm-two \(\pi/3\) contact. Each endpoint cost is strictly below two. The corresponding endpoint Schur residuals are \(10/9\) and \(4/3\), so even the individual finite-endpoint positive-definiteness conditions hold. This is not an obstruction arising from non-Coxeter contact magnitudes or unequal heights.

The current ordinary E7 choice also fails. For the affine star `(1,3,3)`, omissions `(3,6)` and rows \(p=q=e_4\), the heights are one, the energy is zero, and the costs are \(1/2\) and one. Thus

\[
W+C=3/2<14/9.
\]

These two-contact configurations were safely outside the small-cost region \(C<3T/4\) in the 4/3 proof. They enter when a proposed unit-charge proof expands that region to \(C<T\). Excluding the earlier unequal-height one-contact examples therefore misses them.

## Complete E8 omission-pair census

The checker reconstructs all nine positive-definite principal deletions, their exact inverses, and all inverse-entry signs. For each of the 36 distinct unordered omission pairs it tests all nine common unit rows \(p=q=e_k\). Of the 324 rows, 179 obey both endpoint cost bounds. Every omission pair has at least one such row with

\[
0<C_p+C_q<16/9.
\]

Each serialized counterexample has \(M=1\) and \(W=0\). The result covers every preset distinct pair; for identical omissions \(i=j\), the row \(p=q=e_i\) gives \(C=W=0\) and is an immediate counterexample.

For each distinct pair the checker selects the least-cost legal common row. The minimum, over all pairs, of the resulting necessary coefficient

\[
\frac{16/9}{C_p+C_q}
\]

is **\(8/7\)**, attained only by omissions `(1,6)` among this census. Therefore every preset pair requires a coefficient at least \(8/7\), already at balanced heights. This is a lower bound on any proposed uniform local charge; it is **not** a proof that coefficient \(8/7\) suffices for arbitrary rows.

## Scope and remaining possibility

The examples satisfy the local hypotheses used in the energy accounting: nonnegative nonzero contact rows, the gap between zero and one, positive equal affine heights, and individual endpoint costs at most two. They refute a uniform local lemma on that domain. They do not assert simultaneous realization of all factors and both terminating walls in a finite-volume polytope.

The result leaves open global estimates using surplus from other factors, restrictions coupling the endpoint rows across factors, or omission choices that depend on the eventual rows. It also leaves open proving a smaller balanced coefficient above one. No such positive result is established here.

For scale only, a hypothetical balanced \(8/7\) upper lemma, combined with the existing 4/3 bound outside \(M\le21/20\), would yield the accounting constant \(18(8/7)+19=277/7\). The first integer \(N=n-1\) for which \((N-277/7)/9\ge14\) is 166, suggesting a conditional target of 166. This arithmetic is recorded solely to distinguish a possible future target from a proved dimension bound. The frozen reviewed 169 result is unchanged.

```sh
python scripts/fresh_nikulin_balanced_height_20261015.py
python scripts/fresh_nikulin_balanced_height_20261015.py --verify
```
