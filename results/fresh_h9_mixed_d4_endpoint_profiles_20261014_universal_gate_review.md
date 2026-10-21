# Review of the universal mixed-D4 finite-endpoint gate

The frozen [source](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py) and [certificate](fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.json) pass this mathematical/source review. Their conclusion covers all1,825 finite profiles and72 continuous dihedral tails on the52 saved old finite-corner edges. Every valid cap in these families pairs strictly negatively with the saved rational spacelike point. This is not a global obstruction to other caps, affine endpoints, newly created edges or changed cores.

**Orientation and completeness.** For a genuine added facet acute to the old complete finite corner, the forward coefficient is positive. Normalize its future corner vector U to square−1 and let ω_i be the nine incident dual roots. An acute normal has the form `r=aU−Σc_iω_i`, with c_i≥0. Any strict future point of the old corner sector has the form `z=bU−Σd_iω_i`, with b>0 and d_i>0. The incident positive definite acute Gram matrix has entrywise nonnegative inverse J. If a≤0, then `r·z=−ab+cᵀJd≥0`; such a normal cannot define a facet of a chamber with strict interior in that sector. Hence a>0. The positive residual branch meets the outgoing half-edge, and the opposite branch meets the backward half-edge. This argument allows the common interior point to change after adding caps; it does not assume the original center remains interior.

Thus the finite profiles have the necessary form

`r=p+hU+sqrt(h²+rho)V`, with h>0 and rho>0.

Including h=0 is a safe relaxation. For tails, the closed interval `sqrt3≤k≤2`, with `p=kp0` and `rho=2−k²/2`, is also safe. At k=2 it retains the nonzero null residual `h(U+V)`; it does not lose a limit by dividing by sqrt(rho). This boundary is included only to enlarge the finite-endpoint parameter set, not to claim completeness for arbitrary affine endpoints.

**Unbounded parameters and interval endpoints.** For a scalar pairing `b+ah+c sqrt(h²+rho)`, the certified conditions a≤0 and a+c<0 imply monotonic decrease in h≥0 for either sign of c. At rho=0 the formula is simply linear on h≥0, including its right-hand value at zero. For the escape pairing the stronger signs a<0<c yield the linear upper bound `b+c sqrt(rho)+(a+c)h`. The recorded integer H therefore bounds every possible escape cutter, uniformly over the whole tail interval. This also controls the joint limit k→2, h→infinity.

Each parameter box is closed. The two children overlap at their splitting point, and rejection requires either a strictly negative escape upper bound or strict disjointness from the allowed pairing relaxation. Singleton boundaries are retained. The relaxation keeps the isolated values for labels2 through15 and the entire closed ray below the label16 value; it includes every larger finite label, parallel contact and dotted contact. There is no maximum-label assumption.

**Independent and shared checks.** The [separate review checker](../scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate_review.py) imports no producer arithmetic. It checks the exact one-to-one profile/tail coverage, all source manifests, valid outside-wall rejection indices and the binary-tree count identities for all9,941 boxes. It independently encloses π using `4atan(1/2)+4atan(1/3)`, then uses rational endpoint Taylor sums to obtain tighter bounds for every cosine used by the primary. These bounds lie inside the primary's Machin/Taylor enclosures. Generation and saved review replay pass.

The primary's degree128 cosine polynomial has its next omitted term at degree130; `2^130/130!` is a valid error bound because every argument is below2. Its dyadic arithmetic rounds all endpoints outward. The fresh full primary saved replay also passes, using that same interval engine. This review does **not** claim a second implementation of all quartic coefficients or all individual box decisions.

The full source proof, including the finite-union infinite-volume consequence, has no identified gap under its stated scope. A finite coherent collection from these families leaves the spacelike point strict; together with a common strict future interior point this gives an open ideal-boundary patch. Later endpoint supports or unrelated facets are outside this theorem.

```sh
python scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py --verify
python scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate_review.py --verify
```
