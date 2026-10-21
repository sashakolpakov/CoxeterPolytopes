# Independent arithmetic for the product-cusp bound

The [standalone checker](../scripts/fresh_nikulin_product_independent_20261015.py) and [certificate](fresh_nikulin_product_independent_20261015.json) independently verify the affine constants and scalar estimates in the [conditional bound of 225](fresh_nikulin_narrow_20261015_products.md). They import no producer arithmetic. Generation and saved replay pass.

The exceptional \(\widetilde A_1,\widetilde E_6,\widetilde E_7,\widetilde E_8,\widetilde F_4,\widetilde G_2\) energies are reconstructed by inverting a proper finite subdiagram and pinning a null coefficient. This differs from the producer's modified pseudoinverse. Their weighted energies are respectively \(1/2,4/3,2,2,1,2/3\).

Independent Euclidean coordinates check all standard and modified deletion choices for the classical B, C and D factors at ranks 10, 12 and 17. Contact-simplex extreme points verify the claimed constant-coordinate intervals; these finite checks supplement the all-rank coordinate proof, rather than replace it. Symbolic parity and rank tails verify the scalar inequalities for every larger rank. At total cusp rank \(N=225+t\), the final lower bound is exactly \(14+t/9\).

The geometric common-height identity, endpoint wall selection and narrow-facet hypothesis require the accompanying mathematical proof and its hostile review. This checker alone is not an unconditional dimension theorem.

```sh
python scripts/fresh_nikulin_product_independent_20261015.py --verify
```
