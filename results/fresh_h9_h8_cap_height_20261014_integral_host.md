# An integral arithmetic host for the maximal-height seed

The maximal-height chamber gives an explicit crystallographic reflection module over `O=Z[φ]`, where `φ=(1+√5)/2`. This supplies a concrete host for an arithmetic root search, while leaving reflectivity open.

Scale the bottom normal by two and leave the other fourteen normals unchanged. The bottom then has norm `7+√5=6+2φ`; every other normal has norm two. All entries of their Gram matrix and all reflection coefficients

`C_ij=2(r_i,r_j)/(r_j,r_j)`

belong to `O`. In particular, the bottom/cap coefficients are `−(3+φ)` and `−1`. Thus the `O`-span `L` of these fifteen normals is preserved by every generating reflection, since `s_j(r_i)=r_i−C_ij r_j`. It has field rank ten. Scaling the bottom normal changes no hyperplane.

The [exact certificate](fresh_h9_h8_cap_height_20261014_integral_host.json) stores all 225 Gram entries and all 225 reflection coefficients as integer pairs in the basis `1,φ`. The [checker](../scripts/fresh_h9_h8_cap_height_20261014_integral_host.py) verifies them with Fraction arithmetic and separately with exact symbolic matrix products. It also checks full rank and signatures `(9,1)` and `(10,0)` at the two real embeddings.

This distinguishes the maximal-height host from the earlier small-cap form, whose nontrivial conjugate was Lorentzian. It does not show that these fifteen walls constitute the full reflection chamber, that `L` is reflective, or that the reflection subgroup has finite index. The chamber still has its certified strict escape. An arithmetic continuation must enumerate additional roots and prove global termination or a valid obstruction; no unimodular classification is applied here.

```sh
python scripts/fresh_h9_h8_cap_height_20261014_integral_host.py --verify
```
