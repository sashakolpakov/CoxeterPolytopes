# Two norm classes suffice for the arithmetic root search

For the [explicit free lattice](fresh_h9_h8_cap_height_20261014_host_basis.md), every primitive lattice-preserving reflection direction has a representative of norm **2 or `2(3+φ)`**, after rescaling its normal by a unit. This is a complete norm alphabet, not a bounded sample. It does not enumerate the roots in either class.

Write `O=Z[φ]`, `β=3+φ`. The lattice is even, and its Gram determinant is `−φβ`. The element `β` is totally positive and has prime field norm eleven. The [certificate](fresh_h9_h8_cap_height_20261014_norm_pool.json) checks these exact algebraic inputs.

First, `O` is norm-Euclidean. Rounding the two coefficients of any field element leaves an error `a+bφ` with `|a|,|b|≤1/2`, whose norm satisfies `|a²+ab−b²|≤3/4<1`. Hence `O` is a principal ideal domain and every rational reflection line has a primitive integral representative.

Let `r` be such a primitive vector, with positive norm `q`. Evenness gives `q=2b`, with `b` in `O`. Reflection integrality first gives `r_i(Gr)_j/b` in `O` for every pair of coordinates. A Bezout combination of the primitive coordinates `r_i` then shows that `G r` belongs to `b O¹⁰`. Multiplication by the adjugate gives `(det G) r` in `b O¹⁰`; a second Bezout combination shows that `b` divides `det G`. The only possibilities are therefore `b=u` or `b=uβ`, with `u` a unit.

The norm `q` is totally positive: it is positive at the hyperbolic embedding, and the conjugate form is positive definite. Since `β` is also totally positive, `u` is a totally positive unit.

For completeness, the positive units are powers of `φ`. There is no unit `v` strictly between 1 and `φ`: if `Norm(v)=1`, its integer coefficient `(v−1/v)/√5` lies strictly between 0 and 1; if `Norm(v)=−1`, the coefficient `(v+1/v)/√5` also lies strictly between 0 and 1. Both are impossible. Multiplying any positive unit by a suitable power of `φ` puts it in `[1,φ)`, proving the claim. Since `Norm(φ)=−1`, the totally positive units are exactly `φ^(2k)`. Rescaling `r` by `φ^(−k)` gives norm 2 or `2β` as asserted.

For enumeration, every integral norm-two vector gives an integral reflection. A norm-`2β` vector must additionally satisfy `Gr∈βO¹⁰`; the norm equation alone does not suffice.

The [checker](../scripts/fresh_h9_h8_cap_height_20261014_norm_pool.py) verifies the finite algebraic inputs and stores the two norm classes. The Euclidean, adjugate and unit arguments above supply the universal deduction. The scope is reflection directions preserving this lattice; a general nonarithmetic geometric completion is not restricted to this alphabet. Reflectivity and any finite-volume completion remain open.

```sh
python scripts/fresh_h9_h8_cap_height_20261014_norm_pool.py --verify
```
