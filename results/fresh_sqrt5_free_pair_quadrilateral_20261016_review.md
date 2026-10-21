# Independent review of the fixed-edge free-pair quadrilateral obstruction

**PASS.** The continuous exclusion and its stated local scope survive an independent reconstruction. The [review checker](../scripts/fresh_sqrt5_free_pair_quadrilateral_20261016_review.py) imports no project or producer arithmetic. It builds the literal rational `D5+E8+E8` Gram, inverts it by Fraction LDL and triangular solves, and uses the basis `1,√5` for exact signs. Its [certificate](fresh_sqrt5_free_pair_quadrilateral_20261016_review.json) binds the frozen [primary source and data](fresh_sqrt5_free_pair_quadrilateral_20261016.json).

Generation and a fresh saved replay pass all **210 pairs, 420 determinant endpoints, and 420 strict endpoint-sector inequalities**. The split is144 different-component pairs and66 same-component pairs. Both physical and conjugate signs of all21 fixed-edge residual diagonals are checked. The review reconstructs the specific `T=2φ⁸` edges; it does not repeat the predecessor's universal census of other vertex types or edge parameters.

## Why the interval test proves the claimed obstruction

Write `C=A⁻¹`, where `A` is the old positive-definite rank21 Gram. Fixing both new roots' contacts with all old roots fixes their spatial projections to `s_i=e_i−Tω_i` and `s_j=e_j−Tω_j`. Their remaining diagonal entries are

`p_i=2T−T²Cii<0`, `p_j=2T−T²Cjj<0`.

The exact spatial mutual product is `b=Aij+T²Cij`. Hence the full23-root Gram is congruent to `A` plus

`[[p_i,x−b],[x−b,p_j]]`.

A fourth finite corner must contain both new roots. Its two-root positive-definiteness and acuteness require `−2<x≤0`. The primary relaxes this to `−2≤x≤0`, which strengthens the exclusion rather than discarding any allowed boundary.

For each pair, the independent arithmetic confirms `D(−2)>0` and `D(0)>0`, where `D(x)=p_i p_j−(x−b)²`. The identity

`D(−2+2t)=(1−t)D(−2)+tD(0)+4t(1−t)`

proves positivity throughout that interval. Since both diagonals are negative, this block is negative definite. The full matrix consequently has inertia `(21,2)` and rank23, incompatible with the22-dimensional Lorentz space for H21. This needs no integrality, conjugate condition or finite-angle alphabet.

The smallest endpoint determinant is exactly `20169647+32635176φ`, at ports2 and3 with `x=−2`, as reported by the primary. The inequalities `p_i p_j>(b+2)²` and `p_i p_j>b²` independently place the two rank-closure values strictly outside the allowed interval: one below−2, the other above0.

## The real neighboring sectors are present

With a unit future timelike vector `t₀` perpendicular to the old roots, let `a_i=√(−p_i)` and set

`new_i=e_i−Tω_i+a_i t₀`,

`v_i=t₀−a_i ω_i/(TCii−1)`.

The denominator is positive. Direct substitution gives square `−1/(TCii−1)²`, zero pairing with `new_i` and all unreleased old roots, and strictly negative pairing with the released old root. Its active Gram is literally `A` with one root replaced by its copy. The coefficient of `t₀` is1, so the endpoint is future timelike.

For the other new wall the only potentially positive term is `a_i T Cij/(TCii−1)`. Entrywise nonnegativity of `C`, positivity of the denominator, and the independently checked inequality

`(−p_j)(TCii−1)²−(−p_i)T²Cij²>0`

show that this term is strictly smaller than `a_j`. This proves the other-new-wall inequality strictly, in both directions of every pair. The original vertex and both neighboring vertices therefore lie in one actual real partial chamber, and the two edge segments have the stated active walls. A sufficiently small inward perturbation at the original spherical corner supplies a common strict center.

The negative mutual-product branch is ultraparallel. A face retaining these two consecutive finite edges cannot be a quadrilateral, but further side walls could produce a longer face. This review makes no claim of arithmeticity or a common integral lattice for the declared square-root realization, no completeness claim for an ambient edge star, and no compactness, finite-volume or dimension21 exclusion.

```sh
python scripts/fresh_sqrt5_free_pair_quadrilateral_20261016_review.py --verify
```

All recursive input source hashes are checked. The primary Markdown is source-read but intentionally not hash-bound, allowing a review link to be recorded separately without changing the reviewed mathematics.
