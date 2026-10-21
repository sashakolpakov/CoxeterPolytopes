# Independent review of the compact quadrilateral anchor reduction

**PASS.** The [finite-anchor reduction](fresh_sqrt5_compact_quad_finite_anchors_20261017.md)
has the stated geometric scope and exact arithmetic domains. The independent
[checker](../scripts/fresh_sqrt5_compact_quad_finite_anchors_review_20261017.py)
imports no producer or project arithmetic. It works in the integer coordinates
`(p+q*sqrt5)/2`, with p and q of the same parity, and decides every sign using
integer squares. Its complete rectangular lattice enumeration reproduces all
eight direct products,59 squared products, and the separately retained two
and four rational exclusions. Generation and a fresh replay match the
[saved certificate](fresh_sqrt5_compact_quad_finite_anchors_review_20261017.json).

The physical proof is valid. Positive definite adjacent corner blocks bound
each of the four adjacent residual products strictly above minus two. White
elimination decreases the opposite products and leaves diagonal entries at
most two. If both full opposite products were at most minus six, the sums of
opposite residual normals would span a negative definite plane. This is
impossible in Lorentz signature. The conclusion is that at least one product
is strictly greater than minus six, not that both products are bounded.

The identically-zero specialization argument is also valid. The coefficient
of the other squared product is minus the first opposite block's determinant.
If this vanishes, the positive kernel of that block lifts to a nonnegative
combination of actual wall normals. The common strict center proves this
combination is nonzero and future null; acuteness puts it in every halfspace.
It would be an ideal point, contradicting compactness. Thus the remaining
quadratic cannot be identically zero at a realizable bounded anchor.

The arithmetic assumptions are necessary. A normalized squared product is an
integral cycle invariant. When the product itself is in the ground field, it
is integral because its square is integral and the integer ring is integrally
closed. Equality in the conjugate bound would force physical squared product
four, which compactness forbids. Removing the rational possibilities uses
the specific ADE edge with all other products zero or minus one; the smaller
eight/59 lists are not universal lists for arbitrary finite corner diagrams.

This review establishes the finite reduction for each fixed compact
quadrilateral diagram over Q(sqrt5), under the stated fixed-contact premise.
It does not solve the resulting quadratics or establish a dimension bound.
It is an internal mathematical and computational review, not external
refereeing or formal verification.

```sh
python3 scripts/fresh_sqrt5_compact_quad_finite_anchors_review_20261017.py --verify
```
