# A finite parameter reduction for compact acute quadrilaterals

For a compact acute hyperbolic quadrilateral face, normalize the four framing roots to square two and denote their opposite full products by u and v. Then

`u > −6 or v > −6`.

Over Q(sqrt5), arithmetic integrality and conjugate positivity turn the smaller framing product into a finite exact choice. In the current ADE-corner applications there are only8 choices if that product belongs to the ground field, or59 choices for its square otherwise. The other product must satisfy a nonzero quadratic equation. This is a finite reduction for a fixed quadrilateral diagram, not an independent bound on both products and not a dimension bound.

## The physical bound

Eliminate the common positive-definite white support W. Its Gram is a positive-definite matrix with nonpositive offdiagonal entries, so its inverse is entrywise nonnegative. The four-root residual S has diagonal entries `0<a_i<=2`. Every adjacent two-root corner is positive definite, giving

`−sqrt(a_i*a_j) < S_ij <= 0`, hence `−2 < S_ij <= 0`.

The white projection products q and r at the opposite pairs are nonnegative, so `S_02=u−q<=u` and `S_13=v−r<=v`. Restrict S to the two vectors `e0+e2` and `e1+e3`. If both u and v were at most minus six, the two diagonal entries would be at most minus eight. Their offdiagonal entry is the sum of the four adjacent entries, whose magnitude is strictly less than eight. This restriction would be negative definite, contradicting the ambient Lorentz index one. The strict bound follows.

## Complete arithmetic choices for the smaller product

Write phi=(1+sqrt5)/2. A squared normalized Gram product lies in the integer ring Z[phi]. At the nonphysical embedding, positivity gives `0 <= (g²)' <4`. If g itself belongs to the ground field, it is also integral and `−2<g'<2`.

The complete coefficient boxes are elementary. For `g=a+b*phi` with `−6<g<=0`, subtracting its conjugate gives `−4<b<1`, hence b is one of−3,−2,−1,0. Then `−6<a<6`. For `z=g²=a+b*phi` with `0<=z<36`, the analogous bounds are `−2<b<18` and `−34<a<38`. The checker enumerates safe supersets of these integer boxes and filters each point by exact inequalities.

When `−2<g<=0`, the finite Coxeter-label possibilities over Q(sqrt5) are2,3,4,5,6,10. The parallel value g=−2 is forbidden by compactness. Combining those possibilities with the dotted range leaves10 direct-product choices or63 squared-product choices before using the ADE-edge premise.

In the applications here, all other contacts of the relevant edge Gram are zero or minus one. A rational framing product would make the full edge Gram rational and prevent the required physical/conjugate signature change. If its square is rational but the product is not in the ground field, a fixed label-three path violates the cycle-field condition; without such a path its edge Schur entry is rational and again has the same sign at both embeddings. This removes the two rational direct choices and the four rational square choices, leaving8 and59 respectively. The certificate retains both the removed entries and all surviving entries. The8/59 restriction must not be transferred without this fixed-contact premise.

## Why an identically zero specialization is excluded

The coefficient of v² in `det(S)` is minus the determinant of the opposite residual block on roots0 and2. Therefore a specialization that makes the determinant identically zero as a polynomial in v makes that block singular.

Its diagonal entries are positive and its offdiagonal entry is nonpositive. Singularity thus gives a strictly positive kernel vector y. Lift this vector to the original white-plus-black Gram by adding the white coefficients `−W_inverse*B*y`. These are nonnegative because the white inverse is nonnegative and the contact matrix B is nonpositive.

The resulting nonnegative combination of wall normals is null. It is nonzero because its product with a common strict timelike center is negative, and it lies on the future sheet. Its product with every wall normal is nonpositive by acuteness. It therefore gives a feasible ideal ray, contradicting compactness. Both the nonzero and feasibility steps are necessary; a formal zero Gram determinant alone would not establish an ideal point.

Consequently a fixed bounded anchor leaves a nonzero polynomial of degree at most two for the other framing product. Both choices of which member is bounded must be tested. Physical signature, conjugate positivity, integrality, finite angles and face incidence must still be checked on the resulting roots. This argument does not bound a lone edge's unit-translation families and does not apply unchanged to noncompact faces.

The [source](../scripts/fresh_sqrt5_compact_quad_finite_anchors_20261017.py) and [exact domain](fresh_sqrt5_compact_quad_finite_anchors_20261017.json) pass generation and a fresh saved replay. The source independently verifies the generic determinant coefficient identity and enumerates the domains with the frozen Fraction field engine:

```sh
python scripts/fresh_sqrt5_compact_quad_finite_anchors_20261017.py
python scripts/fresh_sqrt5_compact_quad_finite_anchors_20261017.py --verify
```

The two geometric arguments received a separate mathematical read during development. The forthcoming diagram-by-diagram use of these domains is not part of this package.
