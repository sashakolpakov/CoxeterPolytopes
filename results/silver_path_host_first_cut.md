# Silver path host: exact first cutting packet

The rank-23 `Z[sqrt(2)]` lattice has Gram diagonal
`(2c,2,...,2)`, off diagonal `(-c,-1,...,-1)`, and
`c=2+sqrt(2)`.  It is not an orthogonal sum.  Exact invariants are

- determinant `-40-42sqrt(2)` (field norm `-1928=-8*241`),
- selected inertia `(22,1,0)`, conjugate inertia `(23,0,0)`,
- Weyl square `(-40-42sqrt(2))/23`, and minimal field `Q(sqrt(2))`.

Since the conjugate form is positive definite, a nonzero isotropic vector
over `Q(sqrt(2))` is impossible (conjugate its coordinates).  Thus the host
is anisotropic; if its reflection chamber had finite volume, it would be
compact.

The determinant valuation bound gives twelve totally-positive mirror-norm
classes modulo unit squares.  Complete exact enumeration through

`H=(144+103sqrt(2))/23 ~= 12.5940868228`

has 361 controls and 7 integral candidates.  Chronological Vinberg selection
accepts the endpoint root `e0`, followed by two simultaneous norm-two roots:

```
(1+r, 2+2r, r, -1, ..., -1)
(1+r, 3+2r, 2+r, 1, -r, -1-r, ..., -1-r), r=sqrt(2).
```

They are acute to all 23 simplex walls and have mutual pairing `-2-r`.
Thus the first cutting cone has 25 facets.  Independent exact enumeration in
the simplex dual coordinates gives 103 extreme rays: 6 finite vertices,
0 cusps, and 97 hyperideal rays.  Hence this exact packet is not compact and
not finite volume; this is a bounded result, not a global nonreflectivity
claim.

A separate small-coefficient topology scan covered every admissible integral
`a+b sqrt(2)` with `|a|,|b|<=50` and every attachment node of `A22`, `D22`,
and `B22`.  The only best initial cones are the
four end-attached `c=2+sqrt(2)` cases, each with 3 finite and 20 hyperideal
simplex vertices.  All other 994 admissible cases have 2 finite and 21
hyperideal vertices.  Therefore switching among these tree skeletons does
not improve the exact initial cone.
