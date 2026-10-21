# Independent written review: projection and trace-lift obstructions

20 October 2026. Full source read; paper algebra only. No computation,
enumeration, replay, Lean, or build was performed.

Reviewed source:
`results/first_type_projection_and_trace_lift_obstructions_20261020.md`,
SHA256 `91566bd43e57a8edeef496e8da0e8704d677de96e15d38103522834d65848234`.

**PASS within the stated scope.** Both examples are actual compact arithmetic
chambers before the proposed transfer, and the second lifted system is an
actual integral Coxeter wall system with infinite volume. Neither example
asserts a new dimension bound or a compact example above eight.

I independently checked the intrinsic-normal formula
`mbar=(m+c*n)/sqrt(1-c^2)` and the projected height
`(t_m+c*t)/(sqrt(1+t^2)*sqrt(1-c^2))`. For the compact triangle, the
norm-two Gram determinant is `-2sqrt2`, its conjugate is positive, and the
displayed inverse-Gram row sums give `t^2=1/(8+9sqrt2)`. At the facet with
two3-angles the new squared height is `(sqrt2-1)/3`; its ratio to the old
squared height is `(10-sqrt2)/3>1`. Thus a first-type assumption cannot
restore the right-angled contraction without another hypothesis.

For the quadrilateral I rederived `a3=sqrt2*(a2+a4)-a1`, all four adjacent
products `-1`, the two opposite products, and the strict timelike vector
`v=a1+a3`, with square `-4sqrt2`. The conjugate basis Gram is positive
definite. The four angles are pi/3 and all four corners are finite; this is
a compact chamber, rather than merely a local Gram model. Integrality and
the positive conjugate justify the stated first-type arithmetic scope.

In the trace lattice I checked both reflection coefficients and all types
of normalized products. In particular the cross products for opposite pairs
are `-sqrt2` and `-1/sqrt2`, respectively; the corresponding within-layer
opposite products are `-1`, and the cyclic cross products vanish. The real
point `(q,0)` is strict, so all eight walls are actual by acute projection.

The null vector `z=a1+a3=b2+b4` has exactly the eight displayed signs
`(0,-4,0,-4;-8,0,-8,0)`. Its four incident walls are two orthogonal affine
A1 pairs. Their positive normal sums give the same null direction, and
their horospherical normals span only two dimensions in the four-dimensional
horosphere. Consequently a high cusp section is a rectangle times R2. The
other walls have strictly negative cusp product and disappear above a
sufficiently high horosphere. This proves infinite volume; the argument does
not silently add other mirrors of the full trace lattice.

The inherited inputs are the usual hyperbolic Gram realization/arithmeticity
criterion, the discrete integral orthogonal-group construction, and the
standard horospherical description of a cusp. The source does not claim that
the complete reflection group of the trace lattice is nonreflective. I found
no missing global volume premise in the narrower counterexample it proves.
