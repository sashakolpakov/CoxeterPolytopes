# Independent review of the affine-A quadrilateral argument

**PASS within the stated all-simple, n >= 14 scope.** The
[proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md) pays an actual
quadrilateral with an affine-A corner. It does not settle the other affine
corner types or lower the ambient dimension bound.

Root separately read the proof and producer source. The
[independent checker](../scripts/review_simple_cusp_affine_a_quadrilateral_20261017.py)
uses Fraction pairs for Q(sqrt(2)), imports no producer, and recomputes all
256 height-rectangle values and the four positive kernel equations.
Its [certificate](review_simple_cusp_affine_a_quadrilateral_20261017.json)
binds the reviewed source, data and report. Generation and fresh saved
replay pass. The formal determinant expansion and positive-definite
deletion checks remain the producer's checks rather than independently
duplicated calculations.

The following geometric and analytic points were checked separately:

* Deleting one marked root from the full n-cycle leaves A_(n-1).
  An adjacent actual corner is finite or is a connected n-node affine
  diagram, because every vertex is simple. In rank n >= 14 the exceptional
  types cannot occur. The classical finite/affine diagrams give exactly
  the isolated, endpoint3, endpoint4, next-to-end3 and both-end3 extension
  profiles. Neither a hidden large label nor a disconnected affine
  extension is allowed by this classification.
* When both marked arcs have length at least four, the possible contacts
  lie away from the other old marked root. All stated low-contact
  exclusions are literal induced affine diagrams of order at most nine.
  In particular, a simple leaf on the old cycle supplies the affine-E8
  arms of lengths 1,2,5 without an extra edge between their ends.
  The isolated zero-contact case would be orthogonal to the full cusp
  span and hence null, contradicting the new wall's positive norm.
* The displayed path inverse gives the residual marked Gram, endpoint
  moments and cross term. At arc length four the two next-to-end contacts
  may coincide; the inverse still gives the asserted cross term 4/P.
  Thus this limiting overlap is included. The actual n+2-wall star
  necessarily has zero determinant in ambient Lorentz dimension n+1.
* Each height bound follows from the displayed forbidden small affine
  diagram and the Coxeter contact alphabet. Replacing a strict isolated
  bound by its closed endpoint enlarges the tested rectangle, which is
  safe. The residual polynomial is jointly concave in the reciprocal
  heights because its quadratic part is a negative diagonal form plus a
  negative square. Its minimum on a rectangle therefore occurs at a
  corner.
* Both reciprocal-arc partial derivatives are strictly negative at
  positive reciprocal heights. If a moment vanishes, either the other
  moment or the positive isolated constant supplies strictness. Since
  P+Q >= 14, at least one reciprocal arc is strictly below 1/4. The
  nonnegative rectangle calculation therefore becomes a strict positive
  determinant obstruction at every actual pair of arc lengths.

The inherited global affine exclusion is an explicit input. This proof
uses it only for affine diagrams of orders 3,4 and 9, all below n. It
does not assume that every affine subdiagram of any order is forbidden.
That distinction must be preserved when extending the argument to other
corner types or to quadrilaterals whose corners are all finite.
