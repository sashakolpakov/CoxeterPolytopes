# Cubic-seven full-rank path host in H22: exact seed cone

Let `theta=2*cos(2*pi/7)` and `c=2+theta`.  The rank-23 lattice has
symmetric Gram matrix with diagonal `(2c,2,...,2)`, first off-diagonal
`-c`, and the remaining path off-diagonals `-1`.  This is a genuinely
non-diagonal, non-orthogonal-sum arithmetic host.

## Exact arithmetic

- `N(c)=1`; the nonrational Cartan two-cycle product is `c`, so the minimal
  Vinberg field is exactly `Q(theta)`.
- `det(B)=4-42theta-22theta^2=-2c(11c-23)` and
  `N(det(B))=-12488=-2^3*7*223`.
- Exact inertias at the three real places are respectively `(22,1,0)`,
  `(23,0,0)`, `(23,0,0)`.  Thus the host is admissible.
- Positive definiteness at either conjugate place proves anisotropy over
  `Q(theta)`, hence the ambient arithmetic orthogonal lattice is cocompact.
  This does **not** prove its reflection subgroup has finite index.
- The Weyl vector orthogonal to the 22 spatial path roots is
  `rho=(1,c*22/23,c*21/23,...,c/23)` and
  `rho^2=det(B)/23`.

## Exact seed-cone geometry

The 23 displayed path mirrors form a simplicial cone.  For a simplex, the
ray opposite facet `i` has norm sign equal to `(B^-1)[i,i]`.  Exact inversion
gives:

- finite vertices: omitted facets `0,1,2` (3 total);
- cusps: none;
- hyperideal rays: omitted facets `3,...,22` (20 total).

Therefore this initial cone is neither compact nor finite-volume.  Later
Vinberg walls could cut these rays, so this is not a nonreflectivity proof.

## Small-coefficient optimality and scope

Among all integral `a+b theta+c theta^2` with coefficients in `[-10,10]`
that give selected Lorentzian signature and two positive-definite conjugate
forms, there are 14 choices.  This coefficient is the only one whose seed
simplex has three finite vertices; every other choice has two.  A separate
certified Minkowski-box enumeration shows it is the unique admissible
algebraic integer between the Lorentzian threshold `23/11` and `2+theta`,
so it is the closest admissible integral terminal coefficient to threshold.

No positive-height Vinberg completeness claim is frozen here.  In
particular, this record is a new exact host and seed-cone computation, not a
new compact Coxeter polytope.
