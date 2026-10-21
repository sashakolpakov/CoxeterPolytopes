# Root review: farthest-vertex energy

20 October 2026. Complete independent written review of
`compact_farthest_vertex_energy_obstruction_20261020.md`, SHA256
`50c64e5a521d3e4e0f277a8ece02d3a29aa4ae8b71485cf3088211e587f45bff`.
PASS. No execution, uniform ceiling or formal verification.

I independently obtained the edge formula and maximal-height inequality.
Writing z=E*p+z_T gives x=-z_T/E with norm squared1-E^(-2). Maximality
over ALL actual vertices gives tau_i<=x.u_i simultaneously for every
incident actual edge. Thus every nonnegative local weighting is valid
at this SAME vertex, including weights chosen after the local data.
The dual directions form a basis, making the Cauchy--Schwarz inequality
strict for any nonzero weighting.

The minimum-energy feasible set is nonempty, closed and convex, and its
squared norm has a unique minimizer. Its energy is strictly below one
at the chosen vertex. When K inverse times tau is nonnegative, the
equality vector and its nonnegative dual weights give matching bounds.
The source correctly forbids using that unconstrained expression as
the energy when a coefficient is negative. The rank-two boundary
branches and the sum-based bound check directly. Orthogonal component
energies separate and add exactly.

For A1^n the energy is the sum of half-tangent squares. The independently
reviewed two-short-direction bound gives the stated n>=12 contradiction;
49^2-6*20^2=1 checks the delicate strict endpoint without decimals.

This is a sufficient pointwise route to a global exclusion only after
EVERY possible farthest vertex type is controlled. It removes the
shared-white-weight restriction from that route, while leaving the
earlier undirected identity valid within its own scope. No all-vertex
classification, dimension-11 bound or new example has been obtained.
