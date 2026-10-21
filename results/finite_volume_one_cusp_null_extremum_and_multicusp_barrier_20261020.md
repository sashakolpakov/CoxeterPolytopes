# Null-height selection at one cusp, and its precise obstruction at two

This is a paper lemma for finite-sided, finite-volume acute hyperbolic
polytopes. It supplies an extremal-vertex interface, not a new universal
dimension bound. In particular, it does not replace the cusp-shape or
payment arguments of `docs/finite_volume57_lemma_map.md`.

Use Lorentz signature `(n,1)` and future unit points of square `-1`.

## The one-cusp selection lemma

Suppose P has exactly one ideal vertex, represented by a future null
vector xi. Then

`E(q) = -(xi,q)`

attains a positive maximum at a finite vertex p of P. Write

`xi = E(p) (p-x)`, with `x in p-perp`.

Then `|x|=1`. For an outgoing edge at p use the unit dual direction u_i,
so its finite endpoint is

`p_i = cosh(L_i) p - sinh(L_i) u_i`.

Every finite edge satisfies

`(x,u_i) >= tau_i := tanh(L_i/2) > 0`.                 (1)

If the spherical vertex diagram at p has at least two whole orthogonal
components, **every edge at p is finite**. Consequently (1) holds in all
directions, and the sum of the minimum covering energies of the whole
components is at most 1.

For the maximum assertion, truncate the sole cusp sufficiently high.
The remaining polytope is compact, whereas E tends uniformly to zero
down that cusp: in upper half-space with xi at infinity it is a positive
constant divided by the vertical coordinate. Thus E has a positive
maximum in the finite part. Along each unit-speed geodesic `E''=E>0`.
A point in the relative interior of a positive-dimensional face cannot
be a maximum, so a maximizing point is a finite vertex.

Nullity gives `|x|=1`. Comparing E(p_i) with E(p) and using
`(cosh L-1)/sinh L = tanh(L/2)` proves (1).

If an edge i ends at xi, its null ray is `p-u_i`, so x=u_i. An edge j
from a different whole component has `(u_i,u_j)=0`. It cannot also end
at xi, since there is only one geodesic from p to a specified ideal
point. It is therefore finite, and (1) contradicts `(x,u_j)=0`.
Finally the whole-component tangent spaces are orthogonal, so their
minimum covering energies add and their sum is bounded by `|x|^2=1`.

This differs from an interior-point farthest-vertex argument: its budget
is **at most 1**, rather than strictly below 1. A local obstruction
giving energy strictly greater than 1 remains useful; an argument that
only gives energy at least 1 does not transfer automatically. Local
arguments that use global absence of affine subdiagrams also require
separate checking.

## The common-opposite obstruction survives

Assume the maximizing vertex has two whole components, and let i,j be
directions from different components. No actual outside wall c can be
the opposite wall on both finite edges.

Here is a direct proof allowing arbitrary contacts inside either
component. Let W be the old vertex walls other than i,j, and write their
positive-definite Gram as B_W. Add to c its nonnegative orthogonalizing
combination of W; inverse positivity gives coefficients
`alpha = B_W^{-1} beta_W >= 0`. Call the resulting vector c_F. Since
xi lies in every actual closed half-space, `(xi,c_F)<=0`. If c closes
both edges, orthogonality of u_i,u_j gives

`c_F = A p - (A/r_i)u_i - (A/r_j)u_j`,

where `A=-(p,c)>0` and `r_i=tanh L_i`, `r_j=tanh L_j`.
Thus `(x,u_i)/r_i + (x,u_j)/r_j <= 1`. But (1) gives the strictly
larger lower bound

`tau_i/r_i + tau_j/r_j = 1+(tau_i^2+tau_j^2)/2 > 1`.

No strict-interior sign, compactness, arithmeticity or restriction to
one-contact wall profiles was used here.

## Why a causal linear height cannot select such a vertex at two cusps

For any nonzero future causal vector zeta, the height
`H_zeta(q)=-(zeta,q)` is bounded along a ray tending to a future null
ray xi only if `(zeta,xi)=0`. Indeed its leading term is a positive
multiple of `-exp(t)(zeta,xi)`. In Lorentz space, two future causal
vectors are orthogonal only when they are proportional null vectors.
Thus boundedness forces zeta to be proportional to xi.

For two distinct cusps xi and eta, no nonzero future causal zeta is
proportional to both. Every such height is therefore unbounded at at
least one cusp. Positive combinations of null cusp vectors do not
repair this problem. A spacelike height can sometimes be bounded, but
its decomposition at a positive-height point instead gives
`|x|^2=1+(zeta,zeta)/H_zeta(p)^2>1`; the unit energy budget is lost.

The missing general finite-volume ingredient is consequently a genuine
multicusp selection or boundary-loss estimate. Neither this lemma nor
the existing finite-vertex energy exclusions supply that estimate.

**Status:** written source, awaiting independent review. No scientific
execution, finite census, build or formalization was performed.
