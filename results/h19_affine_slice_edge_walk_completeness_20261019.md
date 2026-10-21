# A bounded affine slice and a criterion for a complete H19 edge walk

19 October2026. Root paper lemma; no enumeration or proof program ran.
This gives a sufficient completeness criterion for a future exact walk.
It does not assert that the current28-wall walk has closed, or that the
chamber is compact. The standard Euclidean polytope convex-hull theorem
and connectedness of its vertex graph are explicit external inputs.

## The retained cone has a bounded x=1 slice

Use the native physical coordinates(x,y,z) and Gram from the H19 reports.
Let C0 be the cone cut out by any collection of walls retaining old0,
the eight negative E8 simples, old9 and r. All current chambers do so.
The inequalities from old0 and r imply

    phi*x+2y>=0,  y<=0,

because r's covector is(0,b;0) and b>0. Thus x>=0. If x=0 then y=0.
The eight simple-wall inequalities say C_E8*z>=0 coordinatewise.

The native old9 root is

    old9=(1,0;phi*(2,4,6,5,4,3,2,3)).

This is the literal retained input, not an inferred generic positive
root. At x=y=0 its inequality becomes

    phi*(2,4,6,5,4,3,2,3) . (C_E8*z) <= 0.

Every coefficient is strictly positive and every entry of C_E8*z is
nonnegative. Thus C_E8*z=0, whence z=0. Therefore

    C0 intersect {x=0} = {0}.                       (1)

The slice S=C0 intersect{x=1} is closed and is a Euclidean polyhedron.
It is bounded: otherwise an unbounded sequence in S, normalized by
its Euclidean norm, has a subsequence converging on the unit sphere
to a nonzero vector of C0 with x=0, contradicting(1). Thus S is a
bounded Euclidean polytope. Its strict interior is supplied by the
existing timelike centre after scaling to x=1.

## What a complete edge walk would have to certify

At a finite vertex the nine active facet normals are independent.
Indeed they lie in the positive-definite orthogonal space of the
vertex. Splitting a dependence into disjoint positive combinations,
acuteness makes their common square nonpositive, hence zero; their
products with a strict interior point rule out a nonzero positive
combination equal to zero. There can therefore be exactly nine
active facets at a finite vertex, and it is simple.

In the x=1 slice, write N for the9-by-9 matrix of active covectors
restricted to the other nine coordinates. It is invertible: a vector
in its kernel with dx=0 would be orthogonal to all nine normals and
therefore proportional to their unique timelike vertex direction,
whose x coordinate is nonzero. Each edge direction away from this
vertex is consequently obtained from

    dx=0,  N*d=-e_i,  1<=i<=9.

For every other wall with positive slope, its exact slack divided by
that slope gives a positive intersection parameter. The least one
is the actual next Euclidean edge endpoint. Boundedness rules out an
edge direction with no such intersection. Every inequality and the
active set at each resulting endpoint must be checked, not inferred
from a numerical proposed neighbour.

Starting from p, suppose a finite exact walk records all nine edge
endpoints at every discovered vertex, each endpoint has negative
Lorentz square, and its pending queue eventually empties. Then its
finite vertex set is a nonempty component of the FULL Euclidean vertex
graph of S: no edge from one of its vertices leaves the set. The
external connectedness theorem implies that every vertex of S has
been reached. A stopped walk with a nonempty queue, a vertex cap,
timeout or unknown incidence does not establish this conclusion.

## Why complete finite closure would imply compactness

Every reached timelike point of the x=1 slice is future relative to
p, since(p,v)=-D<0. A convex combination of future timelike vectors
is again future timelike. By the external convex-hull theorem, if
all Euclidean vertices were reached and timelike then every point
of S would be future timelike. On the compact set S the Lorentz
square would be uniformly strictly negative. Renormalizing S to the
hyperboloid would therefore give a compact hyperbolic polytope.

Conversely, a feasible spacelike endpoint together with a strict
timelike centre supplies an open ideal-boundary patch by the previous
escape argument and proves infinite volume. A feasible ideal endpoint
requires its own analysis. These alternatives must be kept separate
from mere finite queue exhaustion or an incomplete diagnostic.

The current28-wall candidate has not yet passed this complete-walk
criterion. No new finite-volume or compact dimension-nine example
is asserted by this lemma.
