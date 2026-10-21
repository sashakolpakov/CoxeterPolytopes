# A new escape after the25-wall discriminant19 polar completion

19 October 2026. Exact paper successor to the
[reviewed polar-cap construction](compact_h9_discriminant19_polar_caps_research_20261019.md),
SHA256792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b.
That stable source is unchanged. No scientific computation or build ran.

**Result, awaiting independent review.** The proposed25-wall H9 chamber
still has infinite volume. The new escape's polar direction is not the
direction of any reflection preserving the H19 lattice. This obstructs
the direct polar-cap mechanism at this next ray, not all other mirrors,
finite-volume completions, or the full lattice's reflectivity.

Use exactly the native form, wall indices, b=4+phi, D=1+5phi, and new
roots w,a,d,z in the preceding construction. Let

    omega1=C_E8^{-1}e1=(7,14,20,16,12,8,4,10),
    omega1²=14,
    M=8phi-5>0,
    lambda=D/M=(37+33phi)/79>0,
    X=(2,-phi,lambda*omega1).

The displayed integral vector verifies C_E8*omega1=e1 directly, so no
inverse computation is hidden.

## Complete feasibility against all25 walls

For any positive-control wall s=(x_s,y_s,z_s), the binary part of X is
the original timelike control p. Consequently

    (X,s)=-D*x_s+lambda*z_s,1
         =lambda*(z_s,1-M*x_s).                    (1)

Here z_s,1 is the second native E8 coordinate (zero-based index1),
not the second binary coordinate. The following residuals cover every
positive-control wall of the25-wall seed:

|Wall|z_s,1-M*x_s|
|---|---:|
|old9|5-4phi|
|old11|4phi-8|
|old12|1-phi|
|old13|4-3phi|
|old14|-9+5phi|
|old15|1-2phi|
|old16|-4+2phi|
|old17|0|
|old18,old19|-1|
|r|-M|
|k|-2-5phi|
|w|-4-6phi|
|a|6-6phi|
|d|-7-6phi|
|z|7-8phi|

Every nonzero entry is strictly negative, using 3/2<phi<5/3 where
needed. The zero-height wall old0 has product0. The eight E8 negative
simple roots give products -lambda at old2 and0 at the other seven.
Thus X satisfies all25 inequalities exactly. No extremality or census
of other rays is required.

## The escape is spacelike

The useful exact identities are

    M*phi³=11+14phi,
    D*phi³=b*phi^4=11+17phi.

Therefore

    X²=-2D+14lambda²
      =2*(33+32phi)/(11+14phi)²>0.                 (2)

For direct checking of the numerator,

    7*(11+17phi)²-D*(11+14phi)²=33+32phi.

Since the preserved center c is strict timelike, a sufficiently small
perturbation of X toward c is both strict for all25 walls and still
spacelike. The segment from c crosses the future light cone in the
interior of every wall halfspace. An open ideal patch is exposed, so
the chamber has infinite hyperbolic volume.

The25-facet/83-finite-vertex construction remains valid; its unresolved
volume status in the preserved predecessor is now settled negatively.
No new finite-volume or compact H9 example results.

## The polar line contains no integral reflection direction

The primitive reflection norm alphabet for the fixed H19 lattice,
up to unit rescaling, is {2,2b}; its proof is in the unchanged21-wall
report. It uses the even integral lattice, prime determinant ideal of
norm19, and the explicit arithmetic of Z[phi].

For a+bphi, the rational norm is a²+ab-b². Thus

    N(33+32phi)=1121=19*59,
    N(11+14phi)=79,
    N(X²)=4*19*59/79².

If a nonzero multiple u=alpha*X were a lattice reflection direction,
rescale its primitive root to canonical norm2 or2b. Scalar multiplication
changes the rational norm of the square by N(alpha)², a rational square.
But the ratios to the two canonical possibilities are

    N(X²)/N(2)  =19*59/79²,
    N(X²)/N(2b)=59/79².

Neither is a rational square. This is a contradiction. It excludes all
scalar multiples in K=Q(sqrt5), which is enough: a K-coordinate root
parallel to X necessarily has scalar in K because X's first coordinate
is2. No choice of primitive scaling, unit, norm class or height repairs
this polar line.

A different short or long lattice mirror could still cut X. Such a
mirror must be new relative to the current subgroup, pass every retained
wall, and address a larger part of the cone. Those conditions have not
been solved or excluded here. No global nonreflectivity, infinite-order
chamber symmetry, or whole completion obstruction follows from the
prime59 polar-line obstruction alone.

