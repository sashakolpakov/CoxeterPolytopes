# Independent review of the25-wall escape

19 October 2026. **Exact written review: PASS.** No matrix program,
scientific computation, root search, optimizer, Lean job or build ran.
The reviewed source is
`results/compact_h9_discriminant19_twentyfive_escape_20261019.md`, SHA256
`1cc74516c9e22b7150c64a95b4c4251dec9d3c02f238dbb4c946ef937e37c096`.
It and the preceding polar-cap packet were preserved unchanged.

The new source supersedes only the predecessor's historical unresolved
volume status: its25-wall chamber is infinite-volume. The construction of
25 actual facets and at least83 actual finite vertices remains valid.
There is no new finite-volume or compact H9 example.

## Feasibility includes every wall

I checked the displayed native vector

    omega1=(7,14,20,16,12,8,4,10)

against the E8 adjacency01,12,23,34,45,56,27. Its C-image is e1 and its
square is14. For M=8phi-5, direct rationalization gives

    D/M=(37+33phi)/79>0.

The binary part of X equals the old timelike control p, so its product
with any archived or new root is exactly

    (X,s)=lambda*(z_s,1-M*x_s).

I independently substituted EVERY positive-height native row from the
pinned prefix and all six later roots r,k,w,a,d,z. All residuals agree
with the complete table, including the single old17 equality. For
example M phi^2=3+11phi, M phi^3=11+14phi and M phi^4=14+25phi give
the old12/13, old15--19 and k/w/z rows directly. The remaining9 old
height-zero normals give zero, except old2 with product-lambda.
Thus all25 halfspaces, rather than a selected subset, contain X.

## Positive norm really gives infinite volume

The identities M phi^3=11+14phi and D phi^3=11+17phi are correct.
Expanding the latter two squares gives

    7(11+17phi)^2-D(11+14phi)^2=33+32phi,
    X^2=2(33+32phi)/(11+14phi)^2>0.

The preserved common center c is strict timelike by the independently
reviewed polar-cap packet. A small positive perturbation of X towards c
therefore stays spacelike and becomes strict against all25 walls. The
future component is the relevant one: (p,X)=-2D<0, and the segment from
the future timelike center to this point first crosses the future light
cone with every inequality strict. A neighborhood of that ideal point
is feasible, so the25-wall chamber has infinite hyperbolic volume.
No claim that X is an extreme ray or the only escape is needed.

## The exact polar-line obstruction is correctly scoped

The quadratic norm computations give

    N(33+32phi)=19*59,
    N(11+14phi)=79,
    N(X^2)=4*19*59/79^2.

The unchanged H19 reflection-norm theorem permits primitive reflection
norm classes2 and2b, up to unit rescaling. Their rational norms are4
and4*19. The two resulting ratios to N(X^2) are nonsquares, owing to
the odd prime59 valuation (and also19 in the first class). A scalar
in K changes N(X^2) by a rational square. No K-scaling can therefore
give either canonical reflection norm. This includes arbitrary
primitive/unit choices, rather than testing one displayed denominator.

There is no unexamined real scaling issue: since X's first coordinate
is2, every K-coordinate root on its real line has scalar in K. The
exact polar direction contains no H19 lattice reflector.

This uses the earlier complete primitive norm-class argument as an
explicit dependency; no new mixed-norm census was performed. The
conclusion excludes this polar line ONLY. Other lattice mirrors may
cut X, and further finite-volume completions are neither produced nor
excluded. In particular no nonreflectivity or whole-family completion
obstruction follows from the prime59 factor alone.
