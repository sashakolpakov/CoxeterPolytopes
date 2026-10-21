# Independent structural review of the discriminant19 polar-cap seed

19 October 2026. **PASS within the scope below; no structural flaw found.**
This was source reading and hand mathematics only. No arithmetic checker,
search, census, Lean process or build was run. Streaming source hashes
were used only to identify the reviewed bytes.

Reviewed report:
`results/compact_h9_discriminant19_polar_caps_research_20261019.md`, SHA256
`792634e935866ed78d95d0a93c04fdd7a8484f3605b5c50a66c1d376c9599f1b`.

I read the complete report, its 21-wall predecessor, the literal native
prefix coordinates and the archived Bugaenko H8 coordinates/Gram/vertex
list. This review concentrates on Sections 2--5: integral transport,
complete facets, distinct finite vertices and genuinely additional
mirrors. The other independent review owns full algebraic checking of
the entire pair ledger. I did not recompute every old/new ledger entry.

## 1. Integral transport is sufficient, without claiming an ambient lattice isomorphism

In the archived H8 coordinates the eight roots numbered 1,...,7,9 are
negative coordinate vectors in a permutation of the eight spatial
positions. Root 8 has control coordinate 1 and integral spatial part.
These nine roots therefore form an O=Z[phi] basis of the source lattice:
the determinant is a unit, rather than merely a nonzero field element.

The proposed native images have the same nine-root Gram matrix. The
E8 chain is

    old5--old4--old3--old2--old1--old13--old0,

with the branch old3--old8. Root old9 meets only old0, with product
-phi. This is the literal source prefix: chain 1--...--7, branch 3--9,
and the 7--8 product -phi. Thus the map extends uniquely to a Lorentz
isometry of the nine-dimensional spaces.

Every source integral vector is an integral combination of that basis.
Since every native basis image is integral, the images of source roots
0 and 10 are integral native vectors. No assertion that the target
basis is the entire primitive lattice in w-perpendicular space is
needed. Nor is an extension of this isometry to an integral isometry
of the ambient rank-ten lattice needed.

The full native reflection claim has a separate valid justification:
the displayed a,d have norm 2 and integral native coordinates, and
the native Gram is integral. Hence their reflection formulas
v -> v-(v,a)*a and v -> v-(v,d)*d preserve the entire native lattice.
For w of norm 2b, Gw in b*O^10 gives the same conclusion after
division by b; its unit y coordinate establishes primitivity.

I independently checked the displayed sparse covectors Ga,Gd from
the native E8 tree, their norms 2, their products a.w=d.w=a.d=0,
and the norm w^2=2b. I also checked the new strict-center products
for a,d,w by summing their spatial coordinates and using Gp=(-D,0).
These hand checks agree with the report.

## 2. Both compact H8 facets are complete and correctly oriented

For the new w-facet the two missing source roots have exactly the
required additional Gram products: a has product -phi with old5
(source 1), while d has product -1 with old1 (source 5) and -2phi
with old8 (source 9). Their mutual product is zero. Together with
the prefix, this matches the literal archived 11-root Bugaenko Gram.

For the original r-facet, the spatial permutation in the predecessor
takes source root 8 to old9, source root 0 to old11 and source root
10 to k. I checked these coordinate identifications against the
archived positive-height roots. The new w,a,d are strictly
ultraparallel to r, so none introduces a new ridge on that facet.

There is an orientation point worth making explicit. A hyperplane
ultraparallel to a facet hyperplane cannot cut it, but one must still
know which of its halfspaces contains the whole hyperplane. Here
the acute projection of the common strict center onto r or w is
timelike and strict against every other wall. That witness places
the entire connected facet hyperplane in the correct halfspace of
each nonintersecting wall. Thus the new facets are the full compact
11-wall chambers, rather than a subregion or an excluded component.

This conclusion inherits the archived H8 chamber's compactness and
41-vertex completeness. I read its literal data and matching Gram;
I did not rerun its exact Vinberg completion or vertex census.
Its source hash independently matches
`c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c`.

## 3. The lower bound of 83 finite vertices does not double-count

Each intrinsic elliptic eight-wall vertex of an H8 facet, with its
positive facet normal adjoined orthogonally, gives an elliptic
nine-wall ambient vertex. The nonbounding walls are strict on the
whole facet, so these are actual vertices of the seed chamber.

The r,w hyperplanes are strictly ultraparallel: their normalized
product is -phi^5/2<-1. Their two 41-element finite vertex sets are
therefore disjoint. The initial finite vertex p remains feasible,
and its pairings with r,w are respectively -b*phi and -2D*phi^4,
both strict. It lies on neither facet. Thus the claimed lower bound
is exactly the justified sum 41+41+1, not a reused ambient census.
It is not asserted to be a complete vertex list.

## 4. The mirrors are new relative to the specified seed subgroups

For any new wall u with positive norm and the given strict center c,

    c_u=c-((c,u)/u^2)*u

is timelike and lies on u. Because (c,u)<0 and every old product
(u,v)<=0, all old inequalities remain strict at c_u. Hence u's
mirror crosses the interior of the previous seed chamber.

The acute Coxeter seed is a fundamental chamber for the group
generated by its walls, regardless of whether its volume is finite.
A mirror of that generated subgroup cannot cross its chamber
interior. Thus w, then a,d (also mutually orthogonal), really are
additional mirrors relative to the preceding seed groups. This
argument addresses subgroup membership geometrically; it does not
mistake an arbitrary linear combination or an old reflection word
for a new cap. They remain integral reflections of the same ambient
lattice, as intended.

## 5. Volume and verification limits remain correctly stated

The displayed q has square 2*(115-71phi)>0. Its feasibility follows
from the report's complete z-pair ledger and q=phi^(-4)*z. A small
perturbation toward the strict timelike center stays spacelike and
becomes strict on all 24 walls. The segment from the center therefore
meets the future light cone inside every halfspace, giving an open
ideal patch and infinite volume. The 83 finite vertices do not
contradict this escape.

The added 25th polar wall removes q and preserves the two facets by
the same oriented-ultraparallel argument. Nothing reviewed here proves
that the 25-wall chamber has finite volume, that its z-facet is
complete, or that all further escapes are gone. The report explicitly
leaves those questions open. No compact H9 example or record follows.

The structural claims in Sections 2--5 pass this independent read,
conditional on the stated full pair ledger and the inherited compact
H8 certificate. This is a written mathematical review, not a fresh
execution or a replacement for either exact replay.
