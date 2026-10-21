# Independent source review of the unified ideal right-angled four-polytope paper

18 October 2026. **PASS for the written consolidation: no mathematical flaw or
circular dependency found.** This is a full mathematical/source read, not a
proof by agreement between reviewers. I did not build the PDF, run a checker,
execute a finite census, or run Lean during this review.

Reviewed source: `manuscripts/ideal_ra4_volume_and_structure.tex`, all 1036
lines, SHA256
`e503ad2cb8edafd827f735a62568a3372407e9dfdaf3f50d09b05fb3bc75f025`.
I also read the gap, equality, restricted decomposition and path-rigidity
predecessor reports and the existing literature-scope audit. Frozen inputs
were not changed.

## Dependency and geometric checks

The unrestricted gap uses the volume/incidence identities, antiprism bound,
nine-vertex classification and proved opposite-cusp injection. The equality
reduction separately adds the ten-vertex Figure 7 classification and its
quadrilateral matching. It establishes eight D9 facets before invoking the
restricted Oct/D theorem. The N=7 cuboctahedral branch includes the embedded
medial recovery and the nonsquare cusp calculation; its exclusion is not
inferred from the desired decomposition.

The common cutting proof retains the steps most easily lost in a merger:
the belt link is a single circle; the elementary Euler bound gives eight
triangles without assuming a spherical surface; the six-cusp graph and its
eight actual facet triangles give the cube; equal short widths apply to
unit normals; their affine normal identity produces an actual section by
convexity. Outside facets cannot touch that section at a finite point. Its
center is strict on every wall, so both halves have interior. Cutting
preserves entirely ideal right angles and removes exactly eight D9 facets.
The disjointness argument includes ideal vertices, as required by the later
opposite-cut argument.

All-octahedral recognition is proved by full-facet propagation rather than
assuming the restricted conclusion. The D4 common-neighbor argument forces
two cuts in one cell to be opposite: otherwise one original facet contains
three octahedra with disjoint relative interiors, contradicting the allowed
facet volumes. Thus the decomposition tree is a path. The no-twist argument
and the reflected double's two-sided halfspace containment are both present.

The normalizations and counts agree: N=8k counts D9 facets, k counts cuts,
k+1 counts cells, and W=12k. The general face counts specialize to the
double's (42,168,164,38), with 30 ordinary and eight doubled octahedral facets.

## Scope and remaining external dependencies

The source explicitly limits the chain theorem to the Oct/D facet spectrum,
classifies existing chains without asserting every length exists, and makes
no universal 24-cell decomposition or stronger volume-divisibility claim.
It neither claims publication priority nor treats finite model certificates
as a universal geometric proof. The shared cube/section/cut proof is given
once; I found no needless duplicated geometric proof.

Published volume/incidence formulas, the antiprism equality case, the complete
nine- and ten-vertex small-facet inputs, and three-dimensional right-angled
rigidity remain external assumptions. I did not freshly inspect Kolpakov's
Figure 7 in this consolidation review. No end-to-end formal verification of
these hyperbolic arguments is asserted. PDF layout, compilation, current
catalogue execution and Lean interfaces remain outside this source review.
