# Root review: complete rank-at-most-four farthest forest

20 October 2026. Full written review and rederivation; no scientific
worker, enumeration, TeX or Lean run.

Reviewed source:
`results/compact_farthest_rank_four_forest_exclusion_20261020.md`,
SHA256 `68ad62a619603eb7e6c2acc5689c901dde3f9d0897a7c1580b3deb948f5832b6`.
I also read the complete three pinned component sources and the
orthogonal-opposite source; their bytes match the displayed pins.

**PASS for the pointwise theorem.** In dimension at least 15, the
chosen farthest vertex of a compact Coxeter polytope has an irreducible
finite component of rank at least five. This is not a numerical ceiling
for unrestricted compact or arithmetic compact polytopes.

The new dihedral argument is complete with arbitrary spectators. Low
energy forces one actual opposite, so its single-contact inequality is
beta^2(9d+7)/16<=1. This excludes old labels at least 8, leaves exactly
the displayed native cases for labels 4 through 7, and isolates the
old-I2(7)/new-3 exception. The double-contact length bound has the correct
direction and includes every spectator contribution. The two distinct
7 defects have the stated different heights and the same strict mate
energy floor greater than 9/100. All rational and radical comparisons
in this reduction check.

The rank-at-most-four finite list is exhaustive by the displayed cycle,
rank-three determinant and rank-four star/path calculations. No finite
search, hyperbolic polytope classification or field list is hidden in
the completeness assertion.

I checked both propagation interfaces against the full small-A source.
They are universal actual-wall inequalities, although that source uses
them for a smaller whole forest. Under the propagated half-tangent
floor, the new proof correctly uses the maximum of covering and
intrinsic estimates. H4 is the weakest component at rate 11/144 per
rank. This handles both the old-A2/new-7 and the small-height isolated
A1 cases without an unsupported covering claim for exceptional types.

In the old-I2(7)/new-3 case, disjoint spatial supports make its wall
incompatible with every short A1 wall, including J, and with the other
7 defect. Thus there is no uncharged second deficit. All remaining
blocks pay their rank/14 baseline. The final exact lower bounds
713/700 and 71/70, and the separate propagation inequalities, each
strictly exceed one at n=15 and increase thereafter. No n=14 assertion
or local-to-global dimension implication is being made.

Root co-developed the common actual-wall mechanism and the H3/H4
proof, and independently derived the new dihedral reduction and final
assembly alongside the author. This is a full-source audit, not an
independent-origin claim. Outside-author full-source review is a
separate record. The result remains supporting research: larger
A,B,D,E components prevent a uniform dimension conclusion, and no
additional manuscript or Lean expansion is warranted yet.

Batch closeout: the outside-author full review is saved as
`results/review_farthest_rank_four_forest_finite57_20261020.md`,
SHA256 `9fa74c3a42a9dced98b4c2ce7181902e5ee4f80ef6d270b382952085d3f10473`.
Catalogue JSON and all76 research-note paths pass, as do404 local
documentation links and seven current source/review hashes. Three
previously frozen primary TeX sources were also hash-checked unchanged.
There remain12 catalogue research entries and no added primary paper.
The proof lock is available, process inspection finds zero Lean/TeX
processes, and free disk is177.6GiB. No scientific or build job ran in
this batch; there are no new disposable build artifacts to remove.
