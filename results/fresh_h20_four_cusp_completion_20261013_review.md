# Independent source review of the H20 224–226 checkpoint

**PASS.** I reviewed the proof-sensitive source and the [lane report](fresh_h20_four_cusp_completion_20261013.md), then independently invoked all four frozen `--verify` commands. All passed. This is a mathematical/source review and fresh execution of the supplied exact arithmetic, not a second implementation of all field operations or a new construction search.

The [224-wall checker](../scripts/fresh_h20_four_cusp_completion_20261013.py) verifies all1,971 new pairs. At each of the four specified cusps, it covers the entire incident set with strictly positive weights whose vector sum is a positive multiple of the null ray. The actual normal span has rank20, hence quotient rank19. A horospherical recession vector is nonpositive on all these normals; the positive relation forces every pairing to vanish, and the quotient span forces the vector to vanish. This proves boundedness, rather than inferring it from rank alone. The common strict timelike center supplies a full-dimensional link interior.

For each of the nine oriented families `r(t)=-theta+t*delta`, an explicit retained wall forces `t>=2`. The saved spacelike escape pairs strictly negatively with both `r(2)` and `delta`, so it survives every real allowed shift. The stated obstruction therefore covers any finite collection from these nine families retaining the required acute wall relations. It does not cover different cusp profiles or wall replacement.

The [225-wall checker](../scripts/fresh_h20_four_cusp_completion_20261013_gap_exact.py) reconstructs the rational norm-one cap and all224 of its pairs. Its fixed17-block obstruction uses an exact feasible spacelike vector equal, modulo that block, to a strictly positive combination of four retained normals. Every new normal orthogonal to the block and nonpositive against those four walls is nonpositive on this vector. For any proposed finite full-dimensional completion, perturbing toward that completion's own strict interior gives a strict spacelike escape. This correctly avoids requiring the old saved center to remain feasible under arbitrary additions.

The [226-wall checker](../scripts/fresh_h20_four_cusp_completion_20261013_next_exact.py) verifies norm one,19 right contacts, two pi/4 anchors and204 dotted pairs for the changed-block cap. The PD19 block and its orthogonal new normal give an actual rank20 finite vertex: the projected point is future timelike, has precisely the claimed zero inequalities, and is strict on every other wall. Both changed-profile caps are strict at all four completed cusps, so those links are preserved. The final exact strict escape has square `1045865484809280029/10^20>0`; the chamber remains infinite-volume.

The [new-cusp diagnostic](../scripts/fresh_h20_four_cusp_completion_20261013_new_cusps.py) also passed. Both rays are feasible future null vectors, with every226-wall sign checked. Their links have affine/spherical quotient ranks `(1,18)` and `(2,15)`. Explicit positive-norm tangent recession vectors establish unboundedness, including the first case where the quotient already has full rank19. A small move along such a vector and then toward the strict center produces a strict spacelike escape. This is a diagnosis of two specified links, not a complete cusp census.

The actual-facet argument is sound throughout: nonpositive pairings and a common strict timelike center imply that its projection to any retained wall is timelike and strict on all other walls. No finite-volume or new-dimension example is inferred from the bounded links or finite vertex.

Replayed from the repository root:

```sh
python scripts/fresh_h20_four_cusp_completion_20261013.py --verify
python scripts/fresh_h20_four_cusp_completion_20261013_gap_exact.py --verify
python scripts/fresh_h20_four_cusp_completion_20261013_next_exact.py --verify
python scripts/fresh_h20_four_cusp_completion_20261013_new_cusps.py --verify
```
