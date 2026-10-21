# Independent complete review of the49-wall H9 cone

**PASS:** the cone has exactly178 finite rays,0 ideal rays,2,199 hyperideal
rays and0 opposite-sheet rays. Its finite graph has409 edges and784 open
finite-edge germs. A renewed exact strict spacelike escape proves infinite
volume. This agrees with the [primary construction](fresh_nonarithmetic_common_escape_20261016.md).

Artifacts: [independent source](../scripts/fresh_nonarithmetic_common_escape_20261016_cone_review.py)
and [literal certificate](fresh_nonarithmetic_common_escape_20261016_cone_review.json).
The review starts from the [independently complete48-wall cone](fresh_nonarithmetic_joint_regions_20261016_cone_review.md)
and uses the [independently reconstructed49th wall](fresh_nonarithmetic_common_escape_20261016_cap_review.md).
It imports no new producer arithmetic and performs a single exhaustive clip,
rather than repeating the producer's thirteen clips from36 walls.

## Complete clipping and zero incidences

The2,281 initial rays have independently certified coordinate bounds and
active rank9. The new wall has negative pairing with2,180 rays, positive
pairing with94, and zero pairing with7. Ambiguous interval signs are resolved
by exact nested-radical reconstruction through the prior independent clipping
provenance. Zero contacts enlarge the active masks explicitly; they are not
discarded or replaced by a floating tolerance.

For every negative/positive pair, fewer than eight common active walls
cannot define a two-dimensional face. All190 remaining pairs have a certified
nonzero8x8 common-covector minor. The two independent common-null rays give
the matching rank upper bound, so each pair is an actual crossing edge.
No rank candidate is left unresolved. This criterion is independent of the
producer's test that the common-zero face contains exactly those two rays.

Each new ray is the positive convex combination of its endpoints that
annuls the clipping wall. Its inherited zero set is exactly the intersection
of the old zero sets: positive coefficients cannot cancel two nonpositive
wall pairings. The new wall adds an independent equation, preserving active
rank9. The review asserts uniqueness of all resulting masks. Pointedness,
full dimension and the common strict center are inherited and independently
checked in the cap package. These facts make the clip exhaustive.

The final census is therefore

| Type | Count |
|---|---:|
| Finite |178|
| Ideal |0|
| Hyperideal |2,199|
| Opposite sheet |0|

All116,473 final wall/ray incidences are checked. Strict norm bounds classify
every ray, and every final active support agrees with the producer. The
172 old finite vertices survive and six new ones appear.

## Finite graph and renewed escape

Each finite ray is simple with nine active walls. The reviewer independently
recounts finite adjacencies by common rank8, checking a nonzero minor for
each of the409 edges. Consequently the number of open finite-edge germs is
`9*178-2*409=784`.

The replacement escape starts from the hyperideal ray with active walls
`[0,1,2,3,5,6,7,30,33]` and adds `H/67108864`. The reviewer reconstructs its
exact nested-field coordinates from the saved literal vector, checks every
active zero and every other ray inequality, verifies the perturbation formula,
and proves its square positive, center pairing negative and all49 wall
pairings strictly negative. This supplies an open feasible ideal-boundary
patch and proves infinite volume.

Generation and fresh saved replay pass. The old48-wall completeness is an
explicit independently reviewed dependency; the current package independently
checks its next clipping step and every final inequality. The result does
not certify a finite-volume completion or an example in a new dimension.

```sh
python3 scripts/fresh_nonarithmetic_common_escape_20261016_cone_review.py --verify
```
