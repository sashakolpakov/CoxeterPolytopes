# Independent local and complete-cone review of the 53-wall H9 chamber

**PASS.** Both the [local reviewer](../scripts/fresh_nonarithmetic_escape52_20261016_cap_review.py)
with its [certificate](fresh_nonarithmetic_escape52_20261016_cap_review.json)
and the [complete-cone reviewer](../scripts/fresh_nonarithmetic_escape52_20261016_cone_review.py)
with its [certificate](fresh_nonarithmetic_escape52_20261016_cone_review.json)
pass generation and fresh saved replay. They check the
[53-wall construction](fresh_nonarithmetic_escape52_20261016.md) using the
independently reconstructed 52-wall geometry, without producer arithmetic imports.

The local reviewer defines the new normal by a separate reverse interval
Gaussian inversion of the literal nine-root Gram matrix, followed by an
exact projection and positive square-root recipe. Certified nonzero pivots
justify the inverse. Projection identities establish squared norm 2, eight
right contacts and angle pi/4 to wall 30. Strict intervals certify the
remaining 43 dotted contacts, the old escape cut and the retained center.
All 53 facets are actual.

All nine possible supports for a new finite vertex are exhausted. Eight
contain a dotted pair. The remaining support is positive definite and
gives one actual new finite vertex. The new root has only A1 and I2(4)
connected finite-label subdiagrams, so it introduces no affine circuit;
absence of old affine circuits is an explicit inherited input.

The cone reviewer starts with the independently certified 2,457-ray
52-wall cone and performs a complete additional clip. It retains 2,307
rays, removes 150 exterior rays, and constructs 204 crossings. Each
crossing has a certified rank-eight common-covector minor and positive
convex coefficients. The final 133,083 ray/inequality checks give:

| Type | Count |
|---|---:|
| Finite vertices | 179 |
| Ideal vertices | 0 |
| Exterior rays | 2,332 |
| Opposite-sheet rays | 0 |

The old 178 finite vertices persist, and the unique new finite support
agrees with the local reviewer. An exact strict spacelike escape proves
infinite volume. The actual codimension-six face perpendicular to
0,1,2,3,5,6 has thirteen finite and seven exterior rays, with no ideal ray.
Its induced dihedral angles are not assumed to be Coxeter angles.

This is a complete cone calculation for this specified chamber, not a
finite-volume completion or a dimension existence theorem. In particular,
removing 150 old exterior rays introduces 203 new exterior rays, so broad
coverage of old escapes alone is not evidence of closure.

```sh
python3 scripts/fresh_nonarithmetic_escape52_20261016_cap_review.py --verify
python3 scripts/fresh_nonarithmetic_escape52_20261016_cone_review.py --verify
```
