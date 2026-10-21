# Scoped review of the 92-wall coupled lift

**PASS.** The frozen [construction](fresh_h20_coupled_lift_20261014.md), [primary source](../scripts/fresh_h20_coupled_lift_20261014.py) and [certificate](fresh_h20_coupled_lift_20261014.json) have no issue within the reviewed scope. A fresh primary `--verify` passed all 2,416 new pair gates. The chamber has 92 actual facets and sixteen specified actual finite vertices; its exact strict spacelike escape still proves infinite volume.

The separate [Fraction checker](../scripts/fresh_h20_coupled_lift_20261014_review.py) imports no producer or project arithmetic. Its [hash-bound result](fresh_h20_coupled_lift_20261014_review.json) independently checks the following from the frozen literal coordinates:

- All 92 root norms and strict common-center inequalities.
- Four rational positive-definite rank-18 frames, the orthogonal positive vectors `m,w`, both cap norms and the timelike vertex projection.
- All four height-flip form isometries, their old-wall permutations and new-root transports, and all 240 old-wall signs at the four old escape images.
- All 1,472 inequalities at the sixteen new vertices, including their exact twenty-wall zero sets, negative norms and future orientation.
- All 96 strict new-wall inequalities at the three inherited cusp rays.
- The final escape's positive norm, future orientation and all 92 strict inequalities.

The independent generation and saved replay both passed. This is a separate geometry implementation, **not a second full pair census**. The old form signature, old-old Coxeter gates and old bounded-link proofs remain explicit inherited dependencies. The primary replay checks the predecessor manifests; the review certificate binds its own source and the frozen primary inputs.

## Formula, interval and geometric review

The identities `m·w=0` and `d*w²=m²=q>0` prove that `m±sqrt(d)*w` each have square `2q` and are mutually orthogonal. Together with the rank-18 positive support they span a positive rank-20 space. Their common orthogonal projection of the strict center is timelike and satisfies all nonincident inequalities strictly, proving an actual finite vertex rather than merely twenty zero pairings. The source's exact coefficient identities supply the right angles; interval overlap is not used to assert equality.

The strict interval gates are sound. The integer square-root bounds enclose each positive radical, signed coefficient rounding is outward, and interval dot products enclose the exact cross-field pairing. A negative upper bound whose square exceeds the exact product of the two positive root norms proves strict dottedness. The recorded counts are 592 right and 1,824 dotted new pairs. All old and new facets are actual: projection of the common strict timelike point onto any wall remains timelike and stays strict on every other wall because all cross pairings are nonpositive.

The explicit new-cap neighbor lists justify the stated 320 outgoing germs. Every such edge retains a new cap, and every wall outside its starting twenty-wall corner is strictly dotted to a retained cap. Such a wall cannot contain a finite or null point of that edge. Each germ therefore reaches the light cone without another inequality becoming active. It has an ideal-boundary limit, but no second finite vertex or bounded-link ideal vertex. This is not a complete edge or vertex census of the full chamber.

All new caps are strictly negative at the three inherited null rays, independently checked above; hence their previous bounded links persist. The review makes no claim that these are all cusps. The final future spacelike norm is exactly

```text
1411253767129583933/200000000000000000 > 0.
```

Consequently neither the new finite corners nor the repaired old escape images imply finite volume.

## Raw-height family proof

The separate all-real obstruction in the report is also sound with its stated assumptions. For each distinguished short root, its unchanged norm-two long neighbor forces the lifted norm to be one or two. Acuteness and the zero-in-convex-hull condition then force either a single zero height or an antipodal unit pair. The three old distinguished roots are pairwise parallel; finite cross-angle gates force their nonzero height axes to be equal or perpendicular. Some pair therefore has height span of dimension at most one.

Their sum is an old feasible future null ray, strictly negative on every other short direction. A transverse positive unit vector orthogonal to the chosen height span gives a feasible spacelike point after a sufficiently large null shift. Because the height family is finite, one shift works for all its inequalities. Perturbation toward the old strict center makes the remaining equalities strict while preserving positive norm and future orientation. This excludes that raw-height family only: all old long walls are retained, every old short height set has zero in its convex hull, and additional projected directions are not allowed.

The conditional non-quasi-arithmetic statement remains inherited from the predecessor and applies to a hypothetical finite-volume completion retaining its witness facets. This review does not call the current infinite-volume reflection group a lattice or assert a new dimension record.

```sh
python scripts/fresh_h20_coupled_lift_20261014.py --verify
python scripts/fresh_h20_coupled_lift_20261014_review.py --verify
```
