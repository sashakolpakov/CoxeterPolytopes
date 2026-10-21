# Independent review of the Lambert-quadrilateral cap

**PASS:** the new wall cuts the previous strict escape, gives 36 actual
facets, and closes a compact Lambert quadrilateral. The complete finite
census becomes 100 vertices and 266 edges. A different exact strict escape
still proves infinite volume. The complete ambient cone is a separate
calculation and is not assumed by this review.

The [independent checker](../scripts/fresh_nonarithmetic_escape35_20261016_review.py)
and [certificate](fresh_nonarithmetic_escape35_20261016_review.json) use a
separate outer radical class over the frozen independent implementations
used by the [35-wall review](fresh_nonarithmetic_escape34_20261016_review.md).
No new-cap producer is imported. Reverse-order elimination independently
reconstructs the particular solution rather than just accepting its saved
coordinates.

The exact normal is `p-sqrt(10+68sqrt5/15)*k`, where `k` is the final coordinate
vector, `p^2=-5-5sqrt5/2`, `k^2=sqrt5/4` and `(p,k)=0`. Its norm is `2/3`.
This normalization makes the anchor product minus one and the Cartan cycle
three, giving angle `pi/6` without adjoining a separate square root of three.
The other 34 new pairings consist of eight right and 26 dotted contacts.

Direct projection identities verify all 36 relative-interior facet
witnesses. Every old finite vertex remains strict against the new cap.
Of the complete nine candidate supports on the new facet, exactly two
are positive definite and give actual finite vertices. Combining these
with the independently checked predecessor census gives exactly 100
finite vertices, 266 finite edges and 368 open finite-corner germs.

The seven-root white support is `[0,1,2,4,5,6,9]`; boundary walls are
`[7,8,11,35]`. The reviewer projects them into the perpendicular Lorentz
three-space and examines all six pairs of boundary lines. Exactly four
intersections have a feasible orientation. Each is future timelike and
satisfies every retained wall inequality with the required active support.
Both remaining intersections are spacelike and have neither orientation
feasible. The boundary normals span the full three-space, so this is a
pointed cone whose complete extreme rays are the four finite ones. It
therefore defines a complete compact quadrilateral. The projected Gram
cycles at its corners independently give three right angles and one
angle `pi/6`.

A replacement escape is reconstructed from the old support
`[0,1,2,4,5,6,7,10,27]`, normalized against the center and perturbed by
`H/8192`. Its square is positive, its center pairing is negative, and
all 36 wall pairings are strictly negative. Its projection perpendicular
to the seven white roots is future timelike and strictly inside the
compact quadrilateral. Hence a coherent acute cap perpendicular to all
those roots cannot cut this escape: it is negative on the finite vertex
rays and thus on their positive combinations. This is the same
finite-vertex preservation argument proved in the preceding review,
with the current compact face replacing the tetrahedron.

Generation and fresh saved replay passed. The current source and data
hashes are bound. There is no claim of a new finite-volume polytope or
new construction dimension.

```sh
python3 scripts/fresh_nonarithmetic_escape35_20261016_review.py
python3 scripts/fresh_nonarithmetic_escape35_20261016_review.py --verify
```
