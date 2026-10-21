# Independent complete review of the 36-wall H9 cone

**PASS:** the specified cone has exactly **100 finite rays, no ideal rays,
1,068 hyperideal rays, and no opposite-sheet rays**. All four clipping stages
agree with the producer. The adjacency filter and active-mask uniqueness are
justified; no completeness discrepancy was found.

Artifacts: [independent source](../scripts/fresh_nonarithmetic_escape35_20261016_cone_review.py),
[literal certificate](fresh_nonarithmetic_escape35_20261016_cone_review.json),
and [producer report](fresh_nonarithmetic_escape35_20261016_cone.md).

This is a complete polyhedral cone census, with the independently complete
32-wall cone as its explicit starting dependency. The separately reviewed
strict spacelike escape still proves **infinite volume**.

## Independent arithmetic and reconstruction

No producer arithmetic is imported. All 972 initial rays are reconstructed
from their active equations in a separately instantiated quartic primitive
field. Reverse-column exact elimination verifies active rank nine and every
initial wall sign. The starting census is bound to its frozen independent
review and checked literally.

Each ray is normalized by the fixed functional

`f(v) = -sum_{i<32}(v.n_i) = 1`.

This functional is positive on every nonzero vector of the old cone: all
summands are nonnegative, and their simultaneous vanishing would put the
vector in the common kernel of the old covectors. Their exact rank is ten.

All later rays are represented by exact positive-combination expression
trees. Their coordinates and every strict sign are enclosed with an
independently implemented, outward-rounded 512-bit rational interval engine.
The physical octic primitive element and the three positive radical roots
are evaluated directly. No assumption of independence or minimal degree for
the three added radicals is needed.

For a clipping wall with endpoint values `a<0<b`, the new normalized ray is

`w = (b*u-a*v)/(b-a)`.

Its coefficients are strictly positive and sum to one. Its new-wall contact
vanishes algebraically. Its old zero set is exactly the intersection of the
endpoint zero sets, because a positive combination cannot cancel two
nonpositive wall contacts. Thus zero identities are established exactly,
not inferred from intervals containing zero.

## Why the adjacency test and filter are complete

The common strict center makes every successive cone full-dimensional in
the ten-dimensional ambient vector space. Exact rank ten of the old
covectors makes the cone pointed; these properties survive every clipping.

For two distinct extreme rays `u,v`, let `S` be their common active walls.
The sum `u+v` is strictly negative on every wall outside `S`. Therefore it
lies in the relative interior of their common-zero face, whose dimension is

`10 - rank(S)`.

The rays are adjacent precisely when this rank is eight. Fewer than eight
common walls cannot have rank eight, so the producer's cardinality filter
cannot discard an edge.

The review uses this **rank criterion**, independently of the producer's
test based on counting rays in a face. Every accepted crossing pair receives
a certified nonzero eight-by-eight minor, obtained by interval elimination
with complete pivoting. The two distinct common-null rays also give the
matching rank upper bound eight. All 64 remaining candidates occur in the
first clip and are rejected by exact quartic rank calculations.

The producer's criterion is consequently sound as well: in a pointed
polyhedral cone with a complete extreme-ray list, a face contains exactly
two distinct extreme rays if and only if it is two-dimensional.

Halfspace clipping retains the old negative rays, retains any zero rays,
and adds exactly the intersections with crossing two-faces. No old ray lies
on any of the four clipping walls here. This propagates the inherited
complete list through all four steps.

## Why active-mask deduplication is valid

Every initial active set has exact rank nine. For a new ray, the common old
constraints have rank eight, and the clipping equation is independent of
them because it has opposite signs on the two edge generators. The new
active set therefore has rank nine as well.

An active mask hence determines a one-dimensional kernel. Pointedness fixes
its feasible orientation, and `f=1` fixes its vector. Equal masks cannot
represent different rays. The independent code asserts uniqueness; it never
silently drops a duplicate candidate.

## Replayed census

| Added wall | Negative old rays | Positive old rays | New crossing rays | Finite | Hyperideal |
| --- | ---: | ---: | ---: | ---: | ---: |
| 32 | 940 | 32 | 112 | 92 | 960 |
| 33 | 1,047 | 5 | 25 | 95 | 977 |
| 34 | 927 | 145 | 167 | 98 | 996 |
| 35 | 1,045 | 49 | 123 | 100 | 1,068 |

There are 427 certified crossing edges in total. The complete final support
lists agree with the producer, and the 100 finite supports agree with the
separate local geometric review. All 42,048 final ray/wall incidences are
checked, together with strict norm signs and the required time orientation.
The certificate saves final coordinate bounds, norm bounds, all rank minors,
all rejected-rank witnesses, and the exact parent-support structure of every
new ray.

## Validation and scope

Generation and a fresh saved-output replay pass. Recursive source manifests,
direct hashes, source syntax and report links are checked. Replay command:

```sh
python3 scripts/fresh_nonarithmetic_escape35_20261016_cone_review.py --verify
```

The initial 32-wall cone's completeness is inherited from its pinned
independent review; this package reconstructs every starting ray and checks
all subsequent clipping steps. It does not repeat the earlier construction
of that starting cone. The result concerns these specific 36 walls and does
not certify a finite-volume completion or a new-dimensional polytope.
