# Hostile audit of the Jacobian--Kummer `H12` doubling cube

Status date: 2026-08-08.

## Verdict

The four seed facets `31`, `58`, `59`, and `68` independently generate a
faithful `C2^4` Coxeter doubling cube.  Every one of its fifteen nonempty
subset chambers is an exact finite-volume noncompact hyperbolic Coxeter
polytope in `H12`.

The subset census is

| doubled seed walls | seed chambers / volume multiple | facets |
|---|---:|---:|
| `31`, `58`, or `59` | 2 | 296 |
| `68` | 2 | 317 |
| any two of `31,58,59` | 4 | 523 |
| `68` and any one of `31,58,59` | 4 | 559 |
| `31,58,59` | 8 | 941 |
| `68` and any two of `31,58,59` | 8 | 1004 |
| `31,58,59,68` | 16 | 1827 |

Thus the exact facet-count histogram is

```text
296 x 3, 317 x 1, 523 x 3, 559 x 3,
941 x 1, 1004 x 3, 1827 x 1.
```

Distinct facet counts prove at least seven pairwise nonisomorphic types.
The single doublings account for only the `296`- and `317`-facet types, so
the cube contributes five further certified facet-count types:
`523`, `559`, `941`, `1004`, and `1827`.  This is a lower bound: the audit
does not identify chambers merely because they have the same facet count.

## Independence and source reconstruction

The hostile verifier does not import, invoke, or read the production cube
certifier or its archive.  It returns to the pinned GAP computational data
accompanying Dolgachev--Shimada and independently performs the following
reconstruction.

1. It parses `GramS16`, `ample16`, all 316 published `D16WallRecs`, and the
   trace-nine involution in the frozen fold-orbit record.
2. It verifies the involution, ambient lattice isometry, ample-class fixing,
   permutation of all 316 walls, and every one of its 180 wall orbits.
3. It reconstructs a saturated rank-13 fixed basis; an explicit maximal
   minor is a unit.  The resulting lattice has signature `(1,12)`,
   determinant `256`, and the independently rebuilt 172-wall raw Gram matrix
   agrees entry by entry with the frozen fold record.
4. It primitively normalizes and orients all 172 fixed-lattice wall roots by
   the published ample class.  These roots span rank 13.

The resulting four selected primitive roots have exact Gram matrix

\[
 \begin{pmatrix}
 -4&0&0&0\\
 0&-4&0&0\\
 0&0&-4&0\\
 0&0&0&-4
 \end{pmatrix}.
\]

Their complete seed relation rows are

| wall | `m2` | `m4` | parallel | ultraparallel | odd finite labels |
|---|---:|---:|---:|---:|---:|
| `31` | 46 | 20 | 45 | 60 | 0 |
| `58` | 46 | 20 | 45 | 60 | 0 |
| `59` | 46 | 20 | 45 | 60 | 0 |
| `68` | 25 | 4 | 59 | 83 | 0 |

Consequently every selected facet is a valid doubling wall, and they are
pairwise orthogonal.

## Faithful `C2^4` action and the gluing theorem

Let the four integral reflections be `s31`, `s58`, `s59`, and `s68`.  The
audit constructs their 13-by-13 lattice matrices directly.  They commute,
each changes the sign of its own selected normal and fixes the other three,
and all sixteen products are distinct.  This proves faithfulness without a
floating-point rank or word test.

There is also an exact Coxeter-group proof.  Send the four selected seed
generators to the standard basis of `C2^4` and every other seed generator to
zero.  The verifier checks every finite Coxeter relation: every relation
incident to one of the selected walls has even order, and the four selected
walls meet with order two.  Hence this assignment is a surjective Coxeter
homomorphism.

For a subset of `k` selected walls, its kernel has index `2^k`.  The external
wall reflections are precisely the Schreier conjugates of the nonselected
seed generators.  Their fundamental chamber is the convex union of the
`2^k` seed chambers indexed by that cube face.  This proves both that the
displayed boundary is the complete boundary and that its volume is exactly
`2^k` times the seed volume.

## Exact boundary checks

For each of the fifteen subsets, the verifier generates every external wall
root from scratch, removes only exact duplicate hyperplanes, and carries out
the following gates.

- Every boundary root is primitive, has negative square, and spans the full
  rank-13 Lorentzian space.
- Every root reflection is integral on the same full divisor lattice.
- The sum of the `2^k` transported seed interiors is timelike and pairs
  strictly positively with every exterior root.
- Each facet receives a transported timelike witness that vanishes on that
  facet and is strictly positive on every other boundary facet.  The minimum
  off-target slack is exactly `1` for all fifteen chambers.  Thus no listed
  inequality is redundant.
- Every unordered pair of facets is checked by exact integers.  Each pair is
  an allowed Coxeter angle, parallel, or ultraparallel; there are no negative
  oriented products and no unrecognized normalized products.
- The archive freezes independent SHA-256 digests of every ordered primitive
  boundary-root list, every strict-witness list, and every full exact pairing
  matrix.

For scale, the final `H12/1827` chamber has 1,668,051 unordered facet pairs.
Its exact relation inventory is

```text
m2             33,010
m3              1,740
m4             10,108
parallel       25,295
ultraparallel 1,597,898
```

These entries sum to 1,668,051, and all 1,827 reflections are integral.

## Finite volume versus infinite volume

The volume conclusion is not inferred from a small collection of parabolic
subdiagrams.  The seed `JK-H12/172` has two separately pinned exact gates:

- the complete literal Vinberg edge-endpoint and full parabolic-neighborhood
  audit has 910,316 elliptic rank-11 edges, no bad edge, and no invalid cusp
  neighborhood;
- the independent GMP cone enumeration has 113,274 finite vertices, 6,214
  ideal vertices, and zero hyperideal rays.

Both classify the seed as finite-volume noncompact.  The exact finite-index
kernel construction above therefore gives each descendant finite volume.
Because a finite-index subgroup of a nonuniform lattice remains nonuniform,
each descendant is noncompact as well.  This is the relevant distinction
from an infinite-volume chamber with some cusps.

## Bounded lineage and the worldwide caveat

The audit directly scans the pinned CoxIter snapshot at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe`: all 1,214 ordinary records and
all 2,959 records in its index-two doubling archive.

The ordinary `H12` headers have 13, 14, 16, or 18 facets; the doubling lane
has 18, 20, 25, or 35.  None of the seven new cube headers
`H12/296`, `H12/317`, `H12/523`, `H12/559`, `H12/941`, `H12/1004`, or
`H12/1827` occurs in either lane.

This is only grade-`B1` bounded-corpus evidence.  Those CoxIter lanes are not
a publication-complete denominator, so their zero headers do **not** prove
worldwide novelty.  The careful status is: seven exact types within this
construction, including five beyond its single doublings; worldwide
publication priority open.

## Algebro-geometric scope

All fifteen chambers use the same integral rank-13 divisor lattice extracted
from the general genus-2 Jacobian--Kummer chamber.  Thus they retain a weak
full-Picard divisor-reflection carrier.  Doubling across divisor mirrors does
not by itself identify the result with a full nef cone, and the lattice
reflections are not asserted to be surface automorphisms.

## Reproducibility

- [hostile verifier](../scripts/audit_jacobian_kummer_h12_doubling_cube_hostile_audit.py)
- [hostile certificate](../data/enumerations/jacobian_kummer_h12_doubling_cube_hostile_audit.json)
- [regression test](../tests/test_jacobian_kummer_h12_doubling_cube_hostile_audit.py)

Replay with

```bash
PYTHONPATH=.:src python3 -u \
  scripts/audit_jacobian_kummer_h12_doubling_cube_hostile_audit.py \
  --check-archive
PYTHONPATH=.:src pytest -q \
  tests/test_jacobian_kummer_h12_doubling_cube_hostile_audit.py
```

Frozen SHA-256 values:

```text
hostile verifier    b74c7295dc39895e25f1cc762ba8da1e430031acce64adcc825176a93cef9cf1
hostile certificate b0e94b5486e47fccc9586e10c3c67515dc5726e770295ab286d1394f2df62a33
regression test      d59d0c45a762e3133726e48cfd37db5bedceaa75d98233f155bff095788bcf1c
```

The regression suite passes `4/4` tests.
