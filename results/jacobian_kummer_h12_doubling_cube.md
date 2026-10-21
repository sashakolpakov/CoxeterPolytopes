# A four-wall K3-carried doubling cube in `H12`

## Result

The certified Jacobian--Kummer `H12/172` chamber has four mutually
orthogonal doubling facets, numbered

```text
31, 58, 59, 68.
```

Their integral reflections commute and act faithfully as `C2^4`.  Every one
of the 15 nonempty selected-wall subsets therefore gives an exact convex
Coxeter union of `2^k` seed chambers.  The complete marked census is:

| selected walls | marked subsets | facets | volume multiplier |
|---:|---:|---:|---:|
| 1 | 4 | `296` three times, `317` once | 2 |
| 2 | 6 | `523` three times, `559` three times | 4 |
| 3 | 4 | `941` once, `1004` three times | 8 |
| 4 | 1 | `1827` once | 16 |

The seven distinct facet counts prove at least seven pairwise nonisomorphic
diagram types without any graph-canonicalization assumption.  The counts
`523, 559, 941, 1004, 1827` add five types beyond the single-doubling
theorem.  Same-count subsets may split further; no equality of their labelled
diagrams is asserted here.

The production artifacts are the
[`certificate builder`](../scripts/certify_jacobian_kummer_h12_doubling_cube.py),
[`JSON archive`](../data/enumerations/jacobian_kummer_h12_doubling_cube.json),
and [regression test](../tests/test_jacobian_kummer_h12_doubling_cube.py).

## Exact geometric gates

The script starts from the pinned rank-13 full-Picard carrier

\[
 U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle
\]

and its 172 primitive oriented wall vectors.  It verifies that all four
selected reflections are integral, that every selected pair is orthogonal,
and that no selected wall has an odd incident finite Coxeter order.  Hence
each remaining selected wall survives every prior doubling and the ordinary
even-angle doubling lemma applies iteratively.

For each subset, all `2^k` lattice transforms are applied to every seed wall.
Internal selected mirrors are removed and equal unoriented boundary roots
are identified exactly.  The sum of the transformed seed interiors is
strictly positive on every retained root.  Each boundary spans the full
13-dimensional linear space, every root is primitive, and every reflection
is integral on the same Picard lattice.

For one representative of every distinct facet count, the audit goes further.
If `q_r` is the exact timelike seed witness for wall `r`, each cube copy uses
its transformed witness.  Every one of the 296, 317, 523, 559, 941, 1004, or
1827 witnesses is replayed against every boundary inequality.  The minimum
off-target slack is one and the minimum witness square is `17/2` in every
representative.  Thus the largest 1,827-wall list is irredundant, not merely
a raw orbit count.

## Volume, compactness, and surface meaning

The seed has an independent GMP extreme-ray enumeration and literal Vinberg
edge/cusp certificate.  Iterated convex doubling gives disjoint seed-chamber
interiors, so each descendant has exactly the displayed finite volume
multiple.  The descendant group has finite index in the nonuniform seed
reflection group, hence every descendant is nonuniform and its polytope is
noncompact.

This gluing argument is the exact volume certificate; a fresh enumeration of
the enormous descendant vertex sets is neither performed nor needed.  All
reflected boundary classes stay in the same full Néron--Severi lattice, so
the same seven-dimensional very-general K3 family is a weak
divisor-reflection carrier for the entire cube.  Doubling crosses walls: no
full-nef, effective-curve, or surface-automorphism statement follows.

## Novelty boundary

The seven facet counts are new to the project cube ledger, and five are new
beyond the separately frozen single doubles.  The pinned ordinary and
index-two CoxIter `H12` headers stop at 18 and 35 facets respectively, so all
seven are absent by header from those snapshots.  This is bounded evidence,
not a worldwide priority theorem.  Coxeter doubling itself is classical; the
positive project content is the Jacobian--Kummer/K3 seed, the exact four-wall
configuration, the complete subset census, and the irredundant large
boundary systems.

## Independent hostile replay

The separate
[`no-import cube audit`](jacobian_kummer_h12_doubling_cube_hostile_audit.md)
starts from the pinned Dolgachev--Shimada raw data rather than this builder or
archive.  It reconstructs the saturated rank-13 lattice and all 172 walls,
proves that the four selected roots have Gram matrix `-4*I4`, and rebuilds
all 15 boundary systems.  Every boundary has full rank, primitive integral
reflections, strict timelike witnesses with minimum off-target slack one, and
only Coxeter, parallel, or ultraparallel pair relations.  It independently
recovers the exact facet histogram and scans all 1,214 ordinary plus 2,959
doubling records in the pinned CoxIter snapshots.  Thus the internal theorem
is independently frozen; the literature-priority caveat remains unchanged.

## Replay

```bash
PYTHONPATH=.:src python3 -u scripts/certify_jacobian_kummer_h12_doubling_cube.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_doubling_cube.py
PYTHONPATH=.:src python3 -u scripts/audit_jacobian_kummer_h12_doubling_cube_hostile_audit.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_jacobian_kummer_h12_doubling_cube_hostile_audit.py
```

The production script SHA-256 is
`c09f9b326ca1d1df9883a0d5fd8c83d312ba5803fe72a9cf86144adb197c94ed`;
the archive SHA-256 is
`0267e3d56b5b0acaa7015e4cadac1d357d387532831e89485c1846eae28be3f1`.
