# TOOT: all ten rank-nine Roulleau chambers are exact full-nef finite-volume polytopes

## Outcome

All ten Picard-rank-nine rows in Roulleau's K3 atlas have now been expanded
into exact Coxeter-polytope certificates.  In every row the listed walls cut
out the **entire nef cone**, not merely a divisor subsystem, and Vinberg's
literal endpoint/parabolic criterion proves finite volume and noncompactness.

| ID | `NS(X)` | polytope | finite vertices | cusps | elliptic edges | lineage |
|---|---|---:|---:|---:|---:|---|
| `R9-UE7-NEF-H8-9` | `U+E7(-1)` | `H8/9` | 8 | 1 | 36 | JKRT `T8` |
| `R9-UD6A1-NEF-H8-10` | `U+D6(-1)+A1(-1)` | `H8/10` | 12 | 3 | 63 | Tumarkin `12-k=2` |
| `R9-UD4A1x3-NEF-H8-15` | `U+D4(-1)+A1(-1)^3` | `H8/15` | 53 | 11 | 315 | no same header in pinned corpus |
| `R9-UA1x7-NEF-H8-37` | `U+A1(-1)^7` | `H8/37` | 569 | 78 | 4,188 | no same header in pinned corpus |
| `R9-U2A1x7-NEF-H8-240` | `U(2)+A1(-1)^7` | `H8/240` | 17,280 | 2,160 | 207,360 | no same header in pinned corpus |
| `R9-UA7-NEF-H8-9` | `U+A7(-1)` | `H8/9` | 7 | 2 | 36 | JKRT `P8` |
| `R9-UD4A3-NEF-H8-10` | `U+D4(-1)+A3(-1)` | `H8/10` | 17 | 4 | 90 | Tumarkin `09` |
| `R9-UD5A2-NEF-H8-10` | `U+D5(-1)+A2(-1)` | `H8/10` | 18 | 1 | 81 | Tumarkin `05` |
| `R9-UD7-NEF-H8-9` | `U+D7(-1)` | `H8/9` | 8 | 1 | 36 | JKRT `Q8` |
| `R9-UA1E6-NEF-H8-10` | `U+A1(-1)+E6(-1)` | `H8/10` | 14 | 1 | 63 | Tumarkin `13-k=2` |

This closes the finite published rank-nine denominator positively.  It is a
substantial new project-level surface/polytope alignment and volume audit,
but the surfaces and curve configurations were published by Roulleau and
seven abstract diagrams are classical named polytopes.  No worldwide-new
abstract-polytope claim is made for the other three.

The replay is
[`certify_roulleau_rank9_all_full_nef.py`](../scripts/certify_roulleau_rank9_all_full_nef.py),
with frozen data in
[`roulleau_rank9_all_full_nef.json`](../data/enumerations/roulleau_rank9_all_full_nef.json).

## The 240-wall anchor: degree-one del Pezzo geometry

The largest row is especially clean algebraically.  Start with the 240
`(-1)`-curves of a degree-one del Pezzo surface.  In the standard basis
`L,E1,...,E8`, the certificate enumerates them in the classical seven
families

```text
8 exceptionals
28 lines
56 conics
56 cubics
56 quartics
28 quintics
8 sextics
```

Their pullbacks to the K3 surface are the complete 240 irreducible
`(-2)`-curves.  The lattice form is

```text
diag(2,-2,-2,-2,-2,-2,-2,-2,-2),
```

and the anticanonical vector `(3,-1,...,-1)` has square `2` and intersection
`2` with every wall.  Thus one uniform algebro-geometric construction
produces the whole Coxeter matrix, an explicit strict ample orientation, and
full-nef equality.  This is the sort of structured positive object sought by
the surface-to-polytope track, rather than an unconstrained diagram search.

The literal polytope audit finds `17,280` finite vertices and `2,160` ideal
vertices.  Its `207,360` elliptic edges split as

```text
 69,120 with 0 finite + 2 ideal endpoints
138,240 with 1 finite + 1 ideal endpoint
```

There are no missing endpoints.  The `15,120` connected affine components
consolidate into the `2,160` full cusps, and no affine zero-neighborhood
fails the parabolic-rank test.

## Why every row is the full nef cone

Roulleau supplies a concrete K3 construction, the curve configuration, and
the complete number of irreducible `(-2)`-curves in each rank-nine case.  The
certifier reconstructs every wall Gram matrix and checks an integral
rank-nine basis, inertia `(1,8,0)`, full wall span, and a strictly positive
ample class.  Since the archived wall set is the complete irreducible
`(-2)`-curve set, the K3 ample-cone theorem identifies its positive-cone
chamber with `Nef(X)` itself.

The geometric models include double covers of `F4` with prescribed ADE
branch singularities, double planes with controlled singular sextics, and
the degree-one del Pezzo double cover above.  Hence the walls are not merely
formal lattice roots: they are geometrically produced divisor classes on
the surfaces.

## The finite-volume distinction is literal

For an eight-dimensional Coxeter polytope the implemented Vinberg gate checks
both requirements:

1. every elliptic rank-seven edge diagram has exactly two finite or ideal
   endpoint extensions; and
2. every connected affine diagram has a full zero-neighborhood of parabolic
   rank seven.

Both bad counts are zero in all ten rows.  Every row has at least one full
parabolic cusp.  The conclusion is therefore **finite-volume noncompact** in
each case—not compactness and not infinite volume.

The exact verifier is pinned by hash.  Its original 128-facet bit container
is deterministically extended to four 64-bit words for the 240-wall row;
the finite/affine classification and Vinberg logic are unchanged.

## Lineage and novelty boundary

Seven diagrams have literal relation-graph matches, with exact facet
permutations and source-file hashes frozen in the certificate:

| K3 row | exact CoxIter record |
|---|---|
| `U+E7(-1)` | JKRT `T8` |
| `U+D6(-1)+A1(-1)` | Tumarkin `8-tum04_11_12-k=2` |
| `U+A7(-1)` | JKRT `P8` |
| `U+D4(-1)+A3(-1)` | Tumarkin `8-tum04_11_09` |
| `U+D5(-1)+A2(-1)` | Tumarkin `8-tum04_11_05` |
| `U+D7(-1)` | JKRT `Q8` |
| `U+A1(-1)+E6(-1)` | Tumarkin `8-tum04_11_13-k=2` |

The pinned 1,214-record ordinary CoxIter snapshot and 2,959-record doubling
snapshot contain no same-`H8`/facet header for the 15-, 37-, or 240-wall
rows.  This proves only bounded-corpus absence.  Because Roulleau already
published the complete configurations, it is not promoted to a worldwide
novelty statement.

## A source erratum

In the `U+D7(-1)` case, the sentence describing the square-260 ample class
assigns degree `2` to `A2` and then also includes `A2` among the degree-one
walls.  The printed graph and exact lattice computation resolve the clash:

```text
(D260.A1,...,D260.A9) = (1,1,1,2,1,1,1,1,1).
```

Thus degree `2` is on `A4`; the ample class remains strict and the polytope
conclusion is unaffected.

## Replay

```bash
python3 -m scripts.certify_roulleau_rank9_all_full_nef --check-archive
python3 -m pytest -q tests/test_roulleau_rank9_all_full_nef.py
```
