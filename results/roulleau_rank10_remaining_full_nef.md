# TOOT: all seven remaining rank-ten Roulleau chambers are now exact full-nef finite-volume polytopes

## Outcome

The remaining seven rows in Roulleau's rank-ten K3 table have been expanded
from published curve pictures into exact hyperbolic-polytope certificates.
Every row is the **entire nef cone** of the stated K3 surface, and every row
passes Vinberg's literal finite-volume test:

| ID | `NS(X)` | polytope | finite vertices | cusps | elliptic edges | verdict |
|---|---|---:|---:|---:|---:|---|
| `R10-UE7A1-NEF-H9-11` | `U+E7(-1)+A1(-1)` | `H9/11` | 15 | 2 | 80 | finite-volume, noncompact |
| `R10-UD4D4-NEF-H9-11` | `U+D4(-1)^2` | `H9/11` | 16 | 10 | 125 | finite-volume, noncompact |
| `R10-UD6A1x2-NEF-H9-14` | `U+D6(-1)+A1(-1)^2` | `H9/14` | 37 | 7 | 231 | finite-volume, noncompact |
| `R10-U2D4D4-NEF-H9-18` | `U(2)+D4(-1)^2` | `H9/18` | 186 | 71 | 1,840 | finite-volume, noncompact |
| `R10-UD4A1x4-NEF-H9-27` | `U+D4(-1)+A1(-1)^4` | `H9/27` | 313 | 45 | 2,245 | finite-volume, noncompact |
| `R10-UA1x8-NEF-H9-145` | `U+A1(-1)^8` | `H9/145` | 13,824 | 1,249 | 123,520 | finite-volume, noncompact |
| `R10-UA2E6-NEF-H9-11` | `U+A2(-1)+E6(-1)` | `H9/11` | 17 | 5 | 105 | finite-volume, noncompact |

Together with the separately certified `U+E8(-1)` and `U+D8(-1)`
simplices, this closes all nine published Picard-rank-ten rows positively.
The result is a new dedicated project-level AG/polytope alignment and volume
certificate.  The surfaces and curve configurations are published, so this
is not presented as seven worldwide-new abstract polytopes.

The replay is
[`certify_roulleau_rank10_remaining_full_nef.py`](../scripts/certify_roulleau_rank10_remaining_full_nef.py),
with frozen data in
[`roulleau_rank10_remaining_full_nef.json`](../data/enumerations/roulleau_rank10_remaining_full_nef.json).

## Why these are full nef cones

[Roulleau's atlas](https://arxiv.org/abs/2003.08985) gives both the curve
configurations in its rank-ten section and the complete numbers of
irreducible `(-2)`-curves in its final table:

```text
11, 11, 14, 18, 27, 145, 11.
```

The certifier reconstructs every printed curve class and Gram matrix.  In
each case ten named curves form an integral basis with inertia `(1,9,0)` and
the expected determinant and Smith invariants.  The printed ample class is
strictly positive on every wall.  Because the list is the complete set of
irreducible `(-2)`-curves, the K3 ample-cone theorem identifies the chamber
cut out by those walls with `Nef(X)` itself.  These are therefore much
stronger than weak divisor-lattice carriers or unproved wall subsystems.

The constructions are concrete:

- `E7+A1`, `D4^2`, and `A2+E6` are double covers of `F4` with the displayed
  ADE singularities in the branch curve;
- `D6+2A1` is a double plane branched over a line and a nodal quintic;
- `U(2)+D4^2` is a double cover of `F2` with eight even-contact fibers;
- `D4+4A1` is a double plane branched over two lines and a quartic;
- `8A1` is a double plane branched over a nodal sextic with eight
  even-contact lines and sixty-four six-tangent conics.

## The finite-volume test is literal

For a Coxeter polytope in `H9`, Vinberg's criterion used here has two parts.

1. Every elliptic rank-eight edge diagram must have exactly two endpoints,
   each finite or ideal.
2. Every connected affine diagram must have a full zero-neighborhood whose
   parabolic rank is eight.

Both bad counts are zero in every row.  Each row has cusps, so the conclusion
is finite volume and noncompactness—not compactness, and not infinite volume.
For example, the 145-wall chamber has endpoint histogram

```text
28,800 edges with 0 finite + 2 ideal endpoints
65,024 edges with 1 finite + 1 ideal endpoint
29,696 edges with 2 finite + 0 ideal endpoints
```

whose sum is all `123,520` elliptic edges.  It also has `6,632` connected
affine components, which consolidate into `1,249` full cusps, with zero bad
parabolic neighborhoods.  The old exact verifier used two 64-bit wall masks;
the replay deterministically expands only that container to three words and
leaves its finite/affine classification and Vinberg logic unchanged.

## Lineage

Three of the wall diagrams have unique exact same-header matches in the
pinned 1,214-record ordinary CoxIter snapshot:

| K3 row | exact classical record |
|---|---|
| `E7+A1`, `H9/11` | Tumarkin `9-tum04_11_09-k=2` |
| `D4^2`, `H9/11` | Tumarkin `9-tum04_11_06` |
| `A2+E6`, `H9/11` | Tumarkin `9-tum04_11_07` |

The certificate freezes the exact wall permutations and source-file hashes.
The ordinary CoxIter snapshot has no `H9/14`, `H9/18`, `H9/27`, or
`H9/145` header.  The pinned 2,959-record doubling lane has eight `H9/27`
headers, but none has even the same relation histogram; it has no same-header
record for the other three sizes.  This is only a bounded corpus statement.
The larger four diagrams were already printed as complete K3 curve
configurations by Roulleau, so corpus absence is not converted into a
worldwide novelty claim.

## A small source erratum

In the `U+D4(-1)+A1(-1)^4` row, the printed graph and the printed coordinate
vector

```text
D_26 = (3,3,3,3,7,3,1,0,0,0)
```

give `D_26.A6=2` and `D_26.A5=1`.  The following prose says the reverse.
The square remains `26`, and every wall degree is positive, so this label
swap does not affect ampleness or the polytope theorem.  The frozen
certificate records the literal matrix calculation.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_roulleau_rank10_remaining_full_nef.py --check-archive
PYTHONPATH=.:src python3 -m pytest -q tests/test_roulleau_rank10_remaining_full_nef.py
```
