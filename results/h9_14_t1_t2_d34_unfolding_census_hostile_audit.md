# Hostile audit: the bounded `T1/T2 D34` unfolding census

Status: **independent exact audit PASS** (2026-08-08).

## Verdict

No discrepancy was found.  An implementation that imports no production
compiler code and uses no Hermite normal form reproduces the entire claimed
bounded census:

| item | independent result |
|---|---:|
| pinned primary corpus | 6 markings, 4 diagram types |
| `D34` filter | 3 markings, 2 diagram types |
| primitive crystallographic scaling rays | `18+18=36` |
| minimal one/two-sheet-compatible rays | `2+2=4` |
| reflection-stable even fixed overclosures | `21+7=28` |
| weak-compiler survivors | exactly one for `T1`, exactly one for `T2` |

The coordinate-free invariant multiset of every one of the 28 overclosures
agrees with the production archive.  Both survivors have fixed lattice
`U+E8(-2)`, by signature and two-elementary invariants `(10,8,0)`.

## Independent route

The audit deliberately replaces the production arithmetic route.

1. It computes a full Smith decomposition of each degenerate `14x14` facet
   Gram.  The last four Smith columns give the saturated integral radical;
   the first ten give a quotient basis.  It verifies the unimodular
   decomposition and exact labelled identity `Q^T F Q=G`.  No ten facets are
   assumed to be an integral basis.
2. It enumerates scaling rays as integer height functions on the graph of
   `m3`-components, with adjacent `m4` components differing by one.  This is
   independent of the production edge-orientation enumeration.
3. It enumerates reflection-stable totally isotropic submodules of the
   two-elementary discriminant space by reflection-orbit closure.  Each
   overclosure is built with a reduced `GF(2)` pivot basis, not HNF.
4. It searches embeddings of the complete marked long-root class set into
   `E8/2E8`, then builds the rank-18 and K3 graph glues with the same pivot
   construction.  Root-span indices are computed by rectangular Smith form
   of the `18x24` lifted-root matrix.

## The delicate `T1` arithmetic

The independent calculation confirms both index-two phenomena; they are not
silently conflated with index one.

- The compatible `short4_long10` folded root lattice has determinant
  `-1024`.  Its unique surviving fixed overclosure has index `2`, determinant
  `-256`, Smith diagonal `1^2 2^8`, discriminant histogram
  `136*0+120*1`, and therefore delta `0`.
- Its ten long classes span dimension `7` with radical dimension `1`.  The
  graph glue with `E8(-2)` has dimension `7`, hence index `2^7`; the resulting
  rank-18 lattice has determinant `-4`, Smith diagonal `1^16 2^2`, and exact
  invariants `(18,2,0)`.  Its discriminant quadratic histogram is
  `3*0+1*1`, again proving delta `0` directly.
- The 24 square-minus-two roots have off-diagonal histogram
  `242*0+34*1`.  Their rectangular Smith diagonal is `1^17 2`, so they span
  index **two** in the rank-18 lattice.  They form four fixed roots and ten
  exchanged pairs, are acute, and their orbit sums reproduce the labelled
  `T1` Gram exactly.
- An anti-isometry of all four discriminant classes with `U+U(2)` gives an
  index-four graph glue.  The resulting rank-22 lattice is even, unimodular,
  and has signature `(3,19)`; both summands embed primitively.  The sheet
  involution extends with characteristic polynomial
  `(lambda-1)^14*(lambda+1)^8`.

For comparison, the independent `T2` replay gives `(18,4,0)`, determinant
`-16`, root-span index `1`, pairing histogram `244*0+32*1`, and primitive K3
complement `U(2)+U(2)`.

## Exact boundary

- `U` remains open.  This is exhaustive only relative to the pinned primary
  exact-Gram corpus; it is not a census of every abstract `H9/14` diagram and
  does not prove an algebraic-surface existence locus.
- `R0` is independently closed for these two labelled data: lattice/action
  construction and exact fold-back both pass.
- `C` is not promoted here.  In particular, `T1` cone completeness and
  hidden-wall exclusion remain open.
- `V` is not recomputed in this audit.  The stored finite-volume Vinberg
  result is imported only after exact labelled diagram identity.
- This supplies a new weak carrier lineage for the known `T1` polytope, not a
  new abstract hyperbolic polytope.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/audit_h9_14_t1_t2_d34_unfolding_census_hostile.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_t1_t2_d34_unfolding_census_hostile_audit.py
```

The focused suite has **5 passing tests**.  Frozen archive SHA-256:
`e637a48b78d203655c4f8db092e3029e1257408c1bad23c27ddfb32e2ae54f24`.

Artifacts:

- [`audit_h9_14_t1_t2_d34_unfolding_census_hostile.py`](../scripts/audit_h9_14_t1_t2_d34_unfolding_census_hostile.py)
- [`h9_14_t1_t2_d34_unfolding_census_hostile_audit.json`](../data/enumerations/h9_14_t1_t2_d34_unfolding_census_hostile_audit.json)
- [`test_h9_14_t1_t2_d34_unfolding_census_hostile_audit.py`](../tests/test_h9_14_t1_t2_d34_unfolding_census_hostile_audit.py)
