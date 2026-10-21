# Hostile audit: low-`B` star--`P5` exact separation and the fork-`B4` prior

**Status (2026-08-10): PASS/FREEZE.**  This was a read-only, lightweight
audit.  The 22,463-row production comparison and the active 106,385-row
fork-`B4` census were not replayed.

## Exact-comparison audit

The low-`B` certifier rebuilds every selected boundary from the colour word
streamed out of its pinned source archive.  Its ordinary invariant is the
lexicographically sorted multiset of the five ordinary Coxeter-label degrees
at every facet, serialized as fixed-width little-endian `uint16` rows.  An
ordinary Coxeter-matrix isomorphism preserves this byte word.

The theorem-producing comparison is lossless: all 1,770 rebuilt star byte
words are retained, and each of the 20,693 rebuilt `P5` byte words is tested
by literal Python `bytes` membership inside the matching `B`/facet header.
The frozen result is zero literal matches.  SHA-256 is used only to pin files,
link each rebuilt row to its archived digest, and summarize streams; it is not
used in place of cross-shape byte equality.  In particular, no collision
resistance assumption is needed for the zero-match conclusion.

All ten pinned input hashes match the files currently read by the certifier,
and the archive's builder hash matches the current builder.  The focused
tests pass (`3 passed`).  No provenance or literal-comparison defect was
found.

Frozen hashes checked in this audit:

```text
builder  5c13641fca1eef719c917bd1c1b6b62a97176463b61c6a3c1056b2e25af7112e
archive  b73e85b9cd3fa61fc727724dca966e906af33feadcaa5b656c758f6640dee15f
report   de43769144b3dffc6b79d799b532f301e1ac782545e10752ac026e755a852092
tests    a2c29aa683eb99e0cc92ca3cb237c8b7434709602dce2c6a2138ee3838a919cb
```

## Active fork-`B4` prior iterator

Static inspection of
`certify_jacobian_kummer_h12_five_vertex_fork_streaming_census.py` shows that
the `B=4` prior iterator is exhaustive for the only matching facet header:

1. it regenerates all `108,159` `P5/BBBB` orbit representatives, requires
   that denominator exactly, and aligns them by `zip(..., strict=True)` with
   all `108,159` archived degree digests from the pinned gzip transcript;
2. it regenerates all `16,872` star representatives, aligns them strictly
   with the pinned complete-star archive, and retains exactly the `8,985`
   rows having four `B` colours.

Thus the active fork-`B4` comparison has exactly

```text
108,159 P5/BBBB + 8,985 star/B4 = 117,144 prior rows.
```

The active certifier uses archived digests only as a lossless negative filter:
unequal digests imply unequal bytes, while every equal-digest candidate is
rebuilt and compared literally; every literal match is then subjected to the
exact canonical colour-refinement invariant.  Relevant current and pinned
hashes checked here include

```text
P5/BBBB rows  d367b7d208ad8c7fff294978a1857f7557d3755f33df8e4afff6689ad2aec1a7
star archive  21fbd0946b910360cdd4096392e4e629907f971be85359d5b87e312a78d4f309
fork builder  64dd51e4410a379723669eee4002f16029c4cab296ba8693e4483f1eb7f6627e
```

## Union arithmetic and claim boundary

The already proved union is

```text
201,747 complete P5 + 16,872 complete star + 92,977 fork B<=3 = 311,596.
```

The frozen fork denominator gives `106,385` rows at `B=4`.  Therefore, if
the active census proves every one new relative to the `117,144` same-header
priors and internally pairwise nonisomorphic, the complete conditional union
will be

```text
311,596 + 106,385 = 417,981.
```

The number `417,981` is conditional until that active run finishes and its
archive, deterministic transcript, and launch-time builder provenance are
frozen and checked.  Worldwide priority, full-nef effectivity, and surface
action remain open and are not enlarged by this audit.
