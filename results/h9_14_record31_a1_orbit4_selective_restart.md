# Record 31, `a=1`, source orbit 4: partitioned incidence no-go

## Current status

Source orbit 4 (positive-`A4` mask 10252) is now **incidence-impossible** in
the corrected compact `a=1` schedule.  The exhaustive proof is partitioned
across all six canonical depth-one classes, and both pieces have empty
canonical `{2:540}` target-orbit sets.

The first inline production traversal completed root classes 1, 2, and 3
before entering class 4.  Its first post-boundary checkpoint was

```text
79,000,000 nodes
seen = [1,4,51,1749,87940,4760678,49082526,22912]
ridge-valid target orbits = 0.
```

The later 83-million-node class-4 checkpoint is diagnostic only.  SIGINT
left the redirected stdout at zero bytes, so there is no misleading terminal
JSON.  A separately terminal traversal then selected exactly original root
classes 4, 5, and 6.  It closed as follows:

| terminal complementary datum | value |
|---|---:|
| search nodes | 105,066,986 |
| discovered / selected / skipped root classes | 6 / 3 / 3 |
| pure 120-facet targets | 3,184 |
| distinct pure ridge signatures | 818 |
| `{2:540}` target orbits | **0** |

Thus the theorem-facing canonical target-orbit union is `empty ∪ empty` over
classes `1--3` and `4--6`.  The boundary and terminal metadata are frozen in
`data/enumerations/h9_14_record31_a1_orbit4_inline_root_boundaries.json`,
SHA-256
`5644ed22f418c767849a2e9a9c4f57f0b2088e2d2d6aa4eb9a4c66a044c3c209`;
the complementary transcript is
`data/enumerations/h9_14_record31_a1_orbit4_root_classes4_6.json`, SHA-256
`dc463688433307ad8fe04c0a03af73b2080ef0fa49e71732743eaaca80675acc`.

## Selective classes 4--6 restart

`scripts/generate_h9_14_record31_a1_root_class_engine.py` adds the inclusive
option

```text
--quotient --root-classes BEGIN END
```

to the fixed-width inline engine.  It canonicalizes every raw depth-one
branch in the original deterministic order.  A separate root-code set records
skipped classes, so a later raw representative of a skipped class cannot be
misidentified as a new selected class.  The completed complementary command was

```bash
printf '1\n4 10252\n' | \
  caffeinate -i /private/tmp/record31_a1_root_class \
  --quotient --root-classes 4 6 \
  > /private/tmp/record31_a1_orbit4_root_classes4_6.json
```

The generated source SHA-256 is
`b2e34b810432db09fe685b52b39ff97ade0c7c1b3250deea0fbe8904713ba945`;
the currently compiled binary SHA-256 is
`c16296e21f5eddc43911c2b1bcbb5f05bccd8850911c2b057390a53a9d10ee0e`.

The complementary run restarted the deeper quotient memo.  Consequently its
node count, pure-target multiplicities, and internal histograms can overcount
states also reachable from classes 1--3 and must not be added to the
interrupted counts.  Its canonical ridge-valid target-orbit set exhaustively
overcovers classes 4--6.  Because it is empty, no topology, Coxeter-label, or
Gram gate is needed for source 4.

## Exact bounded validation

On complete source orbit 0, the unrestricted range-enabled engine reproduces
the old 231,773-node JSON transcript exactly after removing only its six new
range-metadata fields.  Splitting its four canonical root classes as `1--2`
and `3--4` gives:

| run | nodes | selected classes | target orbits |
|---|---:|---:|---:|
| unrestricted | 231,773 | 4 | 1 |
| classes 1--2 | 66,589 | 2 | 0 |
| classes 3--4 | 199,845 | 2 | 1 |

The fragment target-orbit union equals the unrestricted orbit set exactly.
The split node total 266,434 is deliberately larger than 231,773, providing a
frozen regression witness that fragment internal counts are nonadditive.

The three transcripts, their byte hashes, the generated-source hash, the
source-4 boundary, and the target-orbit union are checked by
`scripts/verify_h9_14_record31_a1_root_class_split.py` and
`tests/test_h9_14_record31_a1_root_class_split.py`.
The terminal source theorem is separately pinned by
`scripts/verify_h9_14_record31_a1_orbit4_partitioned_transcript.py` and
`tests/test_h9_14_record31_a1_orbit4_partitioned_transcript.py`.

This is exactly a source-orbit-4 compact-incidence no-go.  Production source
2 and independent-incidence gaps for sources 2, 3, and 4 remain, so it is not
the full `a=1` schedule theorem.  It makes no finite-volume noncompact
inference.
