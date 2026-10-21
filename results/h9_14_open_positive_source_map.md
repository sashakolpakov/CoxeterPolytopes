# Positive-completion source map for the seven open H9/14 shadows

Status: exploratory exact census, 2026-08-07.  This note selects the next
whole-shadow compact targets.  It is not yet a frozen all-source theorem and
makes no finite-volume noncompact inference.

## Exact source-orbit map

For each open graph shadow, the script
`scripts/explore_h9_14_open_positive_a4.py` enumerates every labelled family
of one or two positive minimal four-nonfaces.  It applies the complete local
vertex/edge/triple intervals, the forced five-shadow, strict link positivity,
and aggregate Macaulay capacity, then quotients by the full missing-graph
automorphism group, including permutations of isolated vertices and of the
new rows.

| record | viable labelled `a=1` | `a=1` source orbits | viable labelled `a=2` | `a=2` source orbits |
|---:|---:|---:|---:|---:|
| 0 | 880 | 4 | 384,780 | 36 |
| 1 | 505 | 3 | 127,260 | 23 |
| 3 | 870 | 3 | 376,575 | 35 |
| 4 | 825 | 5 | 335,580 | 64 |
| 5 | 660 | 6 | 173,031 | 69 |
| 7 | 450 | 7 | 98,361 | 77 |
| 9 | 450 | 3 | 97,065 | 25 |

Thus Record 1, `a=1`, is the first exact completion target, tied with Records
3 and 9 at three source orbits but having the smallest `a=2` continuation.

The optimized `a=2` interval gate is exact: every enumerated four-face already
passes all local intervals by itself, all interval lower endpoints are zero,
and for a two-row family the only newly binding constraints have upper
capacity one.  The implementation compiles precisely those forbidden pairs.

## First Record-1 source

The first `a=1` representative is `{0,3,4,5}`.  Its exact fixed-source model
has:

- 9 forced nonminimal five-faces and 19 additional minimal five-nonfaces;
- 1,278 candidate additions;
- 259 base transversals of order at most four;
- 634 base order-five transversals;
- 1,281 six-set ridge equations; and
- scheduled face-vector endpoint 158.

### Correction of the exploratory target count

The earlier command was run with `--max-models 10`; it left no archived
terminal `UNSAT` transcript.  The previous wording that it gave an exact
nine-target census was false.  In fact, the first frozen target alone has
60,480 distinct labelled images under the fixed-source stabilizer

```text
S2({1,2}) x S3({3,4,5}) x S7({6,...,12}).
```

The graph/source signatures prove that this is the full stabilizer, of order
60,480, and exhaustive streaming over it finds that the representative target
has trivial stabilizer.  Thus orbit--stabilizer gives 60,480 targets before
any other target orbit is considered.  A cutoff, timeout, kill, or missing
terminal status is never interpreted as `UNSAT`.

The first literal target is already closed at the label layer.  Production
and a no-production-classifier replay both reconstruct 158 actual rank-nine
vertices and give the same complete compact-label survivor sequence

```text
108 -> 84 -> 8 -> 0.
```

The independent state hashes are:

```text
857c6268544a6e22216a4068a1642a2e48730ff5288ae810cdaa398a9517f730
de2dcf81429b02fe4970c82cc168428caddf019335da066bed4c228047e8b2e5
a3cd5a07b59bf0c8e0b299bd334049f235dc243abded9be5f31d1ad6f011f270
4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945
```

This label obstruction is invariant under relabelling.  Therefore it excludes
the full 60,480-member labelled orbit of this target, not merely the displayed
representative.  Production and separately derived hostile scripts agree on
the action-stream hash, the trivial target stabilizer, all 158 rank-nine
vertices, the ridge histogram `{0:2292, 2:711}`, and the label-state sequence.
This is one target orbit only: it does not complete the symmetry-quotient
incidence census, exclude another target orbit, close the other two `a=1`
source orbits, close another schedule, or say anything about finite-volume
noncompact polytopes.

## Completion path

1. Run a symmetry-quotiented production incidence census for source orbit 0
   through a fresh terminal `UNSAT`, with streaming checkpoints.
2. Replay the quotient target list with an independently derived CP-SAT or
   lossless combinatorial model and compare canonical orbit keys literally.
3. Run production and hostile compact-label audits on every resulting target
   orbit; the frozen 60,480-member orbit is already closed.
4. Repeat for source orbits 1 and 2.
5. If all three sources die, issue a complete Record-1 `a=1` certificate; if a
   label state survives, pass it immediately to exact Gram/signature and
   Vinberg compactness gates.

Heavy completion solves are temporarily paused to protect the live
Record-31 search with more than 370 million exact nodes.  The finite target
audits and code hardening remain safe light work.

## Replay

```text
PYTHONPATH=src python3 scripts/explore_h9_14_open_positive_a4.py \
  0 1 3 4 5 7 9 --max-a 2 --summary
PYTHONPATH=src:. python3 scripts/explore_h9_14_open_positive_completion.py \
  1 1 0 --max-models 10
PYTHONPATH=src:. python3 \
  scripts/independent_h9_14_record1_a1_target0_label_audit.py
PYTHONPATH=src:. python3 \
  scripts/certify_h9_14_record1_a1_orbit0_target_orbit.py
PYTHONPATH=src:. python3 \
  scripts/independent_h9_14_record1_a1_orbit0_target_orbit_audit.py
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_h9_14_open_positive_completion.py \
  tests/test_h9_14_record1_a1_orbit0_target_orbit_no_go.py
```
