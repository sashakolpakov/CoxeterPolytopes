# Record 31, `a=1`, source orbit 2: first matrix-root tranche

## Exact outcome

The first of the six canonical depth-one classes in the independent direct
binary-matrix quotient is terminal.  It contains eight nonisomorphic
incidence targets with the required pure 120-facet blocker and ridge
histogram `{2:540}`.  All eight pass the connected mod-two homology-sphere
audit, and both unrestricted compact Coxeter-label joins leave zero states on
every target.

This excludes **root class 1 only**.  Source orbit 2 and the full Record-31
`a=1` schedule remain open until root classes 2 through 6 are exhausted.  No
timeout, memory-limit disposition, or finite-volume noncompact inference is
used.

## Terminal incidence data

The source is the positive A4 `{2,3,11,12}` (mask `6156`).  The complete
root-class fragment has:

| datum | exact value |
|---|---:|
| search nodes | 11,053,340 |
| canonical depth-two classes below root 1 | 17 |
| pure 120-facet targets | 1,407 |
| ridge-valid literal targets | 8 |
| ridge-valid target orbits | 8 |

The observed guarded replay completed normally in 131.258 seconds at a peak
sampled RSS of 0.414 GiB.  These resource figures are operational metadata,
not mathematical evidence.

Every one of the eight target families has literal blocker histogram
`{5:120}`, ridge histogram `{2:540}`, the scheduled face vector

```text
(14,87,319,770,1281,1477,1142,540,120),
```

and connected mod-two homology-sphere links.  The production support order
and a hostile fixed support order independently give terminal compact-label
state count zero.  The hostile replay regenerates all 108 order-four and 420
order-five classical Lanner placements.

## Lossless partition calibration

The new engine can isolate both canonical root classes and canonical
depth-two classes.  Skipped class codes are recorded before returning, so a
later raw representative cannot acquire a new class id.  Restarting deeper
quotient tables makes internal counts nonadditive, but the union of canonical
target-family keys is exhaustive.

This was replayed on completed source orbit 0.  Its four root fragments have
the same one-element canonical target union as the ordinary 231,773-node
quotient.  The unique target-bearing root has 22 depth-two classes; the union
of all 22 one-class fragments again has exactly the same one-element target
union.  Some fragments return different literal representatives of that one
orbit, which the complete `C2 x S9` family canonicalization correctly
identifies.  The frozen equivalence-certificate digest is

```text
b201093cde22f9063254bed49c384f6007834f0957f778e0d6269fc173cf0c38.
```

## Artifacts and scope

- terminal incidence transcript:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root1.json`;
- production-order downstream audit:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root1_downstream.json`;
- hostile-order downstream audit:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root1_independent_downstream.json`;
- depth-two generator:
  `scripts/generate_h9_14_record31_a1_depth_two_class_engine.py`;
- partition equivalence archive:
  `data/enumerations/h9_14_record31_a1_depth_two_partition_equivalence.json`;
- transcript verifier:
  `scripts/verify_h9_14_record31_a1_orbit2_matrix_root1.py`;
- partition certifier:
  `scripts/certify_h9_14_record31_a1_depth_two_partition.py`.

The next exact step is to run the remaining five whole-root matrix fragments
sequentially under the same RSS ceiling, union their canonical target keys,
and apply both compact-label joins.  This report does not exclude source 2,
does not complete the `a=1` schedule, constructs no Lorentzian Gram matrix,
and says nothing about noncompact finite-volume chambers.
