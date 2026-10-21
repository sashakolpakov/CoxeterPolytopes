# Record 31, `a=1`, source orbit 2: matrix root class 2 closed

## Outcome

Canonical root class 2 is now completely excluded from the compact Coxeter
problem.  Its 26 depth-two classes were exhausted in thirteen disjoint
two-class fragments.  Seven literal fragment targets collapse to exactly
three `C2 x S9` target orbits.  All three pass the connected mod-two
homology-sphere audit, but production-order and hostile-order unrestricted
compact Coxeter-label joins both leave zero states.

This closes root class 2 only.  Together with the earlier root-1 tranche, two
of the six canonical root classes of source orbit 2 are closed.  Root classes
3 through 6 remain open, so neither source 2 nor the full Record-31 `a=1`
schedule is yet excluded.  Nothing here concerns noncompact finite-volume
chambers.

## Exact partition

The unpartitioned root-2 attempt was stopped by the 0.700 GiB watchdog at 15
million nodes and produced no terminal JSON.  That disposition is not used as
evidence.  The lossless depth-two partition was then run under a stricter
0.60 GiB ceiling:

| depth-two classes | nodes | pure 120-facet targets | literal ridge targets |
|---|---:|---:|---:|
| 1–2 | 4,707,403 | 2 | 0 |
| 3–4 | 4,503,548 | 0 | 0 |
| 5–6 | 1,085,753 | 0 | 0 |
| 7–8 | 5,013,325 | 0 | 0 |
| 9–10 | 6,735,369 | 0 | 0 |
| 11–12 | 5,133,384 | 0 | 0 |
| 13–14 | 2,245,049 | 18 | 1 |
| 15–16 | 3,440,668 | 1 | 0 |
| 17–18 | 7,662,165 | 2 | 0 |
| 19–20 | 848,958 | 70 | 2 |
| 21–22 | 5,313,578 | 477 | 3 |
| 23–24 | 5,335,661 | 24 | 1 |
| 25–26 | 3,019,403 | 0 | 0 |

Every row terminated normally.  Sampled peaks ranged from 0.030 to 0.361
GiB.  Fragment node counts, pure-target counts, and literal multiplicities
are nonadditive under restarted deeper quotient tables and are not combined
into a theorem-facing total.  The depth-two class ranges themselves are
disjoint and cover exactly `1,...,26`.

The seven literal ridge targets canonicalize to three families.  Direct
blocker/ridge replay gives `{5:120}` and `{2:540}` for each.  All three have
the scheduled face vector

```text
(14,87,319,770,1281,1477,1142,540,120)
```

and pass every connected mod-two link audit.  The two unrestricted label
joins independently leave terminal state count zero; the hostile replay
reconstructs the 108 order-four and 420 order-five classical Lanner
placements.

## Artifacts and limits

- fragments:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_d2_*.json`;
- exact canonical union:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_target_union.json`;
- production downstream:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_downstream.json`;
- hostile downstream:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root2_independent_downstream.json`;
- assembler:
  `scripts/assemble_h9_14_record31_a1_orbit2_matrix_partition.py`;
- verifier:
  `scripts/verify_h9_14_record31_a1_orbit2_matrix_root2.py`.

The canonical-union archive has SHA-256
`4dcd0af6e502275535d7790dc4eac8e535ec2d2b833891e35351e56cf17dad93`.
This is a terminal exact root-class theorem, not a timeout promotion, a Gram
realization, a full source-orbit theorem, or a noncompact statement.
