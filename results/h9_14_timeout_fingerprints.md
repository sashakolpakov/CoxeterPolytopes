# Exact fingerprints of the 23 H9/14 timeout shadows

## Scope

This note ranks the 23 theorem-undecided graph shadows present immediately
after the record-23 certificate.  The calculation is reproduced by

```sh
PYTHONPATH=src python3 scripts/explore_h9_14_timeout_fingerprints.py
```

It is a search map, not a non-existence certificate.  Every entry is exact,
but the rank is only a heuristic ordering for the next blocker proof.

The 23 indices are

```text
0, 1, 3, 4, 5, 6, 7, 9, 12, 14, 16, 17, 18,
31, 32, 34, 35, 37, 43, 44, 46, 47, 74.
```

At that snapshot they were the 23 records remaining from the frozen global
facet archive's 27 bounded `unknown` records after removing the explicit
product selections 30 and 42 and the independently excluded records 92 and
93.

## Fingerprints

For each graph `H`, the script reconstructs

- every inclusion-minimal graph cover and its size distribution;
- every Dehn--Sommerville schedule and the exact double-counted sum of
  vertex-link `g4`;
- every independent A4, its orbit under the full automorphism group of the
  non-isolated graph together with the symmetric isolate block, and the exact
  vertex/edge/triple A4 interval gate;
- for the `a=1` branch, the complete forced five-shadow, strict vertex-link
  `g4` gate, and aggregate Macaulay A5-incidence interval; and
- the hit-pattern alphabet of all independent A5s against the minimal graph
  covers.

In the table, `A4` is

```text
all one-A4 orbits -> local-interval orbits -> shadow-strict orbits,
```

`covers` is the size histogram of inclusion-minimal graph covers, and
`patterns` is the number of nonempty A5 hit-pattern cells.  The ranking is
lexicographic in the number of strict-sum schedules, shadow-surviving A4
orbits, the `a=0` A5 count, pattern alphabet, and total A4 orbit count.

| rank | record | edges | strict-sum `a` | `a=0` A5s | covers | patterns | A4 |
|---:|---:|---:|---|---:|---|---:|---|
| 1 | 74 | 6 | 0,1 | 4 | `3:4` | 5 | `2 -> 1 -> 0` |
| 2 | 46 | 6 | 0,1 | 4 | `3:4` | 7 | `8 -> 3 -> 0` |
| 3 | 44 | 6 | 0,1 | 4 | `3:4` | 8 | `6 -> 3 -> 0` |
| 4 | 47 | 6 | 0,1 | 4 | `3:4` | 11 | `7 -> 3 -> 0` |
| 5 | 43 | 6 | 0,1 | 4 | `3:4` | 10 | `11 -> 5 -> 1` |
| 6 | 37 | 5 | 0,1,2 | 5 | `3:5` | 10 | `8 -> 4 -> 0` |
| 7 | 35 | 5 | 0,1,2 | 5 | `3:5` | 14 | `15 -> 7 -> 2` |
| 8 | 34 | 4 | 0,1,2,3 | 6 | `3:6` | 12 | `4 -> 3 -> 1` |
| 9 | 32 | 4 | 0,1,2,3 | 6 | `3:6` | 18 | `10 -> 7 -> 4` |
| 10 | 16 | 5 | 0,...,4 | 8 | `2:1,3:2` | 4 | `4 -> 3 -> 2` |
| 11 | 18 | 5 | 0,...,4 | 8 | `2:1,3:2` | 6 | `9 -> 6 -> 3` |
| 12 | 17 | 5 | 0,...,4 | 8 | `2:1,3:2` | 6 | `8 -> 6 -> 4` |
| 13 | 14 | 5 | 0,...,4 | 8 | `2:1,3:2` | 6 | `10 -> 7 -> 5` |
| 14 | 31 | 4 | 0,...,5 | 9 | `2:1,3:3` | 9 | `9 -> 7 -> 6` |
| 15 | 12 | 4 | 0,...,7 | 12 | `2:2` | 3 | `3 -> 3 -> 3` |
| 16 | 9 | 4 | 0,...,7 | 12 | `2:2,3:1` | 5 | `5 -> 4 -> 3` |
| 17 | 7 | 4 | 0,...,7 | 12 | `2:2,3:1` | 6 | `11 -> 9 -> 7` |
| 18 | 5 | 3 | 0,...,8 | 13 | `2:2,3:2` | 9 | `8 -> 7 -> 6` |
| 19 | 6 | 3 | 0,...,10 | 16 | `2:3` | 4 | `2 -> 2 -> 2` |
| 20 | 4 | 3 | 0,...,10 | 16 | `2:3` | 6 | `5 -> 5 -> 5` |
| 21 | 3 | 2 | 0,...,13 | 20 | `2:4` | 9 | `3 -> 3 -> 3` |
| 22 | 1 | 3 | 0,...,13 | 22 | `1:1,3:1` | 3 | `5 -> 4 -> 3` |
| 23 | 0 | 2 | 0,...,16 | 26 | `1:1,2:1` | 3 | `4 -> 4 -> 4` |

The aggregate A5-incidence Macaulay gate leaves exactly the same one-A4
orbit count as the shadow-strict column in every row of this table.

## Consequences for the search

Record 74 was the clear first target. Its graph is `K4 + 10I`; the exact
argument in `results/h9_14_record74_certificate.md` excludes it without an
A5 census. At that historical stage this promoted the theorem-facing
partition to

```text
118 independently impossible + 2 explicit selections + 22 undecided.
```

At that same snapshot, records 46, 44, and 47 were the next natural common
family. In all three,
strict local `g`-vectors eliminate the entire `a=1` branch, leaving only four
minimal A5s in `a=0`.  Their blocker problems have four size-three graph
covers and hit-pattern alphabets of size 7, 8, and 11 respectively.  They
were naturally attacked in that order with a common four-cover
column-pattern compiler.

That programme is now complete for the five highest-ranked rows. Record 46
is independently impossible, as is record 74; records 43, 44, and 47 have
exact ordinary-polytopal multiwedge realizations. Subsequent exact work also
excludes record 37, while record 35 has an exact rational multiwedge
realization. A subsequent exact certificate also excluded record 34.  At
that historical stage the partition was

```text
121 independently impossible + 6 ordinary-polytopal + 15 undecided.
```

Later certificates promote record 32 to ordinary polytopality and exclude
record 16 as a whole graph shadow.  The current partition is therefore

```text
122 independently impossible + 7 ordinary-polytopal + 13 undecided.
```

The table remains the exact fingerprint of the earlier 23-shadow frontier,
not a claim that all 23 are still unresolved.  Among its surviving rows it
now ranks records 18, 17, and 14 first.

Records with two-vertex minimal covers require a different decomposition:
a size-two cover has two extra vertices available before reaching the
forbidden order-four threshold, so the record-23 pair-transversal machinery
is the relevant model.  This explains why a small pattern alphabet alone
does not make records 0, 1, 6, 12, or 16 immediate.
