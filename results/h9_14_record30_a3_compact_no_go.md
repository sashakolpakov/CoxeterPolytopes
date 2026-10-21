# H9/14 Record 30: the compact `a=3` schedule is impossible

## Certified conclusion

For graph-shadow Record 30 (missing-edge graph `3K2 + 8I`), no compact
simple nine-polytope can realize the Dehn--Sommerville schedule `a=3` while
satisfying the already-certified local vertex-link constraints.

This closes the entire `a=3` schedule, not merely one selected abstract
completion.  It is a compact-only incidence theorem: it neither excludes a
finite-volume noncompact Coxeter chamber nor constructs a new hyperbolic
Coxeter polytope.

## Scheduled data

The `a=3` schedule requires

```text
f = (14, 88, 328, 803, 1346, 1554, 1200, 567, 126),
minimal A4 nonfaces = 3,
minimal A5 nonfaces = 26.
```

The frozen positive-A4 quotient has 32 full automorphism orbits.  For every
representative, the audit reconstructs:

1. the A5 shadow forced by its three A4s;
2. the exact vertex-link Macaulay incidence intervals;
3. every admissible additional A5 family eliminating all transversals of
   order at most four; and
4. the number of order-five transversals of the completed clutter.

Once orders at most four have been eliminated, every order-five transversal
is automatically minimal.  Therefore a scheduled completion must have
exactly 126 order-five transversals.

## Exact census

Across the 32 representatives, exactly 1,896,496 addition families survive
the local interval and small-transversal gates.  Almost all occur in five
orbits:

| orbit | additions | surviving families | order-five counts |
|---:|---:|---:|---|
| 19 | 4 | 735,086 | 80--130, but never 126 |
| 20 | 4 | 1,158,881 | 80--130, but never 126 |
| 21 | 3 | 840 | 92, 97, 102, 107 |
| 22 | 3 | 840 | 92--105 or 107 |
| 23 | 3 | 840 | 96, 98--105, or 107 |

The remaining surviving families contribute nine cases, all with 130
order-five transversals.  Thus the required count 126 occurs zero times.
Ridge, topology, and Coxeter-label searches are consequently vacuous.

## Independent replay

Two algorithmically different C++20 searches agree orbit by orbit on the
forced shadows, candidate counts, small-transversal counts, surviving-family
counts, and complete order-five histograms:

- the production route uses a two-pair meet-in-the-middle enumeration for the
  four-addition cases;
- the independent route uses a pivoted exact-cover recursion with canonical
  memoization.

The frozen projection digest is

```text
7afaebc75008dd16c6bf499df24de767c5f0a9cd9bc672b563e175c7de38c5a5
```

and the certificate digest is

```text
e3afaf33526c003307d797bd0789bd2b4275d3d74c83f0a8d9792bb1ad0b9d54
```

## Artifacts

- `data/enumerations/h9_14_record30_a3_compact_no_go.json`
- `scripts/explore_h9_14_record30_a3_additions.cpp`
- `scripts/independent_h9_14_record30_a3_audit.cpp`
- `scripts/certify_h9_14_record30_a3_no_go.py`
- `tests/test_h9_14_record30_a3_no_go.py`

Together with the subsequent independently certified `a=4` census, Record 30
now has schedules `a=3,4,5,6,7` excluded.  The compact schedules `a=0,1,2`
remain open globally; selected `a=0` abstract candidates have already failed
the unrestricted compact Coxeter-label gate, but that is not an exhaustion of
the `a=0` schedule.
