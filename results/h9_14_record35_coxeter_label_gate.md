# Record 35: an unrestricted compact-Coxeter label obstruction

> **Status upgrade.** This report proves the obstruction for the originally
> selected ordinary-polytopal nerve. The later
> `h9_14_record35_whole_shadow_compact_no_go.md` exhausts all five compact
> incidence schedules and excludes the entire record-35 graph shadow.

## Theorem

The ordinary-polytopal `H9/14` record 35 is **not** the nerve of a compact
hyperbolic Coxeter polytope.

This is stronger than a failure of the sparse label ansatz.  The exact
all-right shadow already has a non-elliptic principal diagram inside an
actual rank-nine vertex for every possible choice of the five forced Lanner
blocks.  Perron--Frobenius monotonicity then proves that arbitrary additional
non-right finite Coxeter labels cannot repair it.

The result is combinatorial/local.  It does not assert nonexistence of a
finite-volume noncompact chamber with a different ideal-face structure, and
it is not a global Gram-matrix result for any other nerve.

## 1. Exact record-35 data

There are 14 facets and 102 rank-nine vertices.  The minimal nonfaces are the
five dotted pairs

```text
(0,11), (0,13), (1,12), (1,13), (2,13)
```

and the five order-five supports

```text
(2,3,4,6,8)       (0,5,7,9,10)   (1,5,7,9,10)
(4,5,7,9,10)      (3,6,8,11,12).
```

For a compact Coxeter polytope, every order-five minimal nonface must be a
Lanner diagram and every rank-nine vertex diagram must be elliptic.  The
complete classical order-five table has 420 labelled placements.

## 2. Independent two-support check

The first and last supports meet in the triple `(3,6,8)`.  Joining their two
420-row placement tables gives

```text
420 -> 11,634 compatible states.
```

Their intersections with the 102 actual vertices comprise 10 rank-four, 16
rank-five, and 3 rank-six restrictions.  Exactly 24 states pass them all.
They form one orbit of size 24 under
`S2_{(2,4)} x S3_{(3,6,8)} x S2_{(11,12)}`.  Thus the independently observed
`11,634 -> 24` count is reproduced, but this cluster is not the obstruction.

## 3. Reusable common-four-core gate

Put

```text
K = (5,7,9,10),       (a,b,c) = (0,1,4).
```

The other three supports are `K union {a}`, `K union {b}`, and
`K union {c}`.  Group the 420 Lanner extensions by their labelled restriction
to `K`.  The number `k` of extensions of a fixed core code has histogram

| `k` | core codes |
|---:|---:|
| 1 | 36 |
| 2 | 30 |
| 3 | 40 |
| 4 | 51 |

Consequently

```text
sum k   = 420,
sum k^2 = 1,332,
sum k^3 = 4,620.
```

These are independently the successive exact table-join counts.  The
complete actual-vertex restriction family is archived, rather than inferred
from the support pattern.  Its cumulative gate is

| order | actual restrictions | surviving states |
|---:|---:|---:|
| 4 | 13 | 4,620 |
| 5 | 12 | 96 |
| 6 | 4 | 0 |

The 12 rank-five restrictions have the form

```text
two of {a,b,c} + three of K,
```

and the four rank-six restrictions have the form

```text
{a,b,c} + three of K.
```

For record 35 this restriction family is invariant under the natural
`S3_{a,b,c} x S4_K` action.  The 96 intermediate survivors form exactly two
orbits.

### Orbit of size 24

One representative has non-right labels

```text
(0,10)=3, (1,10)=3, (4,10)=3,
(5,9)=3, (5,10)=3, (7,9)=5.
```

On facets `(0,1,4,5,7,10)` it contains
`affine D4 + A1`: node 10 is joined by simple edges to `0,1,4,5`, while 7 is
isolated.  In that order the simply-laced Cartan null vector is
`(1,1,1,1,0,2)`.  This six-set lies in the actual vertex
`(0,1,2,3,4,5,6,7,10)`.

### Orbit of size 72

One representative has non-right labels

```text
(0,10)=3, (1,10)=3, (4,9)=3,
(5,7)=5, (7,9)=3, (9,10)=3.
```

On facets `(0,1,4,7,9,10)` it is `affine D5`.  In that order its Cartan null
vector is `(1,1,1,1,2,2)`.  The six-set lies in the actual vertex
`(0,1,2,3,4,6,7,9,10)`.

The size-24 orbit fails two of the four rank-six restrictions; the size-72
orbit fails one.  Hence all 4,620 common-core states fail an actual vertex
test.  The two support clusters are pair-disjoint across the decomposition,
so the full five-Lanner sparse family has exactly

```text
11,634 * 4,620 = 53,749,080
```

compatible states and no locally elliptic member.

The exported function
`audit_common_four_core_lanner_shadow(outer_facets, common_core,
vertex_facets)` replays this pattern for another `H9/14` nerve.  Its
conclusion is deliberately tied to the supplied complete actual vertex list;
three abstract common-core supports alone are not asserted to be impossible.

## 4. Why no additional label can rescue it

For a principal vertex diagram use normalized Gram convention

```text
G = I - A,       A_ij = cos(pi/m_ij) >= 0.
```

The diagram is elliptic exactly when `G` is positive definite, equivalently
when the Perron eigenvalue `lambda_max(A)` is less than 1.  The sparse shadow
sets every pair outside the forced Lanner supports to `m=2`, hence to
`A_ij=0`.  Replacing any such entry by an arbitrary finite `m>=3` increases
the symmetric nonnegative matrix `A` entrywise.  Perron--Frobenius
monotonicity gives

```text
lambda_max(A_full) >= lambda_max(A_shadow) >= 1.
```

Thus the corresponding principal vertex Gram remains non-positive-definite.
Edges outside the witness principal set leave that bad principal submatrix
unchanged; edges inside it can only worsen the inequality.  This proves the
unrestricted compact-Coxeter no-go, not merely sparse UNSAT.

As a direct finite check, changing exactly one of the three external pairs
`(0,1),(0,4),(1,4)` to each of `3,4,5` gives
`4,620 * 3 * 3 = 41,580` states and zero survivors.  The monotonic argument
covers every larger shell and every finite label `m>=3` without enumeration.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/certify_h9_14_record35_coxeter_labels.py --verify
PYTHONPATH=src python3 -m pytest -q tests/test_h9_14_record35_coxeter_labels.py
```

Machine-readable certificate:
`data/enumerations/h9_14_record35_coxeter_labels.json`.
