# Record 6: exact incidence no-go

## Result

The H9/14 graph-shadow record

```text
M???????????_?_?_
```

has missing-edge graph \(K_3+11I\).  It cannot be the nerve of a compact
hyperbolic Coxeter 9-polytope.  In fact, the obstruction occurs earlier:
there is no simplicial 9-polytope boundary satisfying the frozen
Coxeter-nerve incidence necessities with this missing-edge graph.

This is a nonexistence theorem, not a new hyperbolic example.

## Three two-graphs

Call the three vertices of the missing triangle *special*.  A facet contains
at most one special vertex.  Its five-element complement therefore has one
of two forms:

- two special vertices and three of the eleven isolates, forming one of
  three families \(Y_s\); or
- all three special vertices and two isolates, forming an edge family \(Z\).

The ridge equations say that every \(Y_s\) is a two-graph.  After choosing one
isolate as switching base, write \(Y_s=\delta X_s\) for an ordinary graph on
ten vertices.  The stronger ridge degree \(0/2\) condition is exactly the
exclusion of triangles and induced \(2K_2\)'s in \(X_s\).  The fourth family
satisfies

\[
\delta Z=Y_0\mathbin\triangle Y_1\mathbin\triangle Y_2.
\]

The complete `geng -q -t 10` census contains 12,172 unlabeled triangle-free
graphs.  Exactly 308 are also induced-\(2K_2\)-free.  Exact pair-coverage,
minimal-nonface-order, ridge-connectivity, \(g\)- and Macaulay tests leave 42
fixed-switching-base link types, with invariant histogram

| link facets | missing triangles \(q\) | \(g_4=r\) | fixed-base types |
|---:|---:|---:|---:|
| 45 | 2 | 1 | 2 |
| 48 | 1 | 1 | 2 |
| 49 | 1 | 2 | 10 |
| 50 | 1 | 3 | 4 |
| 52 | 0 | 2 | 6 |
| 53 | 0 | 3 | 12 |
| 54 | 0 | 4 | 5 |
| 55 | 0 | 5 | 1 |

The only source-link orbits needed by the coupled search are particularly
rigid.  Type \(E=(48,1,1)\) is the transversal hypergraph of a
\(3+4+4\) partition, giving exactly 5,775 labelled copies.  Type
\(H=(45,2,1)\) is the transversal hypergraph of a \(3+3+5\) partition,
giving exactly 4,620 labelled copies.  Direct membership checks prove that
both fixed-base E types and both fixed-base H types lie in these complete
labelled orbits.  The other types occur only as derived targets and are
tested by the full link predicate.

## Why there are exactly ten branches

Let \((y_s,q_s,r_s)\) be the signature of special link \(s\), and let
\(a_0\) count isolate-only missing four-faces.  Then the global number of
missing four-faces is \(a=\sum q_s+a_0\).  Dehn--Sommerville gives

\[
f_8=150-2a,
\]

whereas the link decomposition gives \(f_8=\sum y_s+|Z|\).  Consequently

\[
|Z|=150-2(\sum q_s+a_0)-\sum y_s
     =\sum(q_s-r_s)-2a_0.
\]

Enumerating every unordered triple of the eight admissible link signatures
and every nonnegative \(a_0\), then retaining \(|Z|\ge0\), derives exactly
ten cases:

```text
EEE,z0,a0=0
AHH,z0,a0=0     EEH,z1,a0=0     EFH,z0,a0=0
EHH,z2,a0=0     FHH,z1,a0=0     GHH,z0,a0=0
HHH,z3,a0=0
EHH,z0,a0=1     HHH,z1,a0=1
```

Here \(F=(49,1,2)\), \(G=(50,1,3)\), and \(A=(52,0,2)\).

## Complete relative-orbit search

A common permutation of the eleven isolates sends one E or H source link to
the fixed canonical partition.  The program then enumerates every labelled
copy of the second source link and every labelled \(Z\)-edge support.  Thus
the enumeration covers all relative alignments; it is not a sample of
canonical graph labels.

The residual counts are:

| branch | facet-count candidates | link-compatible | all independent triples covered | forced \(a_0\) matched |
|---|---:|---:|---:|---:|
| EEE,z0 | 54 | 54 | 36 | 0 |
| AHH,z0 | 0 | 0 | 0 | 0 |
| EEH,z1 | 2,790 | 0 | 0 | 0 |
| EFH,z0 | 0 | 0 | 0 | 0 |
| EHH,z2 | 25,086 | 0 | 0 | 0 |
| FHH,z1 | 1,953 | 1,440 | 0 | 0 |
| GHH,z0 | 0 | 0 | 0 | 0 |
| HHH,z3 | 780 | 780 | 0 | 0 |
| EHH,z0,\(a_0=1\) | 60 | 60 | 0 | 0 |
| HHH,z1,\(a_0=1\) | 1,130 | 0 | 0 | 0 |

Every one of the 1,440 FHH residuals, 780 HHH residuals, and 60 EHH
residuals misses exactly one independent isolate triple.  In the EEE branch,
18 of 54 residuals likewise miss one isolate triple.  The remaining 36 cover
all independent triples but each has exactly one isolate-only missing
four-face, contradicting the branch's forced value \(a_0=0\).  Hence all ten
branches are impossible.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record6.py`
- `scripts/certify_h9_14_record6.py`
- `scripts/independent_h9_14_record6_audit.py`
- `data/enumerations/h9_14_record6_incidence_no_go.json`
- `tests/test_h9_14_record6.py`

The internal certificate digest is

```text
cc0c3186b16cafc3c3a9572ea15f6d24880ad7bb3b18da40cf8c635d064f78ff
```

The import-independent hostile replay has summary digest

```text
42685136e71e5964620fdb7cd484784345f9e8bc63c4934d08dcd52e2267635e
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record6.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record6_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record6.py
```
