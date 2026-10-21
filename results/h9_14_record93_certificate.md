# Solver-independent exclusion of H9/14 record 93

> **Current status.** Later record-54/23/74/46/37/34 certificates, together with
> the multiwedge realizations of records 43, 44, 47, and 35, give the current
> partition \(121+6+15\); see `results/h9_14_facet_completion.md`.

## Result

**Record 93 is independently impossible.**  Its graph6 string is

```text
M?????????W?o?o@?
```

and its missing edges are

\[
\{0,11\},\{1,11\},\{0,12\},\{1,12\},
\{0,13\},\{1,13\},\{11,13\}.
\]

The proof uses only integer arithmetic: the vertex-link \(g_4\) double count
from the record-92 argument, Dehn--Sommerville, and an exhaustive blocker
classification. It does not trust the bounded Z3 run.

**Historical stage snapshot.** This second promotion changed the
theorem-facing partition of the 142 graph shadows to

\[
115\text{ independently impossible}
\; +\;2\text{ explicit ordinary-polytopal selections}
\; +\;25\text{ undecided}.
\]

The subsequent record-54 certificate then advanced the partition to
\(116+2+24\). These are proof-history snapshots; the current status is stated
above.

## 1. Only the \(a=0\) schedule can survive

Let \(a\) be the number of missing four-faces.  As proved in the record-92
certificate, every one of the fourteen vertex links has strict
\(g_4\geq1\), while incidence double-counting gives

\[
\sum_v g_4(\operatorname{lk}v)
=5f_4-24f_3+63f_2-112f_1+14\cdot126
=18-6a.
\]

The archived schedules have \(a=0,1,2\), with sums 18, 12, and 6.  Since the
sum must be at least 14, only \(a=0\) remains.

There are 1,008 independent five-sets and the \(a=0\) schedule has
\(f_4=1005\), so exactly three independent five-sets
\(B_1,B_2,B_3\) are minimal nonfaces.  The Lanner order bound excludes larger
minimal nonfaces, and the same schedule requires \(f_8=78\).  Thus a pure
target would have exactly 78 inclusion-minimal transversals, all of order
five.

## 2. Blocker hit patterns

Put \(I=\{2,\ldots,10\}\), the nine vertices isolated in the missing-edge
graph.  Its inclusion-minimal vertex covers are

\[
C_1=\{0,1,11\},\qquad
C_2=\{0,1,13\},\qquad
C_3=\{11,12,13\}.
\]

For an independent five-set \(B\), put
\(p(B)=\{i:B\cap C_i\ne\varnothing\}\).  Exactly five patterns occur.

| \(p(B)\) | structure | count |
|---|---|---:|
| \(\varnothing\) | \(B\in\binom I5\) | 126 |
| \(\{1,2\}\) | \(B\subset\{0,1\}\cup I\), meeting \(\{0,1\}\) | 336 |
| \(\{3\}\) | \(B=\{12\}\cup Q,\ Q\in\binom I4\) | 126 |
| \(\{1,3\}\) | \(B=\{11\}\cup Q,\ Q\in\binom{I\cup\{12\}}4\) | 210 |
| \(\{2,3\}\) | \(B=\{13\}\cup Q,\ Q\in\binom{I\cup\{12\}}4\) | 210 |

If a cover \(C_i\) misses no \(B_j\), it is a 3-transversal.  If the
\(B_j\)'s missed by \(C_i\) have a common vertex, adding that vertex gives a
transversal of order at most four.  Absence of this immediate purity
obstruction forces the three hit patterns to be pairwise disjoint.

The only compatible pattern multisets are

\[
(0,0,0),\qquad(0,0,p),\qquad
(0,\{1,2\},\{3\}).
\]

Every \((0,0,p)\) case has a small transversal because its two zero-pattern
members are five-subsets of the nine-set \(I\) and therefore intersect.
There remain two structural families.

## 3. Three zero-pattern faces

Suppose all \(B_j\) are five-subsets of \(I\).  Put
\(R_j=I\setminus B_j\), so \(|R_j|=4\).  Avoiding a small transversal is
equivalent to

\[
B_1\cap B_2\cap B_3=\varnothing,
\qquad\text{or equivalently}\qquad
R_1\cup R_2\cup R_3=I.
\]

Write

\[
x_{ij}=|R_i\cap R_j|,
\qquad t=|R_1\cap R_2\cap R_3|.
\]

The union condition gives

\[
x_{12}+x_{13}+x_{23}-t=3.
\]

The number of isolate pairs contained in at least one \(R_j\) is

\[
18-\sum_{i<j}\binom{x_{ij}}2+\binom t2.
\]

Consequently the number of five-transversals is

\[
N=3\left(\binom92-18+
\sum_{i<j}\binom{x_{ij}}2-\binom t2\right).
\]

The intersection equation permits only correction terms 0, 1, or 3.  Thus

\[
N\in\{54,57,63\}.
\]

The independent enumeration gives four exact complement-intersection
signatures.  In the table, the first entry is the sorted triple
\((x_{12},x_{13},x_{23})\).

| pairwise signature | \(t\) | \(N\) | unordered triples |
|---|---:|---:|---:|
| \((1,1,1)\) | 0 | 54 | 7,560 |
| \((0,1,2)\) | 0 | 57 | 15,120 |
| \((1,1,2)\) | 1 | 57 | 7,560 |
| \((0,0,3)\) | 0 | 63 | 1,260 |

Thus the two distinct \(N=57\) signatures together have multiplicity 22,680.

None has the required 78 facets.

## 4. Pattern \((0,\{1,2\},\{3\})\)

Choose \(B_0\in\binom I5\) and put \(R=I\setminus B_0\), now a 4-set.
The empty-intersection conditions force

\[
B_3=R\cup\{12\},
\]

and \(B_{12}\) is any five of the six elements
\(R\cup\{0,1\}\).  Hence there are

\[
\binom95\cdot6=756
\]

families.

The missing-edge graph has 3, 5, and 1 active covers of orders 3, 4, and 5,
so it has

\[
3\binom92+5\binom91+1=154
\]

five-element vertex covers.  As in record 92, the disjoint-from counts for
\(B_0,B_{12},B_3\) are 39, 15, and 25, and those forbidden classes are
pairwise disjoint.  Every such family therefore has

\[
154-39-15-25=75
\]

five-transversals, again not 78.

This exhausts all choices.  Any choice omitted by the classification already
has a transversal of order at most four and fails purity; every remaining
choice has the wrong number of five-transversals.  Record 93 is impossible.

## 5. Independent verification and exact scope

The pure verifier reconstructs all source data and enumerates:

- 1,008 independent five-sets and their exact pattern table;
- 12,605,586 pairwise-pattern-compatible triples;
- 32,256 triples with no transversal of order at most four; and
- the histogram
  \(\{54:7560,57:22680,63:1260,75:756\}\).

As a defense-in-depth check, it separately reconstructs the four
complement-intersection signatures above and checks that all family
multiplicities sum to 32,256.

The schedule-split backend had returned UNSAT only for \(a=0\), with timeouts
for \(a=1,2\).  The new proof independently excludes all three branches; the
backend archive is recorded only as untrusted corroboration.

Artifacts:

- `src/coxeter_polytopes/combinatorics/h9_14_record93.py`;
- `scripts/certify_h9_14_record93.py`;
- `data/enumerations/h9_14_record93_certificate.json`; and
- `tests/test_h9_14_record93_certificate.py`.

Reproduce with

```text
PYTHONPATH=src python3 scripts/certify_h9_14_record93.py --verify-only
PYTHONPATH=src pytest -q tests/test_h9_14_record93_certificate.py
```

This excludes one additional graph shadow.  It does not by itself prove the
full compact H9/14 non-existence theorem.
