# Record 32, the distinct `a=0` completion

## Result

There is a second explicit completion over H9/14 missing-edge record 32.  Its
minimal nonfaces are

\[
\begin{aligned}
 &(0,11),(1,12),(0,13),(2,13),\\
 &(4,5,7,8,9),\ (1,2,3,6,10),\ (2,3,4,6,10),\\
 &(3,4,6,10,11),\ (5,7,8,9,12),\ (5,7,8,9,13).
\end{aligned}
\]

Its blocker consists of exactly 114 five-sets.  Taking their complements as
the rank-nine facets gives

\[
f=(14,87,318,762,1254,1428,1092,513,114).
\]

The exact global audit finds 513 ridges, each in two facets; a connected dual
graph; the required Dehn--Sommerville schedule; spherical Euler
characteristic for every link; and the mod-2 homology of the expected sphere
for every link.  This freezes an explicit abstract homology-sphere candidate.
It does **not** assert that the complex is the boundary of an ordinary convex
polytope.

The same selected nerve is nevertheless impossible as a compact hyperbolic
Coxeter nerve.

## The common-four-core obstruction

Use just the three missing five-faces

\[
 K\cup\{4\},\qquad K\cup\{12\},\qquad K\cup\{13\},
 \qquad K=\{5,7,8,9\}.
\]

In a compact Coxeter polytope, each is an order-five Lanner diagram.  The
complete classical table contains 420 labelled placements on a fixed
five-set.  Joining the three tables by equality on their shared core gives

\[
             420\longrightarrow 1332\longrightarrow 4620
\]

compatible labelled triples.  This is an exact table join, not a sampling or
a choice of a preferred golden decoration.  The extension multiplicities of
the possible labelled core diagrams are

| extensions of a core | number of labelled cores |
|---:|---:|
| 1 | 36 |
| 2 | 30 |
| 3 | 40 |
| 4 | 51 |

Intersecting all 114 actual rank-nine vertices with
\(K\cup\{4,12,13\}\) produces 13 distinct restrictions of rank 4, 12 of
rank 5, and 4 of rank 6.  Exact finite-Coxeter classification gives the
cumulative census

| restrictions imposed | surviving triples |
|---|---:|
| all rank 4 | 4620 |
| rank 4 and 5 | 96 |
| rank 4, 5, and 6 | 0 |

Thus every compatible Lanner triple contains a non-elliptic principal
subdiagram inside an actual polytope vertex.

This sparse failure is unrestricted in all unenumerated labels.  Write a
normalized Gram matrix as \(G=I-A\), where \(A\) is the nonnegative cosine
adjacency matrix.  Changing a right angle to any non-right finite Coxeter
label increases \(A\) entrywise.  Perron--Frobenius monotonicity cannot lower
its spectral radius below one, hence cannot turn an affine or indefinite
principal diagram into an elliptic one.  Labels incident to facets outside a
failed restriction are irrelevant because every principal subdiagram of an
elliptic vertex diagram must itself be elliptic.

Consequently:

> The displayed `a=0` selected completion of record 32 is not the nerve of a
> compact hyperbolic Coxeter polytope, for any assignment of its remaining
> finite labels or dotted Gram entries.

The scope is deliberately narrow.  This rules out this selected nerve, not
every completion sharing graph6 shadow `M?????????O?O?g??`; it proves no
uniqueness statement for the record-32 completion problem.

## Reproducibility

The frozen files are

- `src/coxeter_polytopes/combinatorics/h9_14_record32_a0.py`,
- `scripts/certify_h9_14_record32_a0.py`,
- `scripts/independent_h9_14_record32_a0_audit.py`,
- `data/enumerations/h9_14_record32_a0_selection.json`, and
- `tests/test_h9_14_record32_a0.py`.

Run

```console
PYTHONPATH=.:src python3 scripts/certify_h9_14_record32_a0.py
PYTHONPATH=.:src python3 scripts/independent_h9_14_record32_a0_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record32_a0.py
```

The frozen critical digest is

```text
14c3f9f5a425eb6934d7f53a7a7f071467eba7cf2e4a64b7cf9a194d3490d0bb
```
