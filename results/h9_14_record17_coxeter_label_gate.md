# Record 17: an unrestricted compact-Coxeter no-go

## Result

The exact ordinary-polytopal selected H9/14 record-17 nerve is **not** the
nerve of a compact hyperbolic Coxeter polytope.

This selected-nerve theorem is now subsumed by the stronger
[all-orbit theorem](h9_14_record17_all_coxeter_label_gate.md), which excludes
all eight residual incidence-complex orbits and therefore the whole record-17
graph shadow on the compact track.  The proof below remains an independent,
smaller certificate for the displayed nerve.  Neither theorem excludes a
finite-volume noncompact chamber, whose ideal-face structure would require a
separate Vinberg finite-covolume analysis.

## Decisive common-core cluster

Three order-five minimal nonfaces share the four-facet core

\[
K=\{3,4,5,6\}:
\]

\[
K\cup\{0\},\qquad K\cup\{1\},\qquad K\cup\{2\}.
\]

In a compact Coxeter realization each must support an order-five Lanner
diagram.  The complete classical table has 420 labelled placements on each
support.  Exact agreement on the common core gives the successive join
counts

\[
420,\qquad 1332,\qquad 4620.
\]

The selected nerve has 110 actual rank-nine vertices.  Restricting their
diagrams to the seven-facet union \(\{0,1,2,3,4,5,6\}\) and applying the exact
finite Coxeter classification gives:

| restriction rank | distinct restrictions | surviving joined states |
|---:|---:|---:|
| 4 | 9 | 4620 |
| 5 | 16 | 80 |
| 6 | 4 | 0 |

Thus every compatible placement of the three forced Lanner diagrams makes
at least one actual vertex diagram non-elliptic.

## Why the obstruction is unrestricted

No dotted pair lies in the decisive seven-facet union.  First put label 2 on
every meeting pair in that union not already fixed by a Lanner support.  Write
the normalized Coxeter Gram matrix as \(G=I-A\), with
\(A_{ij}=\cos(\pi/m_{ij})\geq0\).  Replacing any remaining label 2 by an
arbitrary finite label \(m\geq3\) increases \(A\) entrywise.  Perron--Frobenius
monotonicity cannot lower its spectral radius, so it cannot turn any certified
non-elliptic principal diagram into an elliptic one.

Consequently no enumeration of dotted Gram values or global inertia solve is
needed: compactness already fails at a finite actual-vertex obstruction.

## Independent replay

The hostile verifier does not import the production record-17 gate, the
reusable common-core gate, or the project's Lanner/finite-Coxeter tables.  It
regenerates the 420 placements from five classical representatives, checks
their Lanner minimality, rejoins them, reconstructs all vertices directly
from the thirteen minimal nonfaces, and reruns the finite classification.  It
reproduces the same census.  Among the 80 rank-five survivors, 72 fail one of
the four rank-six restrictions and 8 fail two.

## Frozen artifacts

- `src/coxeter_polytopes/combinatorics/h9_14_record17_coxeter_labels.py`
- `scripts/certify_h9_14_record17_coxeter_labels.py`
- `scripts/independent_h9_14_record17_coxeter_labels_audit.py`
- `data/enumerations/h9_14_record17_coxeter_labels.json`
- `tests/test_h9_14_record17_coxeter_labels.py`

Certificate digest:

```text
5ae3eb4bb3dabeb0438fc56e8c92a228d88220a9e946040fc5f55a650e933437
```

Replay with:

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_record17_coxeter_labels.py --verify
PYTHONPATH=.:src python3 scripts/independent_h9_14_record17_coxeter_labels_audit.py
PYTHONPATH=.:src pytest -q tests/test_h9_14_record17_coxeter_labels.py
```
