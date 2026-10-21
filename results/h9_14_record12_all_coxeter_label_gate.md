# Record 12: all 371 residual nerves are compact-Coxeter impossible

## Theorem

Under the frozen compact-H9/14 reductions, second-layer record 12 is not the
missing-edge graph shadow of a compact hyperbolic Coxeter polytope.

The corrected whole-shadow incidence certificate leaves exactly 371 abstract
residual orbits.  Complete order-four/order-five Lanner placement joins and
every actual rank-nine vertex restriction leave zero label states on every
orbit.  This is a **whole-shadow compact-Coxeter no-go**, not a result about
one selected ordinary polytope.

Consequently the compact-facing H9/14 frontier becomes 125 of 142 graph
shadows excluded and seventeen open.  The incidence partition itself remains
`123 impossible + 19 with an ordinary-polytopal selection`.

## Exact input

The input is the corrected residual archive
`data/enumerations/h9_14_record12_residual.json`:

```text
file SHA-256        faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05
certificate SHA-256 50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4
residual orbits      371
schedule a=0,...,4  48, 95, 144, 63, 21
```

For every orbit the gate reconstructs the selected facet complements from the
two archived `5 x 5` matrices.  It independently recovers the complete
minimal-nonface clutter and checks that its enumeration of all nonface-free
nine-subsets equals the direct facet list.  Thus the actual rank-nine vertex
counts are checked twice.  They are `130,128,126,124,122` in schedules
`a=0,1,2,3,4`, respectively.  The rank-two minimal nonfaces are checked to be
exactly

```text
{0,12}, {1,12}, {0,13}, {1,13}.
```

## Exhaustive finite-label input

[Burcroff, EJC 120 (2024), Lemma 10.10](https://doi.org/10.1016/j.ejc.2024.103957)
proves globally that a hypothetical compact Coxeter 9-polytope with fourteen
facets has no multi-multiple edge.  With Definition 4.1, every meeting pair
therefore has

```text
m in {2,3,4,5}.
```

This is not a sparse-label assumption.  It is the exhaustive finite alphabet
for this compact H9/14 target.  Dotted pairs remain continuous, but cannot lie
in a finite vertex diagram.

Every order-four/five minimal nonface must support a classical Lanner diagram.
The gate uses all nine order-four and five order-five isomorphism types, giving
all `108` and `420` labelled placements.  Every actual rank-nine vertex must
be one of the finite `A/B/D/E/F/H/I_2(m)` diagrams in the globally reduced
alphabet.

## Deterministic all-orbit transcript

No support sequence is hand-picked per orbit.  The production rule starts
with archived support zero, then maximizes overlap with the cumulative cluster,
minimizes new vertices, and breaks ties by archived index.  The rule only
chooses the order in which necessary conditions are imposed; every chosen
support receives its complete Lanner placement table.

The exact stopping histogram is

| selected supports before zero | orbit count |
|---:|---:|
| 2 | 49 |
| 3 | 33 |
| 4 | 152 |
| 5 | 129 |
| 6 | 8 |

There are 1,498 cumulative stages.  The largest surviving state set at any
stage has size 864, and all 371 terminal state sets have size zero.  The
archive records at every stage the support, cluster, actual-vertex restriction
histogram, raw and compatible extension counts, survivor count, and SHA-256
of the full survivor-state set.

The combined production orbit transcript is

```text
29e7904bd2af75c692f24b0a4c502426c5316a8b2d1926ca4b07d3d7711f810a
```

## Why incremental sparse pruning loses no labeling

On an actual vertex restriction normalize the Coxeter Gram matrix as

\[
G=I-A,
\]

where `A` is symmetric and entrywise nonnegative.  The vertex is elliptic
exactly when `rho(A)<1`.  Giving any not-yet-forced meeting pair a non-right
finite label replaces `A` by an entrywise larger matrix `A'`.  Perron--Frobenius
monotonicity gives `rho(A') >= rho(A)`.  Hence a sparse restriction already
rejected as non-elliptic cannot be repaired by any omitted label.

Since every orbit reaches zero at this necessary local gate, continuous dotted
entries and the global rank-ten Lorentzian Gram equations are never reached.

## Two hostile replays

The primary hostile script imports no project module.  It independently
reconstructs all 371 nerves from matrix bits, regenerates the `9+5` Lanner
types and `108/420` placements, implements a separate finite Coxeter
classifier, and reproduces every production support choice, stage count, and
state hash.  Its full replay digest is

```text
f1b14eaab1be6cc4ac7d1eb49cd79ed1401d6f14cd9ab8ee81efd18b8b863eb3
```

A second standard-library-only replay deliberately follows a different
touching-first support order and uses no label-state symmetry quotient.  It
also reconstructs the incidence data twice and reaches zero on all 371
orbits: 187 after two supports, 164 after three, fifteen after four, and five
after five.  Its 951-stage transcript has maximum transient state count 3,481
and digest

```text
06a0d502d2f739190762f5b7e2e6442a2cadb35a4c4a96577ee3885d9d822258
```

This alternate replay is intentionally slow and is provided as optional
corroboration, not as a runtime benchmark.

## Scope

The 371 residuals are abstract incidence and mod-two homology candidates; the
theorem neither assumes nor proves that each is a PL sphere or a convex
polytope.  Such a failure could only remove candidates, so it cannot weaken
the compact no-go.

The result is compact-only.  It does not exclude a finite-volume noncompact
chamber with different ideal-face and parabolic data.  Such a claim would
require a completed Lorentzian diagram and Vinberg's full parabolic
zero-neighborhood criterion.

## Artifacts and reproduction

```text
src/coxeter_polytopes/combinatorics/h9_14_record12_all_coxeter_labels.py
scripts/certify_h9_14_record12_all_coxeter_labels.py
scripts/independent_h9_14_record12_all_coxeter_audit.py
scripts/independent_h9_14_record12_coxeter_probe.py
data/enumerations/h9_14_record12_all_coxeter_labels.json
tests/test_h9_14_record12_all_coxeter_labels.py
```

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_record12_all_coxeter_labels.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_record12_all_coxeter_audit.py
PYTHONPATH=src:. python3 scripts/independent_h9_14_record12_coxeter_probe.py --all
PYTHONPATH=src:. pytest -q tests/test_h9_14_record12_all_coxeter_labels.py
```

Frozen digests:

```text
certificate       8e06523196d30b344a0b7cec5771d968d75ab253e5c4a9ac493bdf0dec21d3d1
archive           fe6c0ee7f494a67ecac88b51aed301a4fdfe26fbddbf68034916f7b7aecc97eb
production stages 29e7904bd2af75c692f24b0a4c502426c5316a8b2d1926ca4b07d3d7711f810a
primary hostile   f1b14eaab1be6cc4ac7d1eb49cd79ed1401d6f14cd9ab8ee81efd18b8b863eb3
secondary hostile 06a0d502d2f739190762f5b7e2e6442a2cadb35a4c4a96577ee3885d9d822258
```
