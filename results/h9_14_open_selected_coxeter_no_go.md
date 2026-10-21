# Seven additional selected H9/14 nerves are compact-Coxeter impossible

## Result

The exact ordinary-polytopal selections already frozen for H9/14 graph-shadow
records

\[
0,1,3,4,5,7,9
\]

do **not** admit compact hyperbolic Coxeter labelings.  This is an unrestricted
finite-label obstruction: every complete order-four/order-five Lanner
placement compatible with the selected nerve fails ellipticity on an actual
rank-nine vertex.

This result is deliberately scoped to the seven displayed nerves.  It does
not exclude another completion over any of the seven graph shadows, and it
does not test a changed ideal-face structure for a noncompact finite-volume
polytope.

## Exact census

Unlisted intersecting pairs are set right angled at every incremental stage.
The survivor sequences are:

| record | rank-nine vertices | minimal nonfaces by order | decisive support indices | elliptic states after each join |
|---:|---:|---|---|---|
| 0 | 168 | `2:2, 4:6, 5:7` | `0,7,9,2,1` | `108 -> 84 -> 8 -> 16 -> 0` |
| 1 | 152 | `2:3, 4:4, 5:10` | `0,4,1,5,6,8` | `108 -> 84 -> 4 -> 28 -> 16 -> 0` |
| 3 | 162 | `2:2, 4:4, 5:9` | `0,7,4,5,8` | `108 -> 84 -> 24 -> 20 -> 0` |
| 4 | 146 | `2:3, 4:2, 5:10` | `0,3,10,4` | `108 -> 84 -> 8 -> 0` |
| 5 | 136 | `2:3, 4:3, 5:4` | `0,4,5` | `108 -> 84 -> 0` |
| 7 | 130 | `2:4, 5:12` | `0,7,11,1,8` | `420 -> 96 -> 144 -> 8 -> 0` |
| 9 | 126 | `2:4, 4:2, 5:8` | `0,1,5` | `108 -> 48 -> 0` |

Survivor counts need not decrease at every join because a state is a partial
label map: one old state can have several compatible extensions.  The final
zero is the relevant conclusion.

## Why the pruning is unrestricted

For each actual-vertex restriction, write the normalized Coxeter Gram matrix
as

\[
G=I-A,
\]

where (A) is symmetric and entrywise nonnegative.  The restriction is
elliptic exactly when \(\rho(A)<1\).  Changing an unassigned right angle to a
finite non-right angle produces (A'\geq A) entrywise, so
\(\rho(A')\geq\rho(A)\).  A restriction already nonelliptic in the sparse
state therefore cannot be repaired by assigning additional labels outside the
joined Lanner supports.

The order-four and order-five placement tables are complete (108 and 420
labelled placements respectively), and the compact H9/14 finite label alphabet
is exactly \(\{2,3,4,5\}\).  Thus every possible finite labeling is rejected
before global Lorentzian Gram rank and inertia become relevant.

## Independent replay

The production calculation reconstructs each nerve from its certified target
clutter and facet complements.  A separate hostile replay instead reads the
seven frozen polytopality archives, reconstructs all rank-nine faces directly
from their minimal nonfaces, regenerates the classical Lanner tables from
representatives, and uses an independently implemented finite Coxeter
classifier.  It matches every stage record and every surviving-state digest.

The hostile replay summary SHA-256 is

`29c5bc8f79024d92c33b2d2d9fb4f561600658322475e91e0d710ec9d9ea4b26`.

The canonical certificate SHA-256 is

`c0f693a5ec92abb0f146d7919d7c2a22910fa07c0576d402a9700b18ef142ed7`.

The certificate file SHA-256 is

`0dd28552e3b0cea7211b98c8d629618dc7a08b47bd884c8731c7df8cf269243b`.

Artifacts:

- `data/enumerations/h9_14_open_selected_coxeter_no_go.json`
- `src/coxeter_polytopes/combinatorics/h9_14_open_selected_coxeter_no_go.py`
- `scripts/certify_h9_14_open_selected_coxeter_no_go.py`
- `scripts/independent_h9_14_open_selected_coxeter_audit.py`
- `tests/test_h9_14_open_selected_coxeter_no_go.py`

Replay commands:

```bash
PYTHONPATH=src:. python3 scripts/certify_h9_14_open_selected_coxeter_no_go.py --verify
PYTHONPATH=src:. python3 scripts/independent_h9_14_open_selected_coxeter_audit.py
PYTHONPATH=src:. python3 -m pytest -q tests/test_h9_14_open_selected_coxeter_no_go.py
```

The focused regression suite passes `4/4`.
