# Record 31: selected-nerve compact-Coxeter no-go

## Theorem and scope

The explicit ordinary-polytopal record-31 nerve is **not** the nerve of a
compact hyperbolic Coxeter 9-polytope.

This is an unrestricted label theorem for that one selected nerve: all
allowed finite labels `2,3,4,5` and all order-five Lanner placements are
covered. It does not exclude another completion of graph-shadow record 31.
It also does not exclude a finite-volume noncompact realization: no
parabolic cusp or Vinberg finite-covolume endpoint audit is performed.

## Disjoint-Lanner reduction

Use the two disjoint missing five-faces

```text
L = {0,5,6,9,10},    R = {2,3,4,7,8}.
```

In a compact Coxeter realization both induced diagrams must be Lanner. Each
support has 420 labelled placements. The four independent clone swaps

```text
(3 4), (5 6), (7 8), (9 10)
```

preserve the complete selected nerve and its 122 actual rank-nine vertices.
Their actions on `L` and `R` factor, reducing each 420-state table to 110
representatives with orbit histogram

```text
2:10, 4:100.
```

Every pair between `L` and `R` is a meeting pair. With all cross pairs
right-angled, the two Lanner blocks are orthogonal and already give two
negative directions. Thus some finite non-right cross edge is necessary.

## Complete cross-edge census

For each of the `110 * 110` support states, all 25 cross pairs were first
tested at the smallest non-right label `m=3` against the complete actual
vertex list. The exact census is

```text
302,500 trials -> 1,932 locally surviving one-edge states.
```

Increasing the surviving edge to `m=4` or `m=5` leaves zero states. Any
configuration with at least two non-right cross edges either contains an
edge already rejected at `m=3`, or belongs to one of 36 residual pairs of
individually surviving edges. Direct actual-vertex tests reject all 36
residual pairs. Perron--Frobenius monotonicity for the normalized Gram
matrix promotes these minimal-label failures to all larger labels and all
larger cross-edge sets.

The 1,932 one-edge survivors split into four rooted weighted-block classes:

| rooted class | representative trials | disposition |
|---|---:|---|
| `A5_A5_type_0_1` | 36 | residual Lanner CSP |
| `A5_A5_type_1_1` | 900 | negative index two |
| `A5_A5_type_1_1_opposite_roots` | 36 | residual Lanner CSP |
| `A5_A5_type_4_1` | 960 | negative index two |

For the two large classes, exact leading principal minors have sign word

```text
++++-----+
```

including two sign changes after `Delta_0=1`; exact Jacobi/LDL therefore
gives inertia `(8,2,0)` for the ten-wall principal block. Cauchy
interlacing contradicts the ambient compact `H^9` Gram inertia `(9,1,4)`.

The other two classes account for 72 representative trials. For each one,
an exact minimum-remaining-values backtrack uses the complete 420-placement
tables on all seven other missing five-faces and enforces equality on every
shared Coxeter pair. It finds zero simultaneous Lanner completions. The
complete search-node histogram is

```text
1:50, 4:12, 23:2, 38:2, 63:2, 68:2, 113:2.
```

These alternatives exhaust zero, one, and at least two non-right cross
edges, so the selected compact realization is impossible.

## Independent replay and artifacts

The hostile verifier imports neither the production record-31 audit nor the
generic disjoint-Lanner gate. Starting from the frozen ordinary-polytopal
selection, it independently reconstructs the 16-element symmetry group,
both 110-state quotients, all 302,500 actual-vertex trials, the rooted block
classes, the 36 two-edge trials, the 72 Lanner-CSP searches, and both exact
inertia calculations.

- `src/coxeter_polytopes/combinatorics/h9_14_record31_coxeter_labels.py`
- `scripts/independent_h9_14_record31_coxeter_audit.py`
- `scripts/certify_h9_14_record31_coxeter_labels.py`
- `data/enumerations/h9_14_record31_coxeter_labels.json`
- `tests/test_h9_14_record31_coxeter_labels.py`

The independent summary and frozen certificate digests are

```text
4cab6a9ed8ad1ae6e921a91c41e5de00e65495b627c627f7529e6444af5b8c05
92aacb8b7d4ab3ca9faa0eefeba972d1736e6e919376bd77409af64424857b97
```

Replay with

```bash
PYTHONPATH=.:src python3 scripts/independent_h9_14_record31_coxeter_audit.py
PYTHONPATH=.:src python3 scripts/certify_h9_14_record31_coxeter_labels.py --verify
PYTHONPATH=.:src pytest -q tests/test_h9_14_record31_coxeter_labels.py
```
