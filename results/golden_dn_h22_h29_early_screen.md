# Golden `D_n` / `D24+` early Vinberg screen

The norm-two search for

\[
\langle-2\varphi\rangle\oplus D_n,qquad 22\le n\le29,
\]

was run exactly at controls `x=1,phi,phi^2`, with exact incremental cone rays
after all four accepted positive walls.  The D24+ run includes all three
nonzero spinor residues in `O/2O`: `1`, `phi`, and `1+phi`.

The accepted pattern stabilizes by zero padding in every tested dimension:

- at `x=1`: `phi(e1+e2)`;
- at `x=phi`: `phi(2e1)` and `phi(e1+e2+e3+e4)`;
- at `x=phi^2`: `(1+phi)(e1+e2)+phi(e3+e4+e5+e6)`.

At the first cap, the D24 affine alcove has exactly 25 rays, classified as
`(finite,cusp,hyperideal)=(3,0,22)`.  Through `x=phi^2`, the complete current
cone census is:

| dimension | walls | total rays | finite | cusps | hyperideal |
|--:|--:|--:|--:|--:|--:|
| 22 | 26 | 254 | 28 | 0 | 226 |
| 23 | 27 | 268 | 28 | 0 | 240 |
| 24 | 28 | 282 | 28 | 0 | 254 |
| 25 | 29 | 296 | 28 | 0 | 268 |
| 26 | 30 | 310 | 28 | 0 | 282 |
| 27 | 31 | 324 | 28 | 0 | 296 |
| 28 | 32 | 338 | 28 | 0 | 310 |
| 29 | 33 | 352 | 28 | 0 | 324 |

All cones are irredundant, but none has finite volume.  D24+ has the same
H-through-`phi^2` count as D24: no accepted early wall lies in a spinor glue
coset.  Hyperideal rays rise sharply rather than falling, so the requested
trajectory gate does not justify extending this branch to H30.

This archive is a norm-two early screen.  For the nonunimodular root lattices
`D_n`, a complete local mirror-norm audit remains necessary before the
prefix can be called exhaustive over all reflections.  Later controls or
additional norm classes could cut current rays; no nonreflectivity or novelty
claim is made.

- Backend: `src/coxeter_polytopes/search/golden_dn.py`
- Generator: `scripts/search_golden_dn_early.py`
- Tests: `tests/test_golden_dn.py`
- Archive: `data/enumerations/golden_dn_h22_h29_through_phi2.json`
