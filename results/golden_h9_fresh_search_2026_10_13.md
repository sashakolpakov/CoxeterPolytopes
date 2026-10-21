# Fresh compact H9 searches, 13 October 2026

No new compact polytope has been certified. The useful finding is that the
repository's infinite-volume proof for the **norm-two** chamber of
`<-2phi> + E8 + <2>` does not decide the chamber of all integral reflections.
The omitted norm-four mirrors cut that chamber, beginning at very low height.
Here `phi=(1+sqrt(5))/2` and the spatial sign convention is positive.

## The omitted Pell mirror

The integral vector

`r=(phi-1, 0,0,0,0,0,0,0,0, phi)`

has norm four. Its reflection on the `(x,u)` plane is

`[[phi,-1],[phi,-phi]]`.

This is exactly the Pell involution used as a symmetry of the complete
norm-two chamber in `og10_golden_h9_infinite_stabilizer.md`. In the full
reflection arrangement it is a mirror, not an obstruction to finite volume.
Its Vinberg height is `(phi-1)/2`.

For this lattice every primitive integral reflection has a representative
of norm two or four. Indeed `2 B^-1` is integral, so reflection integrality
and primitivity give `q | 4`; evenness gives `2 | q`. The prime `(2)` is
inert in `Z[phi]`, and totally positive units are squares. At norm four,
reflection integrality is equivalent to even E8 pairings (and hence even
E8 coordinates), because E8 is unimodular.

The height-zero system is exactly `E8+A1`. To check this without assuming it,
write an E8 coordinate as `a+phi b`. Its trace norm is

`2 ||a+b/2||² + (5/2) ||b||²`.

If `b != 0`, this is at least five; if `b=0` and `a!=0`, it is at least four.
At norm two the trace budget forces a rational E8 root or the A1 root.
At norm four, E8 coordinates have the form `2z`, so their trace contribution
already exceeds the available budget unless `z=0`; the remaining equation
`u²=2` has no solution in `Q(sqrt(5))`.

## Exact bounded results

The control-height convention is `H=2phi*x²/q`. Controls are enumerated in
both real embeddings, sorted by height, and each new wall is required to
pair nonpositively with all earlier walls. All root equations, divisibility
conditions, comparisons, and cone certifications use exact arithmetic.

| host | root norms | completed height | walls | result |
|---|---|---:|---:|---|
| `<-2phi>+E8+<2>` | 2,4 | 30 | 17 | 86 finite vertices, 0 cusps, 196 hyperideal rays |
| `<-2phi>+E8+<2>` | 2,4 | 300 | 53 | exact hyperideal escape; full ray census deliberately stopped |
| `<-2phi>+E8+<2>` | 2,4 | 1000 | 215 | exact hyperideal escape; no full ray census requested |
| `<-2phi>+E8+<1>` | 1,2 | 30 | 17 | 26 finite vertices, 0 cusps, 149 hyperideal rays |
| `<-2phi>+E8+<1>` | 1,2 | 300 | 184 | exact hyperideal escape; no full ray census requested |

For the 53-wall cone the escape has exact norm `68/9 - (14/3)phi > 0`.
Numerical linear programming suggested its nine active facets; a separate
exact kernel reconstruction verified the norm and **every** wall inequality.
The numerical solver does not decide the verdict. The certificate is
`golden_h9_full_norms_h300_20261013_escape.json`.
The odd-host 184-wall escape has norm `(61-37phi)/49 > 0` and is recorded
in `golden_h9_odd_a1_full_norms_h300_20261013_escape.json`. A standalone
verifier using only Python `Fraction` arithmetic, importing neither the
search nor any project field/matrix module, independently passes both
certificates. It also checks integral reflection coefficients, chronological
wall heights, and a strict timelike interior point.

The natural host's final continuation through height 1000 completed in
772.923 seconds. Its 215-wall cone still contains an exact spacelike ray of
norm `233/2 - 72phi > 0`, independently verified by the same standalone
Fraction checker. The complete root prefix and escape certificate are
`golden_h9_full_norms_h1000_20261013.json` and
`golden_h9_full_norms_h1000_20261013_escape.json`. No higher bound was run.

The odd final summand is an independent nearby host. Its complete candidate
norm alphabet is `{1,2}`: `q | 4` as above; a norm-four integral reflection
would force E8 and final coordinates even, and its norm equation then forces
the time coordinate even as well, contradicting primitivity.

These results prove that the displayed finite prefixes have infinite
volume. They do not prove global nonreflectivity of either host and do not
exclude a later finite-volume completion.

## Independent path-host screen

A second search examined rank-ten tridiagonal lattices with diagonal
`(2c,2,...,2)`, first off-diagonal `-c`, and all other consecutive
off-diagonals `-1`. The field is again `Q(sqrt(5))`. Primitive mirror-norm
candidates were obtained by exhaustively listing canonical totally positive
divisors of `2 det(B)`. The initial A9 chamber is only a seed here; no global
height-zero completeness claim is made for this family.

| c | height | norm candidates | walls | finite vertices | hyperideal rays |
|---|---:|---:|---:|---:|---:|
| `1+phi` | 30 | 4 | 21 | 17 | 636 |
| `2+phi` | 15 | 12 | 12 | 6 | 42 |
| `2+2phi` | 30 | 8 | 19 | 4 | 283 |
| `3+2phi` | 30 | 8 | 30 | 3 | 964 |
| `4+3phi` | 30 | 8 | 34 | 6 | 824 |

Every row has zero cusps and zero opposite-sheet rays. The corresponding
dimension-four `c=1+phi` control recovers the known compact `[5,3,3,3]`
simplex with five facets and five finite vertices. The naive height-100
path runs and the height-50 `c=2+phi` full cone census were interrupted;
neither is recorded as an exclusion or a completed bound.

## Reproduction and validation

```
PYTHONPATH=src:. python3 scripts/search_golden_h9_full_norms_20261013.py \
  --height 300 --skip-geometry \
  --output results/golden_h9_full_norms_h300_20261013.json
PYTHONPATH=src:. python3 scripts/find_golden_h9_prefix_escape_20261013.py \
  results/golden_h9_full_norms_h300_20261013.json
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_golden_h9_mixed_norms_20261013.py \
  tests/test_og10_golden_h9_vinberg_through_3_plus_6phi.py
```

The existing exact norm-two 71-wall replay passed unchanged (one test).
Four additional tests passed, checking the new norm-four wall and its actual
reflection matrix, the old norm-two default and height convention, rejection
of unsupported norms, and an actual integral norm-one wall in the odd host.
