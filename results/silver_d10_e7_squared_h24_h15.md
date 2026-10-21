# Silver-unit `N(D10 E7^2)` exact binary-glue screen

No compact or finite-volume prefix is found for

\[
 \langle-(1+\sqrt2)\rangle\oplus N(D_{10}E_7^2)
 \subset\mathbf H^{24}.
\]

The root discriminant group is `(Z/2)^4`.  With `D10` classes
`0,v,s,c`, the full index-four Niemeier glue is

\[
 \{(0,0,0),(s,1,0),(c,0,1),(v,1,1)\}.
\]

Twice an explicit HNF basis is stored in the backend.  Its Gram matrix is
integral, even, and determinant one.  The height-zero census has exactly
three norm-two dominant orbits, hence root system `D10+E7+E7` with 432
roots.  Even unimodularity gives the complete primitive mirror norms `1,2`.

Through height 15 the two nonempty packets are:

| height | candidates | rejected | accepted | glue | all three factors |
|:---|---:|---:|---:|---:|---:|
| `2+2sqrt(2)` | 3 | 0 | 3 | 0 | 0 |
| `7+5sqrt(2)` | 16 | 0 | 16 | 7 | 3 |

All same-height pairings are nonpositive.  The first affine product has
704 rays: 4 finite, no ideal, and 700 hyperideal.

The ray with time coordinate one and only the first `D10` Dynkin weight
`5/3` has norm

\[
 \frac{25}{9}-(1+\sqrt2)=\frac{16}{9}-\sqrt2>0.
\]

It pairs strictly negatively with all 19 positive-height walls; the closest
pairing is `-1/3-(1/3)sqrt(2)`.  Thus the complete two-packet prefix has
infinite volume.  This is a bounded result, not global nonreflectivity.

Reproduce with

```console
PYTHONPATH=src python3 scripts/audit_silver_d10_e7_squared_h24.py
PYTHONPATH=src pytest -q tests/test_silver_d10_e7_squared.py
```

The focused suite passes `3/3` tests in 13.60 seconds.  Certificate SHA-256:
`ad1a793a2350c0518c6140a07b1660a646db2f194f58d044d994853d14acf9d2`.

