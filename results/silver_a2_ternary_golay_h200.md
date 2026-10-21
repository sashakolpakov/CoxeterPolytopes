# Silver `N(A2^12)` necessary frontier through real height 200

## Result

The exact ternary-Golay block search for

\[
 \langle-\sqrt2\rangle\oplus N(A_2^{12}),\qquad
 \langle-2\sqrt2\rangle\oplus N(A_2^{12})
\]

has now been carried through every arithmetic control of real height at most
200.  After the affine packet, exactly eight levels have any vectors obeying
the affine caps and both ternary glue conditions:

| height | complete candidate count |
|---|---:|
| \(8+6\sqrt2\) | 3,960 |
| \(24+17\sqrt2\) | 4,733,256 |
| \(41+29\sqrt2\) | 1,724,736 |
| \(48+34\sqrt2\) | 720,156,144 |
| \(60+43\sqrt2\) | 27,720 |
| \(63+45\sqrt2\) | 27,720 |
| \(80+57\sqrt2\) | 45,875,104 |
| \(96+68\sqrt2\) | 103,188,357,648 |

Every other control in the cutoff is exactly empty.  These are complete
first-cap **necessary supersets**; later chronological walls can only delete
vectors from them.

The first two nonempty packets were treated explicitly.  All 3,960 walls at
the first level are mutually acute.  The next level contains a certified
primitive accepted wall, and its 264-element zero-glue orbit is mutually
acute.  Full chronological enumeration of the much larger later candidate
sets is unnecessary for the finite-volume obstruction below.

## Uniform obstruction

For the coefficient `a=sqrt(2)`, put Dynkin weight \(c=1/4\) in both
coordinates of every `A2` block.  The Lorentz norm is

\[
 24c^2-\sqrt2=\frac32-\sqrt2>0.
\]

For `a=2sqrt(2)`, use \(c=7/20\), with norm

\[
 24c^2-2\sqrt2=\frac{147}{50}-2\sqrt2>0.
\]

The contingency dynamic program records the exact maximum of the total
Dynkin-weight sum at every control.  Substitution shows that both points pair
strictly negatively with every vector in all eight complete candidate
supersets, as well as every affine and height-zero wall.  The closest
pairings still have negative exact sign; no floating-point comparison is
used.

Consequently both complete bounded prefixes contain an open hyperideal patch
and have infinite volume.  This does **not** prove eventual nonreflectivity:
walls above height 200 may cut these points.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_a2_ternary_golay_h200.py
PYTHONPATH=src pytest -q tests/test_silver_a2_ternary_golay.py
```

The exact archive is
`data/enumerations/silver_a2_ternary_golay_h200.json`.
