# Silver Coxeter--Todd-squared host through the first cutter

## Outcome

For

\[
L=\langle-(12\sqrt2-13)\rangle\oplus K_{12}\oplus K_{12}
\subset \mathbf H^{24},
\]

the old `z=(1,(9/5)h,0)` certificate first becomes cut at

\[
H=393 + 278*sqrt(2).
\]

An explicit primitive norm-two wall at that height is accepted by every
earlier packet and pairs with the old witness by
`8/5 - 1/5*sqrt(2)>0`.  Replacing it by

\[
z'=(1,(7/4)h,0)
\]

gives the strict hyperideal norm
`205/12 - 12*sqrt(2)>0` and survives every possible wall through
the same height.  Thus the complete bounded prefix still has infinite volume;
this is **not** a compact or finite-volume example and is not a termination
claim.

## Complete chronology through the bound

The 16 certified mirror-norm classes give exactly
**128 controls at 87 heights**.
There are 76 controls strictly
between the 613,872-wall first packet and the cutting height.  Only the
following necessary Gram rows survive:

| height | norm | Gram types | directional pairs | max `h.y` | old margin | replacement margin |
|:---|:---|---:|---:|:---|:---|:---|
| `67 + 48*sqrt(2)` | `1` | 1 | 5 | `6 + 4*sqrt(2)` | `9/5 - 14/5*sqrt(2)` | `3/2 - 3*sqrt(2)` |
| `67 + 48*sqrt(2)` | `2` | 2 | 10 | `8 + 6*sqrt(2)` | `-28/5 + 9/5*sqrt(2)` | `-6 + 3/2*sqrt(2)` |
| `356 + 252*sqrt(2)` | `6 + 2*sqrt(2)` | 2 | 42 | `36 + 26*sqrt(2)` | `-56/5 - 56/5*sqrt(2)` | `-13 - 25/2*sqrt(2)` |
| `362 + 256*sqrt(2)` | `10 - 4*sqrt(2)` | 1 | 19 | `30 + 22*sqrt(2)` | `-4 + 8/5*sqrt(2)` | `-11/2 + 1/2*sqrt(2)` |
| `393 + 278*sqrt(2)` | `1` | 4 | 140 | `15 + 11*sqrt(2)` | `-2 + 4/5*sqrt(2)` | `-11/4 + 1/4*sqrt(2)` |
| `393 + 278*sqrt(2)` | `2` | 11 | 715 | `22 + 16*sqrt(2)` | `8/5 - 1/5*sqrt(2)` | `1/2 - 1*sqrt(2)` |

The exact fixed-coordinate formula

\[
\min\{x^2:x_1=k\}=\frac34k^2+\rho(k\bmod4),\qquad
\rho=(0,13/4,1,13/4),
\]

is replayed by PARI on the four base cosets.  It turns the apparent earlier
`H=356+252sqrt(2)` cutter into the forbidden combination of norm eight and
dual coordinate three.  The old witness therefore remains strict at both
intermediate necessary heights.  At `H=393+278sqrt(2)`, the norm-one row is
still strict, while the norm-two row has positive old margin.

## Explicit accepted cutter

Set `D=U-V` and `E=2V-U`.  The displayed vectors in the JSON certificate obey

\[
\operatorname{Gram}(D,E)=
\begin{pmatrix}28&46\\46&84\end{pmatrix},\qquad
\operatorname{Gram}(U,V)=(380,204,278),
\]

with `(h.U,h.V)=(22,16)`.  Hence

\[
r=(10+7\sqrt2;\ U+\sqrt2V,0),\qquad r^2=2.
\]

Exact `qfminim` shell maximization gives `36+26sqrt(2)` against the norm-four
shell and `40+30sqrt(2)` against the norm-at-most-eight shell.  Both are
strictly below the relevant first-packet thresholds.

All three possible intermediate binary Gram types are handled without
enumerating their walls: an exhaustive integral Schur-complement census of
the four-vector Gram matrix gives maximum cutter pairing `0, 0`, and
`-7-4sqrt(2)`.  The four simultaneous norm-one types are strictly negative.
For simultaneous norm-two roots, equal height and the two real Cauchy bounds
leave positive inner products only `1,sqrt(2),2`; the first two contradict the
even minimum-four trace equation and the last identifies the same root.

## Replacement and bounded conclusion

The largest necessary pairing at the cutting norm-two row is `22+16sqrt(2)`.
Thus the replacement margin is

\[
\frac74(22+16\sqrt2)-(38+29\sqrt2)=\frac12-\sqrt2<0.
\]

Every other row in the table is also strict.  Since the basepoint is strict,
the segment to the spacelike `z'` crosses the light cone away from every facet,
so the chamber prefix through `H=393+278sqrt(2)` has an open ideal patch and
infinite volume.  Later walls may still cut this replacement.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_coxeter_todd_squared_h393_replacement.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_h393_replacement.py
```
