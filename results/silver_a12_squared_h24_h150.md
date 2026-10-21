# Silver-unit `N(A12^2)` through height 150

## Outcome

The exact search of

\[
 \langle-(1+\sqrt2)\rangle\oplus N(A_{12}^2)
 \quad\text{in }\mathbf H^{24}
\]

does not produce a compact or finite-volume prefix.  Through real height
150, the complete primitive chronology has 280 accepted positive-height
walls.  An explicit ray remains strictly hyperideal against every one.
This is a bounded negative screen, not a global nonreflectivity theorem.

## Prime-13 glue

Each `A12` factor has discriminant group `Z/13`.  The Niemeier extension is
the graph

\[
 \{(k,8k):k\in\mathbf Z/13\},
\]

equivalently the inverse graph with multiplier five.  The implemented glue
generator is

\[
 g=(\omega_5,\omega_1),\qquad
 g^2=\frac{5\cdot8}{13}+\frac{1\cdot12}{13}=4.
\]

For Dynkin weights `z`, the residue is

\[
 \rho(z)=\sum_{i=1}^{12}i z_i\pmod {13}.
\]

The condition `rho(right)=8 rho(left)` is imposed independently on the
rational and `sqrt(2)` coefficients.  An explicit basis consisting of all
12 simple roots in the first factor, 11 in the second, and `g` has integral
even Gram matrix of determinant one.  The exact height-zero census has two
norm-two dominant roots and therefore exactly `A12^2`, with 312 roots.

## Complete local norm list and chronology

Even unimodularity and the unit coefficient give the complete canonical
primitive mirror norms `1,2`.  The usual pairing-ideal argument at the
ramified prime `(sqrt(2))` proves no other norm ideal is missing, and
division by `sqrt(2)` is tested in the full graph-glued lattice.

The nonempty packets through height 150 are:

| height | candidates | rejected | accepted | prime-13 glue accepted |
|:---|---:|---:|---:|---:|
| `2+2sqrt(2)` | 2 | 0 | 2 | 0 |
| `7+5sqrt(2)` | 4 | 0 | 4 | 0 |
| `14+10sqrt(2)` | 24 | 4 | 20 | 0 |
| `28+20sqrt(2)` | 40 | 10 | 30 | 0 |
| `41+29sqrt(2)` | 180 | 16 | 164 | 0 |
| `56+40sqrt(2)` | 48 | 20 | 28 | 0 |
| `63+45sqrt(2)` | 44 | 12 | 32 | 0 |

The intervening control shells are also enumerated and empty.  Remarkably,
none of the 280 accepted walls has nonzero prime-13 residue: the full glue
decoder is active, but the chronology stays in the root lattice this far.

The first two walls are the factorwise affine `A12` caps.  Their product
has `13^2=169` rays: one finite, no ideal, and 168 hyperideal.

## Strict deep witness

Take time coordinate one, zero weights in the second factor, and the single
first-factor Dynkin weight

\[
 z_1=\frac53.
\]

Since `(C_A12^-1)_{11}=12/13`, its Lorentz norm is

\[
 \frac{25}{9}\frac{12}{13}-(1+\sqrt2)
 =\frac{61}{39}-\sqrt2>0.
\]

All height-zero inequalities are nonpositive.  Exact evaluation on all 280
positive-height walls is strictly negative; the closest pairing is

\[
 -\frac13\sqrt2<0.
\]

Thus the complete height-150 prefix has a strict hyperideal ray and is not
finite-volume.

## Reproduction

- Backend: `src/coxeter_polytopes/search/silver_a12_squared.py`
- Audit: `scripts/audit_silver_a12_squared_h24.py`
- Tests: `tests/test_silver_a12_squared.py`
- Certificate: `data/enumerations/silver_a12_squared_h24_h150.json`

```console
PYTHONPATH=src python3 scripts/audit_silver_a12_squared_h24.py
PYTHONPATH=src pytest -q tests/test_silver_a12_squared.py
```

The focused suite passes `3/3` tests in 147.41 seconds.  The certificate
SHA-256 is
`f5f1125385b1aa11bb1229fd1e094e604c879b4334a8919d0033675974ad3867`.

