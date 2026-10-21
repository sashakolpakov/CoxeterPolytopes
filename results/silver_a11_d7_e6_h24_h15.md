# Silver-unit `N(A11 D7 E6)` through two coupled packets

## Outcome

The exact search of

\[
 \langle-(1+\sqrt2)\rangle\oplus N(A_{11}D_7E_6)
 \subset\mathbf H^{24}
\]

does not produce a compact or finite-volume prefix.  The second nonempty
packet is genuinely glue-rich: 22 of its 35 walls lie outside the root
lattice, and 12 couple all three irreducible factors.  Nevertheless a
strict hyperideal `A11`-axis ray survives all 38 positive-height walls.
This is a complete bounded enumeration through height 15, not a global
nonreflectivity result.

## Cyclic index-twelve glue

The root discriminant group is

\[
 (\mathbf Z/12)\oplus(\mathbf Z/4)\oplus(\mathbf Z/3).
\]

Using compatible generators, the Niemeier glue is cyclic:

\[
 \mathcal G=\langle(1,1,1)\rangle,
 \qquad |\mathcal G|=12.
\]

An explicit generator is

\[
 g=(\omega^{A_{11}}_1,\omega^{D_7}_6,\omega^{E_6}_1),
 \qquad g^2=\frac{11}{12}+\frac74+\frac43=4.
\]

For a triple of Dynkin weights, membership requires the `D7` and `E6`
residues to equal the reductions modulo four and three of the `A11`
residue.  This is imposed separately on both coefficients in
`Z[sqrt(2)]`.  A 24-vector glue basis has integral even determinant-one Gram
matrix.  Height zero has exactly three norm-two dominant orbits, hence
precisely `A11+D7+E6` and 288 roots.

Even unimodularity and the unit coefficient give the complete primitive
mirror-norm list `1,2` by the exact dyadic pairing-ideal argument.

## Exact chronology and geometry

The packet census through real height 15 is:

| height | candidates | rejected | accepted | nontrivial glue | all-three-factor glue |
|:---|---:|---:|---:|---:|---:|
| `(1+sqrt(2))/2` | 0 | 0 | 0 | 0 | 0 |
| `1+sqrt(2)` | 0 | 0 | 0 | 0 | 0 |
| `2+2sqrt(2)` | 3 | 0 | 3 | 0 | 0 |
| `(7+5sqrt(2))/2` | 0 | 0 | 0 | 0 | 0 |
| `7+5sqrt(2)` | 35 | 0 | 35 | 22 | 12 |

Every same-height pairing is nonpositive.  The first packet consists of the
three factorwise affine caps.  Its product has

\[
 12\cdot8\cdot7=672
\]

rays: 3 finite, no ideal, and 669 hyperideal.

## Strict witness

Set time coordinate one and use only the first fundamental `A11` direction,
with Dynkin weight `5/3`.  Since

\[
 (C_{A_{11}}^{-1})_{11}=\frac{11}{12},
\]

the Lorentz norm is

\[
 \frac{25}{9}\frac{11}{12}-(1+\sqrt2)
 =\frac{167}{108}-\sqrt2>0.
\]

Its pairings with all 38 positive-height walls are strictly negative; the
closest equals

\[
 -\frac13-\frac13\sqrt2.
\]

The height-zero Dynkin inequalities also hold, so the ray is strictly
hyperideal and the prefix has infinite volume.

## Reproduction

- Backend: `src/coxeter_polytopes/search/silver_a11_d7_e6.py`
- Audit: `scripts/audit_silver_a11_d7_e6_h24.py`
- Tests: `tests/test_silver_a11_d7_e6.py`
- Certificate: `data/enumerations/silver_a11_d7_e6_h24_h15.json`

The focused suite passes `3/3` tests in 99.28 seconds.  Certificate SHA-256:
`992fcbe04b88ed33d0f3042aff95f550c1613dc492d98c7e7b345a7f0dd7aa6c`.

