# Fixed lattices of the characteristic-three high-dimensional survivors

## Result

The complete 5,944-wall Fermat-quartic stock has four casts which pass exact
global irredundancy, the full Coxeter-angle gate and a complete finite-volume
test.  All four are finite-volume, noncompact Coxeter polytopes.  Their
saturated common fixed lattices have now been computed directly in the
published rank-22 Neron--Severi lattice by Smith kernels.

| cast | dimension | determinant | Smith diagonal |
|---|---:|---:|---|
| cyclic `C2-H11` | `H11` | 1,024 | `1^2 2^10` |
| cyclic `C2-H10` | `H10` | 6,144 | `2^10 6` |
| noncyclic `V4-H11` | `H11` | 768 | `1^6 2^4 4 12` |
| noncyclic `V4-H10` | `H10` | 1,024 | `1^2 2^8 4` |

Every displayed lattice is even and has signature `(1,n)` for the stated
hyperbolic dimension.  The certificate records an explicit integral basis in
the published Neron--Severi basis and the complete integral Gram matrix for
each case.

## Exact identification of the positive 123-wall H10 case

The `V4-H10` 123-wall chamber has 37,451 finite vertices and 1,732 cusps;
literal Vinberg and independent GMP cone enumeration both prove finite
volume and noncompactness.

Its Smith fingerprint agrees with the earlier Jacobian--Kummer fixed lattice,
but the identification is stronger than that.  Exact discriminant-form
reduction gives

\[
A_L\cong(\mathbf Z/2)^8\oplus\mathbf Z/4,
\qquad
q_L\cong q_+(\mathbf F_2^8)\oplus\langle-1/4\rangle.
\]

This is exactly the discriminant form of

\[
U\oplus E_8(-2)\oplus\langle-4\rangle.
\]

Both lattices are even indefinite of signature `(1,10)`, and

\[
\operatorname{rank}L=11=\ell(A_L)+2.
\]

Nikulin's indefinite uniqueness criterion therefore proves an integral
isometry

\[
L\cong U\oplus E_8(-2)\oplus\langle-4\rangle.
\]

The certificate additionally records an explicit unimodular matrix `P`, of
determinant one, and verifies entry by entry that

\[
P^{\mathsf T}G_{V4-H10}P
=G_{U+E_8(-2)+\langle-4\rangle}.
\]

Thus the integral identity does not depend only on invoking uniqueness; the
coordinate isometry is available for the subsequent subgroup comparison.

Under this explicit isometry, the 123 facet mirrors have primitive intrinsic
types

| square and divisibility | mirrors |
|---|---:|
| `(-2,1)` | 10 |
| `(-4,2)` | 31 |
| `(-4,4)` | 82 |

Every corresponding reflection matrix is integral, involutive, and preserves
the target Gram.  Therefore `W_123` is literally a reflection subgroup of
the known maximal group, rather than merely an abstract Coxeter group with a
commensurable quadratic form.

Consequently the full integral reflection group has Vinberg's known 13-wall
`f_2^10` chamber.  The new construction is a finite-index reflection-subgroup
chamber in that known group, not a new lattice, maximal reflection group, or
arithmetic commensurability class.

Solomon's complete sum over all finite special subgroups of the 123-wall
diagram gives

\[
\chi=-\frac{527}{2880},\qquad
[W(f_2^{10}):W_{123}]=65{,}802{,}240,
\qquad
\operatorname{vol}(P_{123})=\pi^5\frac{527}{85050}.
\]

The index is exactly four times the index of the 84-wall Jacobian--Kummer
subgroup.  This numerical relation does not prove containment; an explicit
lattice isometry and subgroup comparison are the next exact lineage gate.

## Exact identification of the positive 196-wall H10 case

The cyclic `C2-H10` chamber has independently passed both volume criteria:
literal Vinberg gives 212,824 finite vertices, 10,560 cusps, and 1,427,720
rank-nine edges, while exact GMP double description gives the same finite and
ideal counts and zero hyperideal rays.

Every entry of its fixed Gram `G` is divisible by two.  The certificate
records a determinant-one integral matrix `P` and verifies

\[
P^{\mathsf T}(G/2)P=\operatorname{diag}(3,-1,\ldots,-1).
\]

Thus the half-scaled lattice is exactly the determinant-three lattice behind
McLeod's `f_3^10` chamber.  Multiplying a form by two changes neither its
integral orthogonal group nor the reflection formula, so the ambient maximal
reflection group is McLeod's known 14-wall group.

In this half-scaled form, the 196 primitive mirrors have types

| square and divisibility | mirrors |
|---|---:|
| `(-1,1)` | 80 |
| `(-2,1)` | 36 |
| `(-6,3)` | 80 |

Every corresponding reflection matrix is integral, involutive, and preserves
the form.  Hence `W_196` is literally a reflection subgroup of
`W(f_3^10)`, not merely a group in the same rational commensurability class.

The complete Solomon sums give

\[
\chi(W_{196})=-\frac{5797}{3780},\qquad
\chi(W(f_3^{10}))=-\frac{341}{16721510400},
\]

and therefore

\[
[W(f_3^{10}):W_{196}]=75{,}202{,}560
=2^{15}3^3\cdot5\cdot17,
\qquad
\operatorname{vol}(P_{196})=\pi^5\frac{46376}{893025}.
\]

This closes the lattice, maximal-group, and commensurability lineage.  What
remains open is priority of the literal 196-wall finite-index subgroup
diagram and its divisor-theoretic realization.

## Mirror lattices of the two positive H11 cases

The cyclic `C2-H11` chamber has 206 facets, 327,566 finite vertices and
18,354 cusps.  A complete literal Vinberg audit proves finite volume and
noncompactness.  All 206 reflections act integrally on its saturated fixed
divisor lattice, whose 2-elementary invariants identify it as

`U + E8(-2) + A1(-2)^2`.

The noncyclic `V4-H11` chamber has 328 facets, 225,408 finite vertices and
11,554 cusps.  Exact double description gives no hyperideal extreme ray, so
it too is finite-volume and noncompact.  On its saturated fixed divisor
lattice, exactly 36 of the 328 reflections have denominator two.  The exact
one-step reflection closure is nevertheless a minimal stable even integral
overlattice of index four, signature `(1,11)`, determinant 48 and Smith
diagonal `1^8,2^3,6`; every facet reflection is integral on it.  Thus the
328-wall group is arithmetic over `Q`.  An explicit unimodular isometry
identifies the stable lattice with Scharlau--Kirschmer List11/Gram11 entry
20, whose full integral reflection chamber has 22 walls; hence the induced
328-wall group is a finite-index reflection subgroup of that known group.
The full mirror-type certificate and exact lineage are recorded in
[`ss_fermat_char3_h11_mirror_lattices.md`](ss_fermat_char3_h11_mirror_lattices.md)
and
[`ss_fermat_char3_v4_h11_lineage.md`](ss_fermat_char3_v4_h11_lineage.md).

## Method and reproducibility

For each subgroup, the fixed-row equations are stacked as
`g_i^T-I`.  A full Smith decomposition gives the saturated integral kernel.
LLL only shortens its basis; an explicit unimodular basis-change matrix is
checked before the fixed Gram is formed.  The computation then verifies
pointwise fixedness, signature, parity, determinant, and Smith invariants.

- certificate:
  [`ss_fermat_quartic_char3_survivor_fixed_lattices.json`](../data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json)
- generator:
  [`certify_ss_fermat_char3_survivor_fixed_lattices.py`](../scripts/certify_ss_fermat_char3_survivor_fixed_lattices.py)
- tests:
  [`test_ss_fermat_char3_survivor_fixed_lattices.py`](../tests/test_ss_fermat_char3_survivor_fixed_lattices.py)
- 196-wall McLeod index certificate:
  [`ss_fermat_quartic_char3_c2_h10_mcleod_index.json`](../data/enumerations/ss_fermat_quartic_char3_c2_h10_mcleod_index.json)
- 196-wall index reproducer and tests:
  [`certify_ss_fermat_char3_c2_h10_mcleod_index.py`](../scripts/certify_ss_fermat_char3_c2_h10_mcleod_index.py),
  [`test_ss_fermat_char3_c2_h10_mcleod_index.py`](../tests/test_ss_fermat_char3_c2_h10_mcleod_index.py)
- H11 mirror and stable-overlattice certificate:
  [`ss_fermat_quartic_char3_h11_mirror_lattices.json`](../data/enumerations/ss_fermat_quartic_char3_h11_mirror_lattices.json)
- H11 mirror reproducer and tests:
  [`certify_ss_fermat_char3_h11_mirror_lattices.py`](../scripts/certify_ss_fermat_char3_h11_mirror_lattices.py),
  [`test_ss_fermat_char3_h11_mirror_lattices.py`](../tests/test_ss_fermat_char3_h11_mirror_lattices.py)

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_ss_fermat_char3_survivor_fixed_lattices.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_ss_fermat_char3_survivor_fixed_lattices.py
```
