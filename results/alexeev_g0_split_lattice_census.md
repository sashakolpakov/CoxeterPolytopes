# A polytope-first split-root obstruction for the nine `r+a=22` chambers

**Status (2026-08-09): exact non-tautological carrier theorem, with a sharp
Type-VI/Type-MII calibration.**  Starting only from Alexeev's nine complete
mixed-reflection wall Gram matrices with `r+a=22`, the split-root condition
excludes **all nine** as exact full-nef strong-Coble chambers over the complex
numbers.  Five fail before any rank bound is used.  The `(16,6,1)` row fails
over every characteristic different from two.  At the positive endpoint,
the `(18,4,0)` chamber itself forces the K3 host to be supersingular in
characteristic three with Artin invariant one.

This is genuine `POLYTOPE -> SURFACE` leverage, but it is a theorem about a
precisely decorated carrier class.  It neither destroys Alexeev's abstract
arithmetic polytopes nor proves that every passing lattice completion is
effective.

## 1. The surface-free input and the exact hypothesis

Let `P_(r,a,delta)` be one of Alexeev's **complete** fundamental chambers for
the mixed reflection group of the even 2-elementary lattice `M` of type
`(r,a,delta)`, where

\[
                         r+a=22.
\]

The raw datum is its primitive wall Gram matrix.  Its leading `r` walls give
an integral basis of `M`; this is checked, not assumed.  The compiler uses no
Kondō type, surface labels, curve list, or known Coble realization.

The decoration is the canonical one forced by an exact strong-Coble quotient
interpretation.  If

\[
 M=\widetilde{\operatorname{NS}}(S)(2),
\]

then a primitive wall has the following meaning.

* A long wall `v^2=-4` of divisibility two must be a downstairs irreducible
  `(-2)`-curve.
* A short wall of square `-2` can instead encode a boundary half or a
  downstairs `(-1)`-curve.

Only the first family is forced through the split-root test.  In particular,
the argument does **not** silently split a stable pullback of a `(-1)`-curve,
and it does not choose between the transverse and tangent alternatives for a
`(-1)`-curve with boundary profile `(2)`.

The theorem therefore has the exact hypothesis

```text
the complete mixed-reflection chamber P, with every long facet interpreted
as an actual downstairs (-2)-curve, is literally the full nef chamber.
```

This is independently checkable from the marked polytope.  It is stronger
than merely asking for the same invariant lattice.

## 2. The forced finite Gram-completion problem

Let `pi:X -> S` be the canonical K3 double cover, in characteristic different
from two.  A downstairs `(-2)`-curve `C_i` is disjoint from the branch divisor.
Since a tame étale double cover of `P1` is trivial,

\[
 \pi^*C_i=c_i^++c_i^-,\qquad
 (c_i^\pm)^2=-2,\qquad c_i^+c_i^-=0,
\]

where the involution exchanges the two curves.  Put

\[
       v_i=c_i^++c_i^-,\qquad d_i=c_i^+-c_i^-.
\]

Then `d_i` is anti-invariant, `d_i^2=-4`, and `d_i.M=0`.  For two distinct
long walls set

\[
 n_{ij}={v_i v_j\over2}=C_iC_j.
\]

There are nonnegative integers

\[
 a_{ij}=c_i^+c_j^+=c_i^-c_j^-,\qquad
 b_{ij}=c_i^+c_j^-=c_i^-c_j^+
\]

with `a_ij+b_ij=n_ij`.  Hence, for

\[
                         H=-\frac12\operatorname{Gram}(d_i),
\]

one must have

\[
 H_{ii}=2,\qquad
 H_{ij}=n_{ij}-2a_{ij}
       \in\{-n_{ij},-n_{ij}+2,\ldots,n_{ij}\}.          \tag{2.1}
\]

The `d_i` lie in the negative-definite algebraic anti-invariant space, so

\[
                              H\succeq0.                 \tag{2.2}
\]

Equations (2.1)--(2.2) are a finite condition on `P` alone.  Here every raw
long-long product is `0`, `2`, or `4`, so a pair has respectively one, two,
or three possible lift allocations.

There is also a saturation condition.  If an integral vector `z` lies in
`ker(H)`, then `sum z_i d_i=0`.  Since the fixed lattice `M` is primitive,

\[
                    \sum_i z_i v_i\in2M.                 \tag{2.3}
\]

The certifier checks (2.3) on a saturated integral basis of every radical.
Thus a singular numerical Gram is not accepted merely because its real
eigenvalues have the right signs.

Over the complex numbers, the algebraic anti-invariant rank is at most

\[
                         20-r.                           \tag{2.4}
\]

This is the bridge from the wall combinatorics to the K3 Picard-rank bound.

## 3. Exact census

| `(r,a,delta)` | `dim P` | walls / long | split completion | minimum anti-invariant rank | complex budget `20-r` | generated `(rank, |det|)` | consequence |
|---|---:|---:|---|---:|---:|---:|---|
| `(11,11,1)` | 10 | `12 / 10` | none | -- | 9 | -- | no split carrier in any `char !=2` |
| `(12,10,1)` | 11 | `14 / 11` | none | -- | 8 | -- | no split carrier in any `char !=2` |
| `(13,9,1)` | 12 | `19 / 10` | none | -- | 7 | -- | no split carrier in any `char !=2` |
| `(14,8,0)` | 13 | `17 / 10` | none | -- | 6 | -- | no split carrier in any `char !=2` |
| `(14,8,1)` | 13 | `24 / 10` | none | -- | 6 | -- | no split carrier in any `char !=2` |
| `(15,7,1)` | 14 | `60 / 10` | `E6(-2)` differences | 6 | 5 | `(21,6)` | complex no-go; positive characteristic required |
| `(16,6,1)` | 15 | `118 / 25` | `(A5+A1)(-2)` differences | 6 | 4 | `(22,12)` | no K3 carrier in any `char !=2` |
| `(18,4,0)` | 17 | `48 / 24` | `(A2+A2)(-2)` differences | 4 | 2 | `(22,9)` | host forced to `p=3`, Artin invariant 1 |
| `(20,2,1)` | 19 | `50 / 5` | `A1(-2)` differences | 1 | 0 | `(21,4)` | complex no-go; positive-characteristic gate remains |

The pair `(rank,|det|)` in the seventh column belongs to the even lattice
generated by the displayed **minimum-rank** completion.  For `(16,6,1)` and
`(18,4,0)` the K3 rank ceiling forces every possible carrier to attain that
minimum, so their determinant conclusions are universal, not merely
properties of a convenient completion.

### The first five rows

Their `n_ij=1` support graph is a forest.  Swapping `c_i^+` and `c_i^-`
changes the sign of row and column `i` of `H`; therefore every signing on a
forest is diagonally congruent to the all-minus representative `2I-A`.
The exact representative inertias and determinants are

| row | inertia `(+, -, 0)` | determinant |
|---|---:|---:|
| `(11,11,1)` | `(9,1,0)` | `-1` |
| `(12,10,1)` | `(10,1,0)` | `-8` |
| `(13,9,1)` | `(9,1,0)` | `-16` |
| `(14,8,0)` | `(9,1,0)` | `-16` |
| `(14,8,1)` | `(9,1,0)` | `-64` |

Thus no signing is positive semidefinite.  This obstruction precedes K3
embedding, periods, characteristic, and effectivity.

### Type VI: the Petersen completion

The ten long walls form the Petersen support graph.  After fixing the signs
on a spanning tree, six cycle signs remain.  The certifier checks all `2^6=64`
switching classes by exact rational congruence inertia.  Exactly one is
positive semidefinite; it has inertia `(6,0,4)`.  Its saturated root lattice
`H` is `E6`, so the differences span `E6(-2)`.

A fixed six-wall principal submatrix is three orthogonal `A2` blocks and is
positive definite.  Consequently rank six is a lower bound independent of
the signing, and the surviving completion attains it.  Every integral
radical relation satisfies (2.3).

This explains, directly from the polytope, why the known characteristic-three
Type-VI carrier cannot move to a complex K3: `r+6=21>20`.  The split lattice
alone does not yet determine the remaining supersingular Picard class or
prove its effectivity, so it does not by itself recover characteristic three.

### The rank-six, rank-four, and rank-one completions

For `(16,6,1)`, the first six long walls contain a forced positive-definite
`A5+A1` block of determinant `12`.  Its nonzero-choice graph is a forest, so
all of its signings switch to the displayed block.  At rank six it spans the
whole anti-invariant difference space.  For each of the remaining nineteen
walls the certifier enumerates every allowed pairing column `p` and solves

\[
                         p^tB^{-1}p=2.
\]

There are exactly two integral solutions, one the negative of the other.
They are precisely the two choices of sheet.  All independent sign choices
remain pairwise compatible.  Hence there is one rank-minimal switching orbit,
with root lattice `A5+A1`.

The same exhaustive argument for `(18,4,0)` starts from two orthogonal `A2`
blocks, determinant `9`.  Every one of the remaining twenty walls again has
exactly one opposite pair of integral solutions.  Thus the minimum is four,
and the unique rank-minimal switching orbit has root lattice `A2+A2`.

For `(20,2,1)`, one `A1` root is already a lower bound and all five difference
classes can be chosen on that line, giving minimum one.

These are exhaustive minimum-rank statements.  The forest normalization
removes every basis signing; once the positive basis spans, pairings with it
determine every other vector, and the finite norm-two search leaves only
independent sheet swaps.

## 4. The glued even lattice and the two strongest consequences

Let `A` be the negative difference lattice, with Gram `-2H`, and define

\[
                 K=\left\langle M,A,{v_i+d_i\over2}\right\rangle. \tag{4.1}
\]

The long walls have divisibility two.  Equations (2.1) and (2.3) therefore
make (4.1) an even integral lattice.  If `s` is the binary rank of the glue
classes `(v_i,d_i) mod 2`, then

\[
 [K:M\oplus A]=2^s,\qquad
 |\det K|={2^a\,2^{\operatorname{rk}A}|\det H|\over2^{2s}}. \tag{4.2}
\]

In all four feasible minimum completions, `s=rk(A)`.  Formula (4.2) gives
the determinants `6,12,9,4` displayed in the table.

### `(18,4,0)` forces the known supersingular host

Here `K` already has rank `22`, signature `(1,21)`, and determinant `-9`.
Any finite-index even overclosure of index `q` has determinant `-9/q^2`.
The only nontrivial numerical possibility is `q=3`, which would be even
unimodular of signature `(1,21)`.  Such a lattice cannot exist because the
signature difference is not divisible by eight.  Thus `K` is maximal even,
and any K3 Néron--Severi lattice containing the split configuration equals
`K`.

A rank-22 K3 surface is supersingular, and its Néron--Severi determinant is
`-p^(2 sigma)`.  Therefore

\[
                         p=3,\qquad \sigma=1.            \tag{4.3}
\]

This recovers from the raw `H17/48` polytope exactly the host class of the
known Type-MII cover, the supersingular `Kum(E x E)` in characteristic three
with Artin invariant one.  The result is necessary, not yet sufficient: an
involution preserving the crystalline period, the ramification frame,
effectivity, and exact projected-nef equality remain separate gates.

### `(16,6,1)` is impossible in every characteristic

Here every K3 carrier would also have to use the minimum rank, because
`16+6=22`.  The generated lattice has determinant `-12`.  For any finite-index
even rank-22 overclosure, integrality of `12/q^2` leaves only

```text
index 1: determinant -12;
index 2: determinant -3.
```

Neither absolute determinant is a prime to an even power.  Hence neither is
the Néron--Severi lattice of a supersingular K3.  Rank 22 is impossible for a
complex or finite-height K3, so the decorated `(16,6,1)` chamber has no K3
cover in any characteristic different from two.

## 5. Why this does not contradict complex Type-I and Type-II surfaces

The `(11,11,1)` no-go concerns Alexeev's complete **12-wall** mixed-reflection
chamber and requires all ten of its long facets to be actual downstairs
`(-2)` walls of that same full nef cone.

The complex one-boundary Type-I and Type-II Coble surfaces have the same
rank-11 invariant lattice type, but their complete full-Néron--Severi nef
chambers are different wall systems: respectively `H10/22` and `H10/30`.
Each contains twelve actual `(-2)` curves together with boundary and
exceptional walls.  They are subgroup/full-nef chambers, not the 12-wall
Alexeev mixed-reflection fundamental chamber.  The invariant lattice alone
does not identify a nef chamber; that distinction is precisely what the
split-root theorem detects.

The same warning applies generally:

```text
same 2-elementary lattice != same effective roots != same nef polytope.
```

## 6. What has and has not been gained

The result supplies a non-singleton, independently defined, surface-free
class on which a polytope calculation yields new algebraic-geometric
restrictions:

* a class-wide complex no-go for exact full-nef Coble realization;
* five characteristic-independent numerical no-gos;
* one all-characteristic K3-lattice no-go at `(16,6,1)`;
* a forced characteristic and Artin invariant at `(18,4,0)`;
* an explanation of the positive-characteristic Picard jump at Type VI.

It remains a **necessary carrier criterion**, not a necessary-and-sufficient
surface compiler.  Passing it does not prove an involution, period, effective
root system, hidden-wall exclusion, or exact return of `P`.  Failure excludes
this Coble/K3-cover realization, not the abstract arithmetic polytope or a
different carrier species.

## 7. Reproducibility

The exact certificate is
[`alexeev_g0_split_lattice_census.json`](../data/enumerations/alexeev_g0_split_lattice_census.json),
generated by
[`certify_alexeev_g0_split_lattice_census.py`](../scripts/certify_alexeev_g0_split_lattice_census.py)
and tested by
[`test_alexeev_g0_split_lattice_census.py`](../tests/test_alexeev_g0_split_lattice_census.py).

The pinned raw source is Alexeev's arXiv v4 ancillary file
`all_cases_in_2el_paper.txt`, expected at

```text
/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt
```

with SHA-256

```text
4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457
```

The certifier refuses to run if the hash changes.  Replay with

```bash
PYTHONPATH=src:. python3 scripts/certify_alexeev_g0_split_lattice_census.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_alexeev_g0_split_lattice_census.py
```

The broader K3 quotient and effectivity gates are separated in
[`coble_k3_cover_polytope_to_surface_bridge.md`](coble_k3_cover_polytope_to_surface_bridge.md),
while the distinct complex Type-I/II full-nef chambers are recorded in
[`coble_one_boundary_full_ns_h10_theorem.md`](coble_one_boundary_full_ns_h10_theorem.md).

Primary mathematical inputs are Alexeev's *Reflective hyperbolic
2-elementary lattices, K3 surfaces and hyperkähler manifolds*
(arXiv:2209.09110v4), the elementary geometry of the tame canonical double
cover, the complex K3 Picard bound, and the Rudakov--Shafarevich form of the
supersingular K3 lattice theorem.
