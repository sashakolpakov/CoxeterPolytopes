# Hostile audit: the `q=2` Hermitian cubic and classical `P^6`

Status date: 2026-08-09.

## Verdict

**Green; no mathematical defect found.**  A second implementation, importing
neither the construction certificate nor the project Vinberg verifier, rebuilds
the entire chain

```text
F4 Frobenius graphs
  -> five transverse general-position centers
  -> smooth Hermitian cubic and all 27 (-1)-curves
  -> literal ERT 6+15+6 wall normals
  -> 72 finite + 27 ideal + 0 hyperideal rays
  -> all 648 Vinberg edge and all 135 affine-neighborhood tests.
```

The result is exactly the classical Potyagailo--Vinberg right-angled `P^6`.
It is a useful positive surface-to-polytope factory calibration, but it is not a
new abstract polytope, a new cubic surface, an arbitrary-polytope inverse
theorem, or an `L4` pushback result.  The `q-1 | 2` statement excludes `q>3`
only from this full-Néron--Severi Möbius-reflection factory.

## Independent finite geometry

Work in

\[
 \mathbf F_4=\mathbf F_2[w]/(w^2+w+1).
\]

On `P1 x P1`, the homogeneous equations

\[
 y_0x_1^2+y_1x_0^2=0,
 \qquad
 x_0y_1^2+x_1y_0^2=0
\]

define the two Frobenius graphs.  Their five displayed intersections are

\[
 ((s:t),(s^2:t^2)),\qquad (s:t)\in\mathbf P^1(\mathbf F_4).
\]

Their divisor product is `(2,1).(1,2)=5`.  On both the finite chart and the
reciprocal chart, the local equations are `b+a^2` and `a+b^2`; in
characteristic two their Jacobian is

\[
 \begin{pmatrix}0&1\\1&0\end{pmatrix}
\]

and has determinant one.  Thus the five intersections are transverse, and
Bézout proves that they exhaust the geometric intersection.

Independent `F4` row reduction gives

```text
all ten triples in H^0(O(1,1))^*:       rank 3
all five four-subsets:                   rank 4.
```

Each triple kernel has a nonsingular `2 x 2` bilinear coefficient matrix and
vanishes at exactly its assigned three centers.  Both ruling projections are
injective.  These are precisely the standard degree-three-del-Pezzo
general-position conditions for five points on `P1 x P1`.

## Independent anticanonical model

The five-center evaluation matrix on `H^0(O(2,2))` has rank five and kernel
dimension four.  Multiplying the four independently recovered sections and
solving for relations gives

```text
degree 1 relation dimension: 0
degree 2 relation dimension: 0
degree 3 relation dimension: 1.
```

The cubic generator is literally

\[
 z_0^2z_3+z_0z_3^2+z_1^2z_2+z_1z_2^2=0.
\]

Its Hermitian coefficient matrix is the rank-four permutation matrix pairing
`z0` with `z3` and `z1` with `z2`.  Its gradient is

\[
 (z_3^2,z_2^2,z_1^2,z_0^2),
\]

which has no projective zero over any extension field.  Hence the cubic is
geometrically smooth.  As a separate finite check, direct enumeration finds
`85` points in `P3(F4)` and exactly `45` points on this cubic, agreeing with
the blowup count `25+5*4=45`.

## The 27 curves and the integral marking

In the basis `F,G,E1,...,E5`, with `F.G=1` and `E_i^2=-1`, the audit constructs

| species | count | class |
|---|---:|---|
| exceptionals | 5 | `E_i` |
| ruling fibres | 10 | `F-E_i`, `G-E_i` |
| projectivity graphs | 10 | `F+G-E_i-E_j-E_k` |
| Frobenius graphs | 2 | `F+2G-sum E_i`, `2F+G-sum E_i` |

All have square `-1` and canonical pairing `-1`.  They are geometrically
realized: the first two species are immediate, the ten projectivity graphs
were checked above, and the last two are strict transforms of smooth rational
graphs separated by the five transverse blowups.

For exhaustiveness, put

\[
 H=3F+G-\sum E_i,
 \quad L_i=F-E_i\ (1\leq i\leq5),
 \quad L_6=2F+G-\sum E_i.
\]

The seven columns `H,L1,...,L6` have determinant `+-1` and Gram matrix
`diag(1,-1^6)`, while `K=-3H+sum L_i`.  If

\[
 D=dH-\sum m_iL_i,
 \qquad D^2=K.D=-1,
\]

then

\[
 \sum m_i^2=d^2+1,
 \qquad \sum m_i=3d-1.
\]

Cauchy forces `d=0,1,2`; exhaustive integer enumeration gives `6+15+6`.
The resulting 27 vectors equal the 27 displayed curves exactly.  Under

\[
 dH-\sum m_iL_i\longmapsto(m_1,\ldots,m_6,d),
\]

their normal multiset is literally

```text
(-e_i,0)                                  6
(e_i+e_j,1), i<j                         15
(sum_(j!=i)e_j,2)                         6,
```

the normal list in ERT Table 2.  This is stronger than a graph or face-count
match.

## Exact cone and literal finite-volume test

The hostile audit sends those 27 outward normals directly to
`cddexec_gmp`, in ERT's Lorentz form `spatial^2-time^2`.  It replays every
rational generator against every inequality, primitive-normalizes it, and
classifies it by exact Lorentz square:

| kind | count | square | active walls |
|---|---:|---:|---:|
| finite | 72 | `-1` | 6 each |
| ideal | 27 | `0` | 10 each |
| hyperideal | 0 | -- | -- |

Every wall is active.  The timelike vector `(1,1,1,1,1,1,3)` has square
`-3` and strict slack one at every wall.

The independent right-angled Vinberg audit then enumerates:

```text
elliptic rank-5 edge diagrams:                     648
elliptic rank-6 finite vertex diagrams:              72
connected affine two-wall components:               135
maximal rank-5 cusp supports:                         27
edge endpoints (finite, ideal)=(1,1):                432
edge endpoints (finite, ideal)=(0,2):                216
bad finite-volume edges:                               0.
```

For every parallel wall pair, the common zero-neighborhood has exactly eight
walls, split into four further parallel pairs.  Together they form ten walls
in five affine two-node components, of total parabolic rank five.  This gives
`135` successful full zero-neighborhood audits.  These connected affine
components are exhaustive: the Coxeter Gram matrix is `I-A`, and any connected
graph on at least three vertices contains a `P3` or `K3`, whose adjacency
spectral radius is greater than one, so its Gram matrix cannot be positive
semidefinite.

As a cross-check between logically separate calculations, the 72 active-wall
sets from finite cdd rays equal all 72 elliptic rank-six diagrams, and the 27
active-wall sets from ideal rays equal all 27 full parabolic cusp supports.
Thus the exact verdict is

```text
finite volume, noncompact.
```

This explicitly rules out the common error of inferring finite volume merely
from a few ideal endpoints while overlooking hyperideal rays or incomplete
cusp neighborhoods.

## The full-NS factory dichotomy

For general prime-power `q`, a Möbius circle has class

\[
 C_A=F+G-\sum_{i\in A}E_i,
 \qquad |A|=q+1,
\]

so

\[
 C_A^2=1-q=-(q-1).
\]

It has divisibility one in the full Néron--Severi lattice, witnessed by
`C_A.E_i=1` for any incident exceptional.  Its reflection is

\[
 s_C(x)=x+\frac{2(x.C_A)}{q-1}C_A.
\]

If `q-1` divides `2`, this preserves the full integral lattice.  Conversely,
apply the reflection to an incident `E_i`; since `C_A` is primitive,
integrality forces `2/(q-1)` to be integral.  Therefore

```text
full-NS integral circle reflection  <=>  q-1 divides 2
                                     <=>  q=2 or q=3.
```

This is a necessary-and-sufficient theorem for this factory only.  It says
nothing against abstract `q>3` polytopes, other lattices, other surfaces, or
other wall classes.

## Primary lineage and novelty boundary

The lineage claims were checked conservatively against primary sources:

- [Everitt--Ratcliffe--Tschantz](https://arxiv.org/abs/1009.3231), Theorem
  2.1 and Tables 2--3: the exact normal vectors, congruence-two chamber,
  `72+27` vertex census, and `2_21` duality.
- [Potyagailo--Vinberg](https://doi.org/10.4171/CMH/4): the finite-volume
  right-angled `P^n` series.
- [Hirschfeld](https://doi.org/10.1017/S144678870002276X): the classical
  `F4` double-six and its 27-line cubic surface.
- [Coxeter](https://doi.org/10.2307/2371466): `2_21` and the classical
  27-line correspondence.

Accordingly, the abstract polytope, cubic, and cubic-line/Gosset
correspondence are classical.  Worldwide priority for the particular
Frobenius-graph presentation is not claimed.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/audit_mobius_q2_cubic_p6_hostile.py --check-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_mobius_q2_cubic_p6_hostile_audit.py
```

Frozen hostile artifacts:

- `scripts/audit_mobius_q2_cubic_p6_hostile.py`
- `data/enumerations/mobius_q2_cubic_p6_hostile_audit.json`
- `tests/test_mobius_q2_cubic_p6_hostile_audit.py`
- `results/mobius_q2_cubic_p6_hostile_audit.md`
