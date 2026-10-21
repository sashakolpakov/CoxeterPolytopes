# A Coxeter polytope determines a Looijenga-pair period problem

**Status (2026-08-08):** proved conditional inverse-carrier theorem for a
fixed lattice-combinatorial Looijenga decoration, plus a complete worked
eight-facet arithmetic `H4` example.  The worked example has an associated
one-dimensional marked Looijenga-pair stratum `C* \ {1}` and its Coxeter cone
is exactly the full nef cone.

This is the direction

```text
arithmetic Coxeter datum P
    -> finitely many Looijenga decorations extracted from P
    -> period subtori with hidden-root hypertori deleted
    -> marked Looijenga pairs
    -> exact full-nef promotion test.
```

It is not a construction of a chamber starting from a favorite surface.  In
the worked case the eight wall vectors recover the Picard lattice, canonical
class, anticanonical cycle, toric blowdown, period lattice, and surface family.

## 1. What is proved, and what is not

Arithmeticity by itself does **not** currently imply that a Coxeter polytope
has a Looijenga-pair carrier.  A direct carrier has extra integral data:

1. a rational-surface Picard saturation;
2. a characteristic canonical vector;
3. an anticanonical cycle with a realizable deformation type;
4. a period point and a choice of nef/Weyl chamber; and
5. equality of the proposed cone with the complete nef cone.

For a bounded class these enhancements can be enumerated from `P`.  For each
fixed enhancement the existence criterion below is necessary and sufficient.
Consequently, after an exhaustive enhancement enumeration, nonemptiness of
the finite union of period strata is necessary and sufficient for `P` to have
a carrier of the prescribed kind.

This direct statement concerns rational (`Q`-form) arithmetic data.  A
totally-real form of degree greater than one would normally need a trace
lattice and a distinguished eigenspace of larger Picard rank.  Classical
Looijenga-pair Torelli theory does not by itself provide that eigencarrier
promotion.

## 2. The lattice-combinatorial Looijenga decoration

Let `P` be an oriented rational arithmetic Coxeter polytope in `H^n`, with
primitive wall vectors `r_j` in a Lorentz lattice.  A **direct Looijenga
decoration** consists of the following data, all in the real span of `P`.

### 2.1 Rational-surface lattice

There is a full-rank integral saturation

```text
wall lattice M(P)  subset  L  subset  M(P)^dual
```

with

```text
L = I_(1,n),                 rho(L)=n+1
```

for `rho >= 3`.  Equivalently, one enumerates the unimodular overlattices of
the wall lattice and retains the odd Lorentzian one.  A smooth rational
surface of Picard rank `rho` has

```text
K^2 = 10-rho = 9-n.
```

Thus the decoration also contains a characteristic vector `k in L` with

```text
k^2=9-n.
```

### 2.2 Anticanonical cycle

For a cycle of length `m >= 3`, choose classes `d_1,...,d_m in L` such that

```text
sum_i d_i = -k,
d_i.d_(i+1)=1,
d_i.d_j=0 for non-neighbours.
```

Then

```text
d_i^2+k.d_i=-2
```

is the adjunction equation for each rational boundary component.  These
numerical equations alone are not sufficient: the same boundary
self-intersection sequence can have more than one deformation type.
Friedman's Theorem 3.1 proves that there are only finitely many deformation
types for a fixed sequence, not that the sequence uniquely selects one.

The safe constructive certificate is a toric model.  It specifies a smooth
complete toric fan with boundary classes `bar d_i`, mutually orthogonal total
exceptional classes `e_ij` of square `-1`, and

```text
d_i = pullback(bar d_i)-sum_j e_ij.
```

Blowing up the corresponding smooth boundary points constructs `(Y,D)` and
automatically gives `D=sum d_i=-K_Y`.  Gross--Hacking--Keel explain this
construction and prove that every Looijenga pair acquires a toric model after
possibly making corner blowups; see Definition 1.2 and Lemma 1.3 of their
[primary paper](https://arxiv.org/abs/1211.6367).  In the direct, fixed-rank
problem we retain an explicit toric model rather than silently increasing the
Picard rank.

### 2.3 Allowed geometric wall types

Every proposed facet must be assigned one of the following roles:

- an interior exceptional curve `e`, with `e^2=k.e=-1`;
- an internal nodal curve `alpha`, with `alpha^2=-2` and
  `alpha.d_i=0` for all `i`; or
- a negative component of the boundary cycle.

Adjunction shows that these are all possible irreducible negative curves away
from the boundary.  In particular, for primitive internal curve walls the
only finite Coxeter labels are forced by intersection theory:

| wall types | intersection | relation |
|---|---:|---|
| `(-1),(-1)` | `0` | `m=2` |
| `(-1),(-1)` | `1` | parallel (`m=infinity`) |
| `(-2),(-2)` | `0` | `m=2` |
| `(-2),(-2)` | `1` | `m=3` |
| `(-2),(-2)` | `2` | parallel |
| `(-1),(-2)` | `0` | `m=2` |
| `(-1),(-2)` | `1` | `m=4` |

Larger intersections give ultraparallel walls.  Thus an ordinary full-nef
Looijenga realization with only primitive internal curve normals cannot
produce a finite `5`-edge.  A `4`-edge has the particularly rigid meaning
`(-1)` wall meeting a `(-2)` wall once.

## 3. Exact period nonemptiness: the genuine `P -> surface` theorem

Fix a deformation type produced by the decoration and put

```text
Lambda = D^perp = {x in L : x.d_i=0 for every i}.
```

The natural period torus for a Picard-marked pair without marked points on
the boundary is

```text
T' = Hom(Lambda,C*).
```

This is the base of the **rigidified** Picard-marked, boundary-unmarked moduli
problem.  The fully boundary-marked version uses
`T=Hom(L,C*)`; quotienting by `Aut^0(D)=(C*)^m` maps it to `T'`.  The
unrigidified boundary-unmarked stack can retain the subgroup of automorphisms
acting trivially on `Pic(Y)`.  Gross--Hacking--Keel make these distinctions
precise in Theorem 6.1.

Let `C_D++` be their deformation-invariant generic cone: start with the
positive cone, impose all oriented numerical-exceptional inequalities, and
then impose the boundary inequalities.  Friedman's/Gross--Hacking--Keel's
root theorem gives the deformation roots exactly as

```text
Phi = {alpha in Lambda :
       alpha^2=-2 and alpha^perp meets interior(C_D++)}.
```

It is essential that `Phi` need not be every square-minus-two vector of
`Lambda`.

### 3.1 Intended root subsystem and hidden-root deletion

Let `R0 subset Lambda` be the sublattice generated by the intended nodal
root subsystem derived from the Coxeter walls.  The exact-root period stratum
is

```text
T_R0^o = {
  phi in Hom(Lambda,C*) :
  phi|R0=1,
  phi(beta)!=1 for every beta in Phi\R0
}.
```

Equivalently, put `Q=Lambda/R0`.  Then

```text
T_R0^o = {
  chi in Hom(Q,C*) : chi(beta_bar)!=1 for every forbidden root beta
}.
```

This gives an exact nonemptiness criterion:

> `T_R0^o` is nonempty if and only if there is a character of `Q` which is
> nontrivial on every forbidden-root image.

For nonprimitive `R0`, this includes a finite check on characters of the
torsion subgroup of `Q`; individual point-separation is not enough if a
finite union of kernels covers every torsion character.  For primitive `R0`,
`Q` is free.  Every forbidden root has a nonzero image precisely when it is
not already forced by `R0`, and a generic character avoids the finite (or
countable) union of proper root hypertori.  Thus in the primitive case the
criterion reduces to:

```text
no root declared forbidden lies in R0.
```

The proof is elementary character theory, while existence of a Looijenga
pair at every surviving period point is Friedman's period-surjectivity
theorem (Theorem 3.17) and the Gross--Hacking--Keel marked-moduli theorem
(Theorem 6.1).

### 3.2 Why a chamber lift is necessary

The equation `phi(alpha)=1` does **not** by itself say that `alpha` is the
class of an irreducible curve.  It says that `alpha` belongs to the
period-trivial root system.  Gross--Hacking--Keel prove

```text
Phi_phi = W_Y . Delta_Y,
```

where `Delta_Y` is the set of actual irreducible internal `(-2)`-curves.
Their marked moduli is therefore a nonseparated lift of the period torus: a
point is a period `phi` together with a component of

```text
C_D++ \ union_(alpha in Phi_phi) alpha^perp.
```

Choose the branch whose closure has inequalities `x.alpha >= 0` for the
intended simple roots.  The root walls bounding that selected nef chamber
are then the actual irreducible `(-2)`-curves—not merely unspecified Weyl
conjugates.  This is precisely the extra datum missing from a bare period
point.

### 3.3 Independent full-nef gate

Period theory controls only the nodal roots.  It does not decide the
exceptional curves.  The candidate Coxeter cone `C_P` is the complete nef
cone if and only if

```text
C_P = C_D++ intersect (the selected simple-root halfspaces).
```

Equivalently, every effective numerical-exceptional and negative-boundary
inequality must occur among the proposed walls or be a certified nonnegative
combination of them, and there may be no additional simple root wall.
Friedman's Lemma 2.3 gives the complete negative-curve list underlying this
test; Gross--Hacking--Keel Lemma 2.15 gives the same nef-cone formula in their
notation.

Combining the character criterion, chamber branch, and independent
`C_D++` equality proves the promised statement:

> **Fixed-decoration inverse-carrier theorem.** A decorated rational
> arithmetic Coxeter datum has a marked full-nef Looijenga-pair carrier if
> and only if its exact-root character stratum is nonempty, its intended
> root halfspaces bound a chamber in the lifted period space, and its
> Coxeter cone equals the independently computed `C_D++` chamber.  The
> resulting pair is supplied by the marked period theorem.  Conversely any
> marked full-nef carrier supplies all three pieces of data.

This theorem is non-tautological: the period and character gates are
algebraic, while the exceptional-cone equality is deliberately retained as
an independent birational-geometric gate.

## 4. The worked `H4` polytope

### 4.1 Coxeter input

Let

```text
L = ZH + ZE1 + ZE2 + ZE3 + ZE4,
H^2=1, Ei^2=-1, H.Ei=Ei.Ej=0 (i!=j).
```

Starting only with this marked integral Coxeter realization, take the eight
walls

```text
E1, E2, E3, E4,
A  = H-E1-E3-E4,
D1 = H-E1-E2,
B  = H-E2-E3,
C  = H-E2-E4.
```

In the displayed order their Gram matrix is

```text
[-1  0  0  0  1  1  0  0]
[ 0 -1  0  0  0  1  1  1]
[ 0  0 -1  0  1  0  1  0]
[ 0  0  0 -1  1  0  0  1]
[ 1  0  1  1 -2  0  0  0]
[ 1  1  0  0  0 -1  0  0]
[ 0  1  1  0  0  0 -1  0]
[ 0  1  0  1  0  0  0 -1].
```

It has rank `5` and exact inertia `(1,4,3)`.  The five walls
`E1,E2,E3,E4,A` form a unimodular basis, so the wall vectors generate the
full lattice `I_(1,4)`, not merely a hidden finite-index sublattice.

The only non-right finite angles are

```text
A.E1=A.E3=A.E4=1,
```

giving three angles `pi/4`.  Six pairs of square-minus-one walls have inner
product one and are parallel; the remaining nineteen pairs are orthogonal.
Thus the exact Coxeter-pair profile is

```text
19 x m2, 3 x m4, 6 x parallel.
```

### 4.2 The polytope recovers the surface data

The equations

```text
K.r=-1 for every square-minus-one wall r,
K.A=0
```

have the unique solution

```text
K=-3H+E1+E2+E3+E4,
K^2=5.
```

It is characteristic.  A Looijenga decoration is then recovered by

```text
d1 = D1             = H-E1-E2,
d2 = A+E1+E4        = H-E3,
d3 = A+E1+E3        = H-E4.
```

Their Gram matrix and sum are

```text
[d_i.d_j] = [-1 1 1; 1 0 1; 1 1 0],
d1+d2+d3=-K.
```

The toric blowdown is literal:

```text
d1+E1+E2 = H,
d2+E3    = H,
d3+E4    = H.
```

Contracting the four disjoint exceptional curves therefore gives `P2` and
sends the boundary to a triangle of lines.  The blowup distribution on its
three components is `(2,1,1)`.  This proves algebraic existence before any
period specialization.

### 4.3 The soft surface family

Write the four blowup points as

```text
q1,q2 on the first boundary line,
q3 on the second,
q4 on the third.
```

The period lattice is the saturated kernel of the boundary-pairing map:

```text
Lambda=D^perp = ZA + Zdelta,
delta=E1-E2,
[A,delta]=[-2 1; 1 -2] = A2(-1).
```

Hence `T'=Hom(Lambda,C*)=(C*)^2`.  Put

```text
u=phi(A),  t=phi(delta).
```

The six deformation roots are

```text
plus/minus A, plus/minus delta, plus/minus(A+delta).
```

The intended subsystem `R0=ZA` is primitive and
`Lambda/R0=Z delta`.  Its exact-root stratum is therefore

```text
u=1,
t!=1,
T_R0^o = C* \ {1}.
```

Geometrically, `u=1` is the condition that `q1,q3,q4` are collinear.  This
is the same collinearity equation in Gross--Hacking--Keel Example 5.14.  The
extra coordinate `t` measures the difference between the two points on the
first boundary component.  Since the line through `q3,q4` meets that
component at `q1`, the condition `t!=1` simultaneously keeps `q2` distinct
and prevents the second collinearity root.

Over this stratum, `Phi_phi={plus/minus A}`.  The lifted period space has two
branches.  Selecting the branch

```text
x.A >= 0
```

makes `A` itself the effective irreducible curve: it is the strict transform
of the line through `q1,q3,q4`.  The opposite branch makes `-A` effective.
Thus the branch choice eliminates the Weyl-conjugacy ambiguity.

This is the desired softening:

```text
one rigid H4 Coxeter polytope
    <-> a one-dimensional marked Looijenga-pair stratum C*\{1}.
```

Here `C*\{1}` is a branch of the rigidified Picard-marked,
boundary-unmarked moduli.  Adding marked points on all three boundary
components gives the corresponding fully boundary-marked cover; those
extra marking directions are not part of the Coxeter datum.

### 4.4 Complete negative curves and full-nef equality

Every numerical exceptional class has the form

```text
e=aH+sum_i c_i Ei,
e^2=K.e=-1.
```

The two equations give

```text
sum c_i = 1-3a,
sum c_i^2 = a^2+1.
```

Cauchy--Schwarz implies

```text
5a^2-6a-3 <= 0,
```

so `a` is exactly `0` or `1`.  Hence the complete list has ten classes:

```text
E1,E2,E3,E4,
H-Ei-Ej  (1<=i<j<=4).
```

On the exact period stratum the irreducible negative curves are precisely

```text
E1,E2,E3,E4, A, D1, B, C.
```

Indeed, three numerical exceptional classes become the reducible divisors

```text
H-E1-E3 = A+E4,
H-E1-E4 = A+E3,
H-E3-E4 = A+E1.
```

The other three pair classes are `D1,B,C`.  There is no other internal
`(-2)`-curve because the exact period-root set is only `plus/minus A` and the
chosen branch makes `A` the simple effective root.  The remaining boundary
classes give redundant inequalities:

```text
d2=A+E1+E4,
d3=A+E1+E3.
```

Thus every omitted `C_D++` or boundary inequality is a nonnegative sum of
the eight wall inequalities.  Friedman's complete negative-curve theorem
then proves

```text
Nef(Y) = {x in C+ : x.r>=0 for all eight listed walls}.
```

This is full-nef equality, not a divisor subsystem.  Each of the eight
inequalities is irredundant: the certificate supplies a positive-square
point on that wall with strictly positive pairing against the other seven.

### 4.5 Exact finite-volume certificate

The cone is full-dimensional (for example
`20H-E1-E2-E3-E4` is strictly interior) and pointed because its walls span
`L`.  Enumerating every rank-four active-wall subset gives exactly nine
extreme rays:

| ray | square | vertex type |
|---|---:|---|
| `H-E1` | `0` | ideal |
| `H-E2` | `0` | ideal |
| `H-E3` | `0` | ideal |
| `H-E4` | `0` | ideal |
| `H` | `1` | finite |
| `2H-E1-E2-E3` | `1` | finite |
| `2H-E1-E2-E4` | `1` | finite |
| `2H-E2-E3-E4` | `1` | finite |
| `3H-E1-2E2-E3-E4` | `2` | finite |

There is no spacelike extreme ray.  Therefore the projectivized cone has
only finite or ideal vertices, so it is a finite-volume Coxeter polytope.
The four null rays prove it is noncompact.  The exact result is

```text
H4, 8 facets, 9 vertices = 5 finite + 4 ideal,
finite volume, noncompact.
```

### 4.6 Exact arithmetic status

For every wall `r`, reflection is

```text
s_r(x)=x-2(x.r)/r^2 r.
```

The seven square-minus-one walls give coefficient `2`, and the unique
square-minus-two wall gives coefficient `1`.  All eight reflection matrices
are integral involutive isometries of `I_(1,4)` with determinant `-1`.
The Coxeter reflection group is therefore contained in
`O^+(I_(1,4);Z)`.  Since the certified fundamental chamber has finite
volume, the reflection group is a finite-index arithmetic subgroup.  Thus
this is an arithmetic hyperbolic Coxeter polytope in the literal sense, not
only a rational polyhedral nef cone.

No surface automorphism realization of the reflections is asserted or
needed.  In particular, reflection in a `(-1)`-curve does not preserve the
canonical class.

## 5. Smallest bounded completion lane

Gross--Hacking--Keel Example 5.14—the blowup of one point on each side of a
triangle—is the rank-four Picard / `H3` seed.  Its period lattice has rank one,
so imposing its collinearity root fixes the entire unmarked period torus.
Adding the fourth blowup above raises the period rank to two and leaves the
one-dimensional exact-root stratum.  It is therefore the smallest member of
this triangular lane that exhibits the user's intended rigid-to-continuous
correspondence.

A finite completion problem independent of surfaces can now be posed:

```text
C_LP(4,8): rational arithmetic H4 Coxeter data with
  at most 8 facets,
  primitive wall squares in {-1,-2},
  pairings in {0,1},
  wall span I_(1,4),
  and finite volume.
```

For every such `P`:

1. enumerate characteristic `K` vectors satisfying the wall adjunction
   types;
2. enumerate anticanonical cycles and toric blowdown certificates;
3. compute `Lambda`, `C_D++`, and the finite deformation-root set;
4. enumerate intended primitive root sublattices;
5. apply the exact character/hypertorus criterion;
6. enumerate the lifted root chambers;
7. enumerate all numerical exceptional classes and test independent
   `C_D++` equality; and
8. replay Coxeter angles and all extreme rays.

Every step terminates.  More generally, when `K^2=9-n>0`, hence `n<=8`,
`K^perp` is negative definite.  Both the numerical-exceptional equations and
the root enumeration are finite after bounded boundary data are fixed.  This
gives a natural finite positive-Looijenga search through `H8`.

At `n=9`, `K` is isotropic, and for `n>=10` it is spacelike; the relevant
Diophantine sets can be infinite.  Looijenga pairs themselves have no K3-like
rank ceiling, so they remain promising in high dimension, but a high-rank
completion needs a new stopping theorem—e.g. a Cox-ring/Mori-dream criterion,
a finite-orbit theorem for exceptional classes, or an independently bounded
fundamental-domain algorithm.  Period surjectivity alone does not solve that
gate.

## 6. Lineage and novelty discipline

The published seed is Gross--Hacking--Keel Example 5.14: blow up one point
on each coordinate line of `P2`; the root
`H-E1-E2-E3` becomes effective exactly when the three points are collinear,
and the two universal-family charts are related by its reflection.

The present `H4` result adds one free blowup on the first boundary component
and freezes:

- the reverse reconstruction from the eight Coxeter walls;
- the exact `C*\{1}` period stratum;
- the branch which makes `A` itself effective;
- the complete exceptional and nodal curve enumeration;
- full-nef equality;
- all integral reflections; and
- all extreme rays and finite volume.

This is a **project-new P-to-moduli reinterpretation and exact certificate**.
It is not currently claimed as a worldwide-new abstract Coxeter polytope,
nor as a new weak del Pezzo surface.  A dedicated catalogue/lineage search is
still required before making either priority claim.

## 7. Reproducibility

- Exact certifier:
  [`certify_looijenga_h4_inverse_carrier.py`](../scripts/certify_looijenga_h4_inverse_carrier.py)
- Frozen ledger:
  [`looijenga_h4_inverse_carrier.json`](../data/enumerations/looijenga_h4_inverse_carrier.json)
- Regression tests:
  [`test_looijenga_h4_inverse_carrier.py`](../tests/test_looijenga_h4_inverse_carrier.py)

Run from the repository root:

```bash
PYTHONPATH=.:src python3 scripts/certify_looijenga_h4_inverse_carrier.py --check-archive
PYTHONPATH=.:src python3 -m pytest -q tests/test_looijenga_h4_inverse_carrier.py
```

The frozen certificate records `8` facets, `9` vertices, `5` finite vertices,
`4` ideal vertices, and decision `finite_volume_noncompact`.

## 8. Primary sources

1. Mark Gross, Paul Hacking, Sean Keel,
   [*Moduli of surfaces with an anti-canonical cycle*](https://arxiv.org/abs/1211.6367):
   Definition 1.5 (period), Theorem 1.8 (Torelli), Theorem 3.1 (root
   criterion), Proposition 3.4 (period-trivial versus nodal roots), Example
   5.14 (the three-point collinearity seed), and Theorem 6.1 (marked moduli as
   the lifted period torus).
2. Robert Friedman,
   [*On the ample cone of a rational surface with an anticanonical cycle*](https://arxiv.org/abs/1207.7012):
   Theorem 1.3 (period surjectivity), Lemmas 2.3--2.4 (complete negative
   walls and effectivity), and Theorem 2.14 (root criterion).  The
   [publisher PDF](https://msp.org/ant/2013/7-6/ant-v7-n6-p07-s.pdf) is also
   openly available.
3. Robert Friedman,
   [*On the geometry of anticanonical pairs*](https://arxiv.org/abs/1502.02560):
   Theorem 3.1 (finitely many deformation types for fixed boundary sequence),
   Theorem 3.17 (period surjectivity), Theorem 6.6 (root/effectivity
   criterion), and Theorem 6.12 (roots as differences of disjoint exceptional
   curves on a very general deformation).

