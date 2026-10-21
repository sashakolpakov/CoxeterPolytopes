# Exchanged-root obstruction to invariant-nef cone completeness

Status: **proved general criterion; exact `D34/T1` application** (2026-08-08).

## Result

An involution-exchanged pair of effective `(-2)` classes contributes one
unavoidable inequality on the invariant divisor space.  If its invariant
sum takes both signs on a proposed Coxeter cone, that cone cannot be the
invariant section of a nef cone.  This is an algebraic-geometric obstruction
to obligation `C` (cone completeness); it is not a failure of the weak
lattice carrier and says nothing by itself about hyperbolic finite volume.

The newly compiled `D34/T1` rank-18 K3 package fails this test.  An explicit
hidden orthogonal root pair cuts its proposed fourteen-facet cone.  Thus the
known `T1 H9/14` Coxeter polytope has this new weak K3 lattice/action carrier,
but this particular carrier cannot realize that cone as its exact invariant
nef section.

## The criterion

Let `X` be a smooth projective surface, let

\[
N=\operatorname{NS}(X)/\text{torsion},\qquad
F=N_{\mathbf R}^{\tau},
\]

where `tau` is an involutive isometry induced by an automorphism, or more
generally just a lattice involution used to define the slice.  Let
`C` be a full-dimensional convex polyhedral cone in `F` on the chosen
positive-cone sheet.  Suppose there are primitive classes
`r_+,r_- in N` such that

\[
r_+^2=r_-^2=-2,\qquad \tau(r_+)=r_-,\qquad
f:=r_++r_-\in N^\tau,
\]

and, for one sign `epsilon in {+1,-1}`, both `epsilon r_+` and
`epsilon r_-` are effective.  If `C` contains two vectors `v_+` and `v_-`
with

\[
(f,v_+)>0>(f,v_-),
\]

then

\[
C\ne \operatorname{Nef}(X)\cap F.
\]

Indeed, for every `x in F`, invariance of the intersection form gives

\[
(x,r_+)=(x,r_-)=\tfrac12(x,f).
\]

Nefness and effectivity therefore force

\[
\epsilon(x,f)\ge 0
\quad\text{for every }x\in\operatorname{Nef}(X)\cap F.
\]

The two signs on `C` show that neither possible oriented half-space contains
`C`.  If the effective orientation `epsilon` is already known, one ray with
`epsilon(f,v)<0` is enough.  The mixed-sign certificate is stronger because
it rules out both orientations at once.

No surface automorphism is needed for this main statement: the proof uses
only the lattice exchange, fixedness of `x`, and independent effectivity of
the two classes.  A geometric action is merely one convenient way to deduce
the paired effectivity automatically.

For a projective K3 surface the effectivity hypothesis is automatic once
`tau` is geometric.  Riemann--Roch gives effectivity of exactly one of
`r_+` and `-r_+`; preservation of effectivity gives the same sign for
`r_- = tau(r_+)`.  This automatic orientation is **not** available for a
bare abstract lattice involution unless its geometric/effective-cone action
has also been proved.

Because `C` is full-dimensional, mixed signs also imply that `f^perp` meets
the relative interior of `C`.  To see this, choose `h in relint(C)`.  If
`(f,h)` is nonzero, add a positive multiple of a cone vector having the
opposite sign; `relint(C)+C` is contained in `relint(C)`, and the multiple
can be chosen to make the pairing zero.  Thus this is a genuine hidden cut,
not merely contact with a boundary ray.  If a lifted effective root is
reducible, the conclusion still holds: its effective components impose the
needed nef inequalities.  Irreducibility is required only to call the new
hyperplane an irredundant geometric wall.

## Exact `D34/T1` witness

Use the fixed basis of the `T1` compiler.  Its fixed lattice has rank ten and
signature `(1,9)`.  The hidden invariant vector is

```text
f  = (0,0,0,1,0,0,-1,0,0,0),       f^2 = -4.
```

Two exact extreme rays of the candidate fourteen-facet cone are

```text
v- = (3,1,2,5,4,6,0,0,0,-4),       v-^2 = 10,
v+ = (5,3,6,4,4,2,1,0,0,-4),       v+^2 = 10.
```

They lie on the same positive-cone sheet because `(v-,v+)=36`.  Their
fourteen facet-slack vectors are respectively

```text
(0,0,10,0,2,0,6,0,8,0,0,4,0,0),
(0,0,0,6,0,10,0,2,8,0,0,4,0,0).
```

All slacks are nonnegative, and the active facet normals have rank nine in
both cases, so these are genuine extreme rays.  Their hidden pairings are

```text
(f,v-) = -6,                         (f,v+) = +6.
```

In the archived rank-18 lattice basis, `f` lifts to the two vectors

```text
r+ = ( 1, 1, 1,1,0,0, 0,0,0,0,-1,-2,-2,0,0,0,0,-1),
r- = (-1,-1,-1,0,0,0,-1,0,0,0, 1, 2, 2,0,0,0,0, 1).
```

They are integral and primitive, satisfy

```text
r+^2 = r-^2 = -2,    (r+,r-) = 0,    tau(r+) = r-,    r+ + r- = E f,
```

and are not among the 24 displayed roots.  On the invariant lifts of the two
rays their individual pairings are

```text
                  r+       r-
v-                -3       -3
v+                +3       +3.
```

Thus whichever of the pairs `{r+,r-}` or `{-r+,-r-}` is effective on a K3
realization, one of the two candidate rays is not nef.  The obstruction is
transported by every marked isometry carrying the complete triple
`(N,tau,C)`; an isometry of `N` alone, without the action and cone marking,
does not carry this conclusion.

The production certificate is
[`h9_14_t1_strong_k3_nef_obstruction.json`](../data/enumerations/h9_14_t1_strong_k3_nef_obstruction.json),
with SHA-256
`4eb338784ec452cb08fb4a75b6f4be4c43155d1e6fce2a2259bbe39f410037c1`,
generated by
[`certify_h9_14_t1_strong_k3_nef_obstruction.py`](../scripts/certify_h9_14_t1_strong_k3_nef_obstruction.py)
and tested by
[`test_h9_14_t1_strong_k3_nef_obstruction.py`](../tests/test_h9_14_t1_strong_k3_nef_obstruction.py).
The case-specific discussion is
[`h9_14_t1_strong_k3_nef_obstruction.md`](h9_14_t1_strong_k3_nef_obstruction.md).

## Minimal computational certificate schema

An exact failure certificate needs only the following data.

1. An integral Gram matrix `G_N`, an integral involution `T` with
   `T^2=1` and `T^t G_N T=G_N`, and a fixed-space embedding `E` satisfying
   `TE=E`.
2. The candidate facet matrix and two rational or integral vectors `v_+,
   v_-`.  Verify all facet slacks are nonnegative, that both vectors are on
   the selected nonnegative-square sheet, and optionally that active-normal
   rank is `dim(F)-1` when they are advertised as extreme rays.
3. Integral primitive vectors `r_+,r_-` with square `-2`, with
   `Tr_+=r_-`, and with `r_++r_-=Ef` for an invariant vector `f`.
4. The two exact inequalities `(f,v_+)>0>(f,v_-)` (equivalently, the
   individual root pairings after applying `E`).
5. A geometric effectivity certificate: explicit effective divisors, or the
   K3 Riemann--Roch argument together with proof that `T` is induced by an
   effectivity-preserving geometric action.

Completeness of the entire extreme-ray list is unnecessary for a failure:
two verified members of the cone suffice.  It is necessary when a one-sign
scan is used to claim that no such obstruction exists.  A strict interior
point `x_0` with positive old-facet slacks and `(f,x_0)=0` is a convenient
optional strengthening, but follows abstractly from full dimensionality and
the mixed signs.

## Consequences for the transfer obligations

- **`U` (universality):** the weak lattice/period carrier may still exist.
  What fails is any strong forward-transfer claim that this compiled package
  automatically produces the prescribed invariant-nef cone.  A different
  package or a different kind of variety is not excluded.
- **`R` (marked recovery):** the obstruction does not attack recovery in
  principle.  For this package `R0` (the exact labelled fold) survives, but
  the production audit also proves that no literal equivariant bridge to the
  original 24-wall source package exists; full geometric marked `R` is not
  claimed.  The obstruction depends on transporting `(N,tau,C)`, not merely
  on matching a determinant, discriminant form, or unmarked lattice.
- **`C` (cone completeness):** closed negatively for this `T1` package.  The
  hidden exchanged pair imposes an extra nef inequality and cuts the
  fourteen-facet cone.
- **`V` (hyperbolic closure):** untouched.  The original labelled `T1`
  system remains the already certified finite-volume noncompact Coxeter
  polytope.  If one refines the cone by the hidden inequality, the refined
  wall system requires a new, independent Vinberg audit.

Accordingly, the strong algebraic stratum demanding exact invariant-nef
equality is empty for this marked package, while its weak K3 period stratum
remains nonempty.  This separation is precisely why `C` must be kept distinct
from primitive embedding, period existence, marked recovery, and finite
hyperbolic volume.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_t1_strong_k3_nef_obstruction.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_t1_strong_k3_nef_obstruction.py
```
