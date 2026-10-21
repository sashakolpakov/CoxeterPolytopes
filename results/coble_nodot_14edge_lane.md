# No finite-volume chamber in the degree-three 14-edge nodot lane

**Status (2026-08-09): exact exhaustive negative theorem.  No compact or
finite-volume noncompact survivor; no surface or K3 gate is reached.**

The executable certificate is
[`certify_coble_nodot_14edge_lane.py`](../scripts/certify_coble_nodot_14edge_lane.py),
with archive
[`coble_nodot_14edge_lane.json`](../data/enumerations/coble_nodot_14edge_lane.json)
and regression tests
[`test_coble_nodot_14edge_lane.py`](../tests/test_coble_nodot_14edge_lane.py).
The frozen archive SHA-256 is

```text
e1a14bfed652de62979baf61afaabf262bd96482e0ec7e3376994267884a03aa
```

## 1. Declared denominator

Take twelve norm-`-2` proposed Coble--Mukai walls.  Require that their
ordinary nonorthogonality support is

1. connected;
2. has exactly fourteen edges;
3. has maximum ordinary degree at most three; and
4. has pairings in `{0,1,2}`, with zero, one, two, or three of the fourteen
   nonzero pairs equal to two.

Every weighted presentation with exact Gram inertia `(1,9,2)` is retained.
Unlike the preceding thirteen-edge theorem, no determinant is selected in
advance.  For every saturated radical quotient the certifier enumerates
**all** even unimodular finite-index extensions using isotropic subgroups of
the complete discriminant quadratic form.

The lane is broader than the determinant-four thirteen-edge slice, but it is
not all fourteen-edge nodot data: ordinary degree above three and four or
more pairing-two pairs remain outside the quantifier.  It is therefore not
the full `D_1_split^nodot` denominator.

## 2. Exhaustive combinatorial and integral census

Nauty's unlabelled graph enumeration and exact arithmetic give:

| stage | count |
|---|---:|
| connected unlabelled fourteen-edge supports | 5,291 |
| weighted trials | 2,486,770 |
| exact labelled Lorentzian corank-two presentations | 236 |
| weighted-diagram isometry classes | 194 |
| classes admitting an even unimodular `E10` extension | 148 |
| distinct integral `E10` extension branches | 169 |

The weighted trial and exact-survivor counts by number of pairing-two pairs
are:

| pairing-two pairs | weighted trials | exact Lorentzian presentations |
|---:|---:|---:|
| 0 | 5,291 | 2 |
| 1 | 74,074 | 12 |
| 2 | 481,481 | 90 |
| 3 | 1,925,924 | 132 |

For an even quotient lattice `L`, even overlattices correspond to isotropic
subgroups of its discriminant quadratic form.  A unimodular extension has
index `k` exactly when

\[
                         |\det L|=k^2.
\]

The certifier enumerates every subgroup of order `k`, keeps precisely the
quadratically isotropic ones, constructs an HNF basis for the corresponding
overlattice, and rechecks integrality, evenness, determinant `-1`, signature
`(1,9)`, and all twelve root coordinates.  Thus every retained extension is
the even unimodular Lorentzian lattice `E10`.

| root-quotient determinant | `E10` branches | extension index |
|---:|---:|---:|
| `-1` | 75 | 1 |
| `-4` | 63 | 2 |
| `-9` | 6 | 3 |
| `-16` | 10 | 4 |
| `-36` | 12 | 6 |
| `-64` | 3 | 8 |
| **total** | **169** | |

All nonsquare determinant classes, and square determinant classes whose
discriminant form has no maximal isotropic subgroup, are rejected before
hyperbolic geometry.

## 3. Complete volume result

Each of the 169 integral branches is sent to GMP cddlib in an integral
`E10` basis.  The complete homogeneous V-representation is replayed with
exact arithmetic: every ray satisfies every wall inequality, every Lorentz
square is computed exactly, all input-facet incidences are recorded, and a
primitive hyperideal witness is frozen for every failed branch.

> **Theorem 3.1 (fourteen-edge gap theorem).**  Every integral `E10` branch
> in the declared fourteen-edge lane has a hyperideal extreme ray.  Hence
> all 169 chambers have infinite volume.  There are no compact and no
> finite-volume noncompact Coxeter polytopes in this lane.

The number of hyperideal extreme rays ranges from **5** to **56**.  Thus no
case is rejected merely because a local parabolic audit was inconclusive:
each has a direct exact escape direction in the complete polyhedral cone.

The certifier's literal Vinberg policy is deliberately ordered:

1. reject immediately if complete cdd exhibits a hyperideal endpoint;
2. on every no-hyperideal branch, enumerate all elliptic rank-eight edges
   and verify their two endpoints;
3. enumerate every connected affine component and its complete zero
   neighbourhood; and
4. verify every possibly disconnected maximal cusp decomposition.

There are zero branches reaching steps 2--4, because all 169 fail step 1.
Consequently there is no finite-volume survivor for which an endpoint or
affine-neighbourhood check is missing.

The complete ledger is therefore:

| geometric status | branches |
|---|---:|
| compact | **0** |
| finite-volume noncompact | **0** |
| infinite-volume, explicit hyperideal ray | **169** |

## 4. Polytope-to-surface meaning

This theorem is upstream of algebraic geometry.  A `D_1_split` carrier must
first supply a finite-volume Coble--Mukai root chamber.  Since no branch in
this lane does so, none receives a cusp lift, anti-Gram split census,
equivariant K3 embedding, period, or Coble return attempt.

That stopping rule matters: failure here is a theorem about the raw
polyhedral input, not a failure to find a surface.  In the primal/dual
programme this deletes a complete fourteen-edge band from the example
search while leaving the surface compiler untouched.

No new hyperbolic polytope or surface is claimed.  The result is a scoped
negative classification, not strict L4 and not an all-nodot theorem.  The
remaining adjacent lanes are:

- fourteen ordinary edges with degree at least four;
- fourteen ordinary edges with at least four pairing-two pairs; and
- the fifteen-edge/no-parallel lane containing Type II.

## Replay

```bash
PYTHONPATH=src:. python3 scripts/certify_coble_nodot_14edge_lane.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_coble_nodot_14edge_lane.py
```
