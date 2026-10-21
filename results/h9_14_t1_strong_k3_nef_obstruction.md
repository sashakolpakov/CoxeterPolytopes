# TOOT: an exact hidden-root obstruction for the compiled `T1` K3 carrier

Status: **exact carrier-specific negative `C` theorem PASS** (2026-08-08).

## The result

The unique `D34/T1` compiler output remains a valid weak K3 lattice/action
carrier:

\[
N=(18,2,0),\qquad \det N=-4,\qquad
N^\tau=U+E_8(-2),\qquad N_\tau=E_8(-2).
\]

It embeds primitively in the K3 lattice with complement `U+U(2)`, and the
involution extends as the identity on the complement.  The exact period

\[
\omega=(1,2\sqrt6,i\sqrt2,i\sqrt3)
\]

has square zero, positive norm `8*sqrt(6)`, and no nonzero rational vector in
`U+U(2)` orthogonal to it.  Thus a projective K3 surface with exactly
`NS(X)=N` and a genuine symplectic Nikulin involution exists.

However, its labelled fourteen-facet `T1` cone **cannot** be

\[
\operatorname{Nef}(X)\cap N^\tau_{\mathbb R}
\]

for any invariant ample orientation of this compiled package.  Five hidden
square-minus-four restricted walls already occur with fixed coordinates in
`{-1,0,1}`.  Each lifts to an integral orthogonal pair of square-minus-two
roots exchanged by the involution, and each takes opposite signs on two exact
extreme rays of the proposed cone.  The first pair alone proves the theorem.

This is a negative promotion theorem for this particular `(18,2,0)` carrier.
It is not a nonexistence theorem for all possible algebraic realizations of
the abstract `T1` Coxeter polytope.

## First exact hidden wall

In the frozen compiler basis of `F=N^tau`, take

```text
f  = (0,0,0,1,0,0,-1,0,0,0),       f^2 = -4
k  = (-1,-2,-2,-1,0,0,0,-1),        k^2 = -4 in E8(-2).
```

The parity/glue condition is exact.  In the rank-18 `N` basis,

```text
f_N = (0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0)
k_N = (2,2,2,1,0,0,1,0,0,0,-2,-4,-4,0,0,0,0,-2).
```

Hence

```text
r+ = (f_N+k_N)/2
   = (1,1,1,1,0,0,0,0,0,0,-1,-2,-2,0,0,0,0,-1)

r- = (f_N-k_N)/2
   = (-1,-1,-1,0,0,0,-1,0,0,0,1,2,2,0,0,0,0,1).
```

Exact arithmetic gives

\[
r_+^2=r_-^2=-2,\qquad (r_+,r_-)=0,\qquad
\tau(r_+)=r_-,\qquad r_++r_-=f_N.
\]

Neither root, up to sign, is one of the 24 compiled roots.

The proposed cone contains the two finite extreme rays

```text
v- = (3,1,2,5,4,6,0,0,0,-4),   (v-)^2=10,   (f,v-)=-6
v+ = (5,3,6,4,4,2,1,0,0,-4),   (v+)^2=10,   (f,v+)= 6.
```

Their fourteen facet-slack vectors are respectively

```text
(0,0,10,0,2,0,6,0,8,0,0,4,0,0)
(0,0,0,6,0,10,0,2,8,0,0,4,0,0).
```

Each has nine active normals of rank nine, so these are exact extreme rays,
not numerical samples.  Also `(v-,v+)=36`, placing them on the same time
sheet.  On the invariant lifts, each of `r+` and `r-` pairs `-3` with `v-`
and `+3` with `v+`.

For any `tau`-invariant K3 Weyl chamber, Riemann--Roch makes exactly one of
`+/-r+` effective.  Since `tau` exchanges the roots, the same sign of `r-`
is effective.  Therefore the invariant nef section must satisfy one of the
two inequalities `+/-f >= 0`.  The two rays above violate the two possible
orientations in turn.  Irreducibility is unnecessary: if either root is
reducible, its effective decomposition into smooth rational curves still
imposes the same necessary nef inequality.

## Five independent crossing normals

The frozen certificate includes the complete lifted roots and exact rays for
all five witnesses.  Their compact summary is:

| number | fixed normal `f` | `f^2` | negative pairing | positive pairing |
|---:|---|---:|---:|---:|
| 1 | `(0,0,0,1,0,0,-1,0,0,0)` | -4 | -6 | 6 |
| 2 | `(0,0,0,1,1,0,-1,0,0,0)` | -4 | -6 | 6 |
| 3 | `(0,0,1,1,0,0,-1,0,0,0)` | -4 | -4 | 6 |
| 4 | `(0,1,1,1,0,0,-1,0,0,0)` | -4 | -4 | 16 |
| 5 | `(1,0,1,1,0,0,-1,0,0,0)` | -4 | -4 | 16 |

The first bounded audit also examined 4,374 square-minus-two fixed roots and
found no crossing one; the obstruction appears precisely in the exchanged
orthogonal-pair channel that produces square-minus-four orbit sums.

## Exact source-bridge boundary

There is an exact labelled identity of the fourteen normalized Coxeter Gram
matrices, so this is genuinely the same labelled abstract `T1` diagram.  A
literal bridge between the two 24-wall source packages is impossible:

| package | ambient rank/determinant | source norm colours | action orbits |
|---|---|---|---|
| original factory `(14,8,1)` | `14 / -256` | `14*(-2)+10*(-4)` | `8` fixed + `8` pairs |
| compiled `(18,2,0)` | `18 / -4` | `24*(-2)` | `4` fixed + `10` pairs |

Thus neither an ambient integral isometry nor an equivariant 24-wall
bijection exists.  The original factory action has only 20 walls in its 14
active orbits (`8` fixed and `6` pairs); the other four source walls occur in
redundant orbits.  This mismatch is independent of the hidden-root
obstruction.

## What remains true

- `R0`: the compiler's 24 roots fold exactly to its labelled `T1` Gram.
- weak K3 carrier: the primitive K3 extension, generic exact period, and
  genuine Nikulin action are valid.
- `V`: the abstract `T1` polytope remains finite-volume noncompact, with
  `112+15` finite/ideal endpoints, by the pinned factory/Vinberg certificate.
- strong `C`: closed negatively for this compiled carrier.
- `U`: still open for any nontrivial variable class and for a different
  algebraic realization of `T1`.

## Replay

```bash
PYTHONPATH=.:src python3 scripts/certify_h9_14_t1_strong_k3_nef_obstruction.py --check-archive
PYTHONPATH=.:src pytest -q tests/test_h9_14_t1_strong_k3_nef_obstruction.py
```

Artifacts:

- [`certify_h9_14_t1_strong_k3_nef_obstruction.py`](../scripts/certify_h9_14_t1_strong_k3_nef_obstruction.py)
- [`h9_14_t1_strong_k3_nef_obstruction.json`](../data/enumerations/h9_14_t1_strong_k3_nef_obstruction.json)
- [`test_h9_14_t1_strong_k3_nef_obstruction.py`](../tests/test_h9_14_t1_strong_k3_nef_obstruction.py)
- [`probe_h9_14_t1_strong_k3_nef_promotion.py`](../scripts/probe_h9_14_t1_strong_k3_nef_promotion.py)
