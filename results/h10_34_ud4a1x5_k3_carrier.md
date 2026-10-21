# A common full-Picard K3 carrier for the known `H10/34` chamber

Last updated: 2026-08-07.

## Positive result

The finite-volume noncompact `H10/34` Coxeter polytope previously obtained in
two different divisor-lattice ways,

- `A15-C2-H10-34`, inside Alexeev's `(15,9,1)` Hilbert-square chamber; and
- `SS2-V4-H10-34`, inside the characteristic-two Borcherds chamber,

has a common realization on the full Picard lattice of a complex projective
K3 surface:

`NS(X) = U + D4(-1) + A1(-1)^5`.

This is Roulleau's actual rank-11 K3 stock with its published complete
90-curve nef chamber, not merely an atlas-header match or an abstract lattice
allowed by signature.

The 34 imported primitive mirrors have types

`8*(-2,div 1) + 16*(-2,div 2) + 10*(-4,div 2)`.

They generate the whole Picard lattice (Hermite index one), and every one of
their 34 reflection matrices is integral, unimodular, involutive, and an exact
Picard-lattice isometry.

## Why both old integral realizations land on this K3

The reflection-stable discriminant denominators are now exhaustive.

For the A15 fixed lattice, of determinant `2048`, the only stable integral
overlattice indices are `1` and `4`.  The index-four closure is even, has
determinant `128` and Smith invariants `1^4 2^7`, and an explicit determinant
`-1` isometry identifies it with Roulleau's displayed
`U+D4(-1)+A1(-1)^5` basis.

For the SS2 `V4` fixed lattice, of determinant `1024`, the stable integral
indices are exactly `1,2,4,8`; the stable even indices are exactly `1,4`.
The unique maximal index-eight closure `M` is odd of determinant `16`, but its
doubled dual

`N = M^vee(2),   Gram(N)=2 Gram(M)^(-1)`

is even of determinant `128` and Smith invariants `1^4 2^7`.  An explicit
determinant-one isometry `N -> M_A15` sends all 34 primitive oriented SS2
mirrors exactly to the A15 mirrors under the already certified 34-facet
permutation.  Thus the agreement is integral and wall-by-wall, not just a
normalized-Gram coincidence.

The common discriminant quadratic-form histogram is

`q=0:36, q=1:28, q=1/2:28, q=3/2:36`.

An independent order-128 discriminant anti-isometry with

`T = U + U(2) + D4(-1) + A1(-1)^3`

glues `NS(X)+T` to an even unimodular lattice of signature `(3,19)`.  The two
embeddings are primitive.  This supplies a direct K3-lattice embedding proof
in addition to the named Roulleau realization.

## An explicit ample class inside the 34-wall chamber

The sum of the 34 transported roots is timelike of square `1648` and has
minimum imported-wall pairing `4`.  Deterministic Weyl reduction through the
published Roulleau curve walls uses only five reflections, with one-based wall
sequence

`20, 4, 14, 16, 3`.

The reduced sum is big and nef: its minimum pairing with the complete 90-wall
nef boundary is zero.  Adding Roulleau's published ample class

`D18=(1,5,2,4,0,0,0,0,0,0,0),  D18^2=18`

gives the explicit ample divisor

`A=(1,125,54,98,-26,-13,-13,-26,-15,0,-11),  A^2=2178`.

Its minimum pairing with all 90 published nef walls is `1`, while its minimum
pairing with the 34 imported mirrors is `4`.  Hence the imported chamber has a
certified nonempty ample interior on an actual K3 surface.

## Literal Vinberg status

Both upstream presentations replay the same literal finite-volume audit:

| datum | count |
|---|---:|
| finite vertices | 1,774 |
| ideal vertices | 146 |
| elliptic rank-9 edge diagrams | 12,792 |
| finite--finite / finite--ideal / ideal--ideal edges | 6,046 / 5,648 / 1,098 |
| connected affine components | 426 |
| bad edges / bad parabolic neighbourhoods | `0 / 0` |

Thus the chamber is finite-volume and noncompact by Vinberg's edge-endpoint
and full affine zero-neighbourhood criterion.

## Scope boundary

This is a positive algebraic realization of an already known abstract
polytope, not a new `H10/34` combinatorial type.  The imported 34-wall chamber
is not asserted to equal the complete 90-wall K3 nef chamber.  The ten
square-minus-four mirrors are not curve walls.  The 24 ample-positive
square-minus-two classes are effective by K3 Riemann--Roch, but are not all
asserted irreducible.  Lattice reflections are not being promoted to surface
automorphisms.

Replayable artifacts:

- [`h10_34_ud4a1x5_k3_carrier.json`](../data/enumerations/h10_34_ud4a1x5_k3_carrier.json)
- [`certify_h10_34_ud4a1x5_k3_carrier.py`](../scripts/certify_h10_34_ud4a1x5_k3_carrier.py)
- [`test_h10_34_ud4a1x5_k3_carrier.py`](../tests/test_h10_34_ud4a1x5_k3_carrier.py)

Canonical payload SHA-256:

`55b990d3002e787f7199ca6577e2e39a05ccf713ca05365e08fafd8212395e88`.
