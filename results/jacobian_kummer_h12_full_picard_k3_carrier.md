# TOOT: the 172-wall Jacobian--Kummer `H12` chamber has a full-Picard K3 carrier

## Positive theorem

The certified finite-volume noncompact chamber `JK-H12-172` admits a
seven-dimensional characteristic-zero family of complex projective K3
surfaces whose very-general full Neron--Severi lattice is

\[
 S=U(2)+D_5(-1)^2+\langle-4\rangle .
\]

In the archive's code notation this is `U(2)+D5(-1)^2+<-4>`.

All 172 oriented primitive wall normals lie in `S`; every associated
reflection acts integrally on all of `S`; and the 172 normals generate `S`
with index one.  Thus this is a full-Picard mixed divisor-reflection carrier,
not merely a rational Lorentzian realization.

The chamber already occurred algebro-geometrically as the trace-nine fixed
section of the complete chamber of a general genus-two Jacobian--Kummer K3
surface `Y16`, where `S` is a saturated fixed divisor sublattice of the
rank-16 Picard lattice.  The new gain is sharper: `S` itself is now the full
Picard lattice of a very-general seven-dimensional K3 family.

## Primitive K3 embedding

The exact carrier invariants are

```text
rank(S)       = 13
signature(S)  = (1,12)
det(S)        = 256
A_S           = (Z/2)^2 + (Z/4)^3
length(A_S)   = 5
q_S           = u(2) + <3/4> + <3/4> + <7/4>.
```

The fixed-lattice certificate reconstructs the saturated integral basis
inside the published rank-16 Jacobian--Kummer Picard lattice.  A unit maximal
minor proves saturation; exact Gram, Smith and finite quadratic-form checks
identify the displayed isometry type.

Nikulin's stable-range primitive-embedding criterion applies directly.  An
orthogonal complement in

\[
 \Lambda_{K3}=U^3+E_8(-1)^2
\]

must have rank 9 and signature `(2,7)`, while `length(A_S)=5`; hence

\[
 9\ge 5+2.
\]

Both complement signs are available, so a primitive embedding
`S -> Lambda_K3` exists.  Period surjectivity and avoidance of the countably
many extra algebraic hyperplanes then give a very-general projective K3 with
`NS(X)=S`.  Its lattice-polarized moduli dimension is `20-13=7`.

The minimal scale for which `S^vee(c)` is even and integral is also recorded:
`c=8`.  The scaled dual has determinant `2^31`, Smith invariants
`2^3,4^2,8^8`, and discriminant length 13.  It is not used for the positive
theorem because the much smaller original fixed lattice embeds already.

## Complete transport of the 172 walls

The published ample class is reconstructed in the saturated fixed basis as

```text
(-6,-6,-6,-6,-6,-6,11,13,11,14,11,13,17).
```

It has square 8 and pairs strictly positively with every primitive wall
normal; the minimum pairing is one.  No post-hoc orientation choice is
needed: all 172 archived orbit sums already have the positive sign.  Their
primitive types are

```text
16 * (-2, divisibility 1)
36 * (-4, divisibility 2)
96 * (-4, divisibility 4)
24 * (-8, divisibility 4).
```

For every normal the certificate checks primitivity, negative square,
strict orientation, integral reflection matrix, involutivity, root negation
and preservation of the full rank-13 Gram matrix.  Smith reduction of the
172-by-13 root matrix is `1^13`: the normals span `S` with index one.

## Finite volume: two independent complete proofs

Finite volume is not inferred from the K3 lattice.  It is replayed by two
independent exact certificates.

The literal Vinberg endpoint/parabolic audit enumerates all 910,316 elliptic
rank-11 edge diagrams.  Every edge has exactly two finite-or-ideal endpoints;
there are 113,274 finite vertices and 6,214 cusps.  Every full affine
zero-neighbourhood is valid, with zero bad or unrecognized edges.  This proves
finite volume, and the cusps prove noncompactness.

Independently, GMP rational double description enumerates every extreme ray
of the full pointed halfspace cone:

| ray kind | count |
|---|---:|
| timelike / finite vertex | 113,274 |
| lightlike / cusp | 6,214 |
| spacelike / hyperideal | 0 |
| **total** | **119,488** |

Every one of the 119,488 incidence sets is distinct and every ray satisfies
all 172 inequalities in exact arithmetic.  The 6,214 lightlike rays agree
exactly with the exhaustive affine neighbourhoods.  Because this is the
complete extreme-ray set, zero hyperideal rays independently closes the
finite-volume gate.

## Exact boundary

This is **not full-nef**.  K3 nef walls come from effective square-minus-two
classes, while 156 of the 172 mirrors have square minus four or minus eight.
No claim is made that those vectors are irreducible curves or that this
Coxeter chamber equals the nef cone of the new very-general K3 surfaces.

The integral Picard reflections are also **not surface automorphisms** and
are not asserted to extend as Hodge isometries.  The abstract chamber and its
original invariant-divisor occurrence were already known in the project; no
new-abstract-polytope or worldwide-priority claim is made.  The positive
advance is the exact full-Picard deformation carrier.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h12_full_picard_k3_carrier.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_full_picard_k3_carrier.py
```

Frozen payload SHA-256:

`07a31cf36fd1e3ca7d548dd89b976d1d6919455243087ec3250c12ea463a1420`
