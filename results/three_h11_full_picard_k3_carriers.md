# TOOT: three `H11` chambers gain full-Picard complex-K3 carriers

## Positive theorem

The already certified finite-volume noncompact chambers

- `F3-C2-H11-206`, with 206 facets;
- `F3-V4-H11-328`, with 328 facets; and
- `SS2-3B-H11-22`, with 22 facets

all admit characteristic-zero algebraic carriers.  More precisely, each is
realized in the positive cone of the full Neron--Severi lattice of an
eight-dimensional very-general family of complex projective K3 surfaces.
Every recorded primitive wall normal belongs to that full Picard lattice,
and its Coxeter reflection is integral on the full lattice.  In all three
cases the complete set of wall normals generates the full Picard lattice,
not a proper finite-index sublattice.

For the 206- and 22-wall cases the original wall lattice is already the K3
carrier.  The 328-wall case has one essential correction: 36 reflections
are nonintegral on the original characteristic-three divisor-fixed lattice,
so the carrier is its certified minimal even reflection-stable index-four
overlattice.

This supplies three new characteristic-zero, full-Picard algebraic
realizations of known abstract Coxeter polytopes.  It is not a claim that
the abstract polytopes themselves are new.

## Exact carrier census

| chamber | full Picard carrier `S` | `det(S)` and `A_S` | primitive K3-embedding proof | wall span |
|---|---|---|---|---:|
| `F3-C2-H11-206` | `U+E8(-2)+<-2>^2` | `2^10`; `(Z/2)^10`, length 10 | explicit embedding | index 1 |
| `F3-V4-H11-328` | minimal stable closure, the negative of Scharlau--Kirschmer List11/Gram11 entry 20 | 48; `(Z/2)^3 + Z/6`, length 4 | Nikulin stable range `10 >= 4+2` | index 1 |
| `SS2-3B-H11-22` | `U+A2(-1)^5` | `3^5`; `(Z/3)^5`, length 5 | Nikulin stable range `10 >= 5+2` | index 1 |

Every carrier has rank 12 and signature `(1,11)`.  Its primitive embedding
in

\[
 \Lambda_{K3}=U^3+E_8(-1)^2
\]

therefore has complement of signature `(2,8)`.  Period surjectivity, followed
by avoidance of the countably many periods with an extra algebraic class,
gives a very-general projective K3 with `NS(X)=S`, rather than merely a K3
whose Picard lattice contains `S`.  The dimension is `20-rank(S)=8`.

For `F3-C2-H11-206`, the discriminant length equals the complement rank, so
the stable-range shortcut is unavailable.  The certificate instead gives a
22-by-12 primitive embedding matrix:

- `U` maps to an ambient `U` summand;
- `E8(-2)` maps diagonally to `E8(-1)+E8(-1)`; and
- the two `<-2>` generators map to `e-f` in the remaining two `U` summands.

The orthogonal complement is `E8(-2)+<2>^2`.  The verifier checks the two
Gram matrices, orthogonality, primitivity, ambient K3 signature and
unimodularity, and the combined gluing index 1,024.

For the other two carriers, Nikulin's sufficient stable-range criterion is
literal: the complement rank is 10, whereas the discriminant lengths are 4
and 5.  Both inequalities have a margin of at least three, and both signs
remain available in the complement.

## Why the index-four closure matters

The original saturated `F3-V4` fixed lattice has determinant 768 and Smith
invariants

```text
1^6, 2^4, 4, 12.
```

Exactly 36 of its 328 wall reflections are nonintegral.  This is a genuine
obstruction to using that original lattice as the divisor carrier.  The
already certified closure computation produces the unique minimal even
overlattice on which all 328 reflections are integral.  It has index four,
determinant `768/4^2=48`, and Smith invariants

```text
1^8, 2^3, 6.
```

The positive statement uses this closure as the full Picard lattice; it does
not silently treat the nonintegral original action as an integral one.

## Complete wall and volume audits

| chamber | primitive mirror types `(square, divisibility)` | strict interior | independent finite-volume certificate |
|---|---|---|---|
| `F3-C2-H11-206` | `12*(-2,1) + 164*(-2,2) + 30*(-4,2)` | square 4; minimum pairing 1 | 327,566 finite vertices, 18,354 cusps, 2,688,957 elliptic rank-10 edges; zero bad edges and zero invalid full affine neighbourhoods |
| `F3-V4-H11-328` | `28*(-2,1) + 24*(-4,2) + 276*(-6,3)` | square 4; minimum pairing 1 | exact GMP cone: 225,408 finite vertices, 11,554 cusps, zero hyperideal rays |
| `SS2-3B-H11-22` | `12*(-2,1) + 10*(-6,3)` | square 14; minimum pairing 1 | 2,304 finite vertices, 72 cusps, 14,580 elliptic rank-10 edges; zero bad edges and zero invalid full affine neighbourhoods |

The archive contains all 556 oriented wall vectors, not only their type
histograms.  Every vector is checked for primitivity, negative square,
strict orientation by the displayed positive vector, integral reflection,
involutivity, root negation and preservation of the carrier Gram matrix.
Smith reduction of each complete root matrix gives rank 12 and index one.

For the two Fermat rows, the displayed positive vectors are transported all
the way back through the certified fixed-lattice bases and equal the
published Fermat ample class.  For the characteristic-two row, the square-14
interior is derived without a hardcoded choice: it is one third of the sum
of all 22 archived primitive roots, pairing by 1 with the twelve
square-minus-two walls and by 3 with the ten square-minus-six walls.  The
identification `U+A2(-1)^5` is independently replayed from the pinned
Scharlau--Kirschmer entry-27 lineage certificate.

Finite volume is kept logically separate from the K3-embedding argument.
For the 206- and 22-wall systems it comes from the literal Vinberg
edge-endpoint test together with every full affine zero-neighbourhood.  For
the 328-wall system it comes from complete exact double description, with no
hyperideal ray.  Thus no finite-volume conclusion is inferred merely from a
root group, a list of cusps, or a partial Coxeter diagram.

## Scaled-dual census

The minimal positive scale `c` for which `S^vee(c)` is even and integral is
also frozen, although none of these scaled duals is needed for the positive
K3 realization:

| carrier | minimal `c` | determinant of `S^vee(c)` | discriminant length |
|---|---:|---:|---:|
| `U+E8(-2)+<-2>^2` | 4 | 16,384 | 12 |
| 328-wall stable closure | 6 | 45,349,632 | 11 |
| `U+A2(-1)^5` | 3 | 2,187 | 7 |

For comparison, the rejected original 328-wall fixed lattice requires
`c=12`; its scaled dual has determinant 11,609,505,792 and discriminant
length 11.

## Exact boundary

These are mixed divisor-reflection carriers, **not full-nef** realizations.
A K3 nef chamber is cut out by effective square-minus-two classes, whereas
the three systems also contain square-minus-four or square-minus-six walls.
No claim is made that those higher-norm walls are irreducible curves or that
the Coxeter chamber equals the K3 nef cone.  The integral Picard reflections
are also **not surface automorphisms** and are not asserted to extend as
Hodge isometries.  No worldwide-priority claim is made.

What is new here is the exact carrier theorem: three certified cofinite
`H11` chambers, including 556 walls in total, now have characteristic-zero
very-general K3 families with the complete chamber-wall system living on the
full divisor lattice.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_three_h11_full_picard_k3_carriers.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_three_h11_full_picard_k3_carriers.py
```

Frozen payload SHA-256:

`58cf2623ba476ec5ee0d19a173864fa5955a9e7a49f8e3948b0a4d87831252a2`
