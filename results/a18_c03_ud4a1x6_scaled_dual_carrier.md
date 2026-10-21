# A scaled-dual full-Picard K3 carrier for `A18-C03 H11/28`

## Positive theorem

The finite-volume noncompact 28-facet `A18-C03` Coxeter chamber in `H11`
has a second full-Picard complex-K3 divisor-lattice realization.

Let `L` be the archived even fixed lattice of rank `12` and determinant
`-256`. Its 28 facet reflections preserve a unique maximal integral
index-four overclosure `M`. This closure is odd of determinant `-16`. Its
doubled dual

```text
N = M^vee(2),             Gram(N) = 2 Gram(M)^(-1)
```

is even of signature `(1,11)`, determinant `-256`, and Smith type
`1^4 2^8`. Its exact 2-elementary invariants are

```text
(r,a,delta) = (12,8,1),
q_N = q_D4(-1) + <-1/2>^6.
```

Thus its standard indefinite representative is

```text
U + D4(-1) + A1(-1)^6.
```

The certificate does more than invoke the classification. It gives an
explicit primitive embedding into the K3 lattice.

## Explicit K3-lattice glue and exact Picard occurrence

Take the signature-`(2,8)` lattice

```text
T = U + U(2) + A1(-1)^6.
```

An eight-bit backtrack exhibits an anti-isometry `q_N -> -q_T` and checks all
`256` discriminant classes, including every quadratic and bilinear value.
The graph glue of `N+T` has index `256`. In the archived basis its Gram
matrix is integral, even, unimodular, and has signature `(3,19)`. It is
therefore the K3 lattice `U^3+E8(-1)^2`. Smith replay of both inclusion
matrices has only unit invariants, so both `N` and `T` are primitive.

K3 period surjectivity now gives an eight-dimensional family of
`N`-polarized projective complex K3 surfaces. Choosing a period outside the
countable union of extra integral Noether--Lefschetz hyperplanes makes

```text
NS(X) = N
```

exactly. This is a full Picard lattice, not merely a divisor sublattice.

## Primitive covectors and all 28 halfspaces

Primitive reduction in the odd closure gives the following exchange:

| primitive mirror in `M` | count | primitive mirror in `N=NS(X)` |
|---|---:|---|
| `(-1, div 1)` | 12 | `(-2, div 2)` |
| `(-2, div 1)` | 8 | `(-4, div 2)` |
| `(-4, div 2)` | 8 | `(-2, div 1)` |

For an oriented primitive closure root `q`, put

```text
d = div_M(q),
p = Gram(M) q / d,
Phi(x) = Gram(M) x.
```

Then

```text
(Phi(x),p)_N = (2/d) (x,q)_M.
```

The multiplier is positive. Seventeen archived primitive presentations have
scale `-1` and eleven have scale `+1`; the certificate applies every recorded
orientation and checks its archived pairing column. Hence all 28 oriented
halfspaces, rather than only the underlying hyperplanes, are transported.

Every target reflection is an integral unimodular involutive isometry of
`N`. It is the transpose of the corresponding reflection on `M`. Pairwise
Gram entries differ only by positive row-and-column factors, so all
normalized Coxeter relations are unchanged.

The 12-by-28 primitive-root matrix has identity Hermite normal form. The
mirrors generate the whole Picard lattice with index one. Their oriented sum
is a strict timelike chamber point of square `608`, with minimum wall pairing
`2`.

## Literal Vinberg finite-volume gate

The inherited literal Vinberg audit gives:

| invariant | value |
|---|---:|
| dimension | `11` |
| facets | `28` |
| finite vertices | `1,252` |
| cusps | `64` |
| elliptic rank-10 edges | `8,851` |
| connected affine components | `186` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Thus every elliptic edge has exactly two finite-or-ideal endpoints and every
irreducible affine component has the required full parabolic
zero-neighborhood. The chamber is finite-volume and noncompact.

## What is new, and what is not

The abstract `H11/28` chamber is not new. It already had a project-level
geometric realization: a symplectic Nikulin involution on a special
Picard-rank-20 K3 surface has fixed Neron--Severi lattice equal to `L`.

The new positive enrichment is a second full-Picard rank-12 K3 carrier:

- the exact scaled-dual lattice `N` occurs as `NS(X)`;
- an explicit discriminant graph glue proves its primitive K3 embedding;
- all 28 oriented mirrors and integral reflections are attached;
- those mirrors generate the complete Picard lattice.

Although `L` and `N` both have determinant `256`, they are arithmetically
different: `L` has Smith type `1^6 2^4 4^2`, while `N` has elementary type
`1^4 2^8`.

The imported arrangement is not the complete nef chamber: eight primitive
mirrors have square `-4`, and further `(-2)` roots may refine it. The twenty
square-minus-two divisor classes are not asserted irreducible. No surface
automorphism preserving the transported chamber is claimed. This
lattice-only certificate originally made no equation claim; a companion
result now supplies an explicit eight-moduli `I0*+6I2+6I1` Weierstrass
family with this exact Picard lattice and a literal 28-wall fiber-basis
marking. There is no new abstract-polytope or worldwide-priority claim. See
the [equation-level carrier](a18_c03_ud4a1x6_weierstrass_carrier.md).

## Reproduction

- [certificate generator](../scripts/certify_a18_c03_ud4a1x6_scaled_dual_carrier.py)
- [exact JSON certificate](../data/enumerations/a18_c03_ud4a1x6_scaled_dual_carrier.json)
- [focused regression test](../tests/test_a18_c03_ud4a1x6_scaled_dual_carrier.py)
- [original special-period K3 realization](a18_c03_c04_special_k3_symplectic_lifts.md)
- [explicit eight-moduli Weierstrass carrier](a18_c03_ud4a1x6_weierstrass_carrier.md)

Replay from the repository root:

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_a18_c03_ud4a1x6_scaled_dual_carrier --check-archive
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_a18_c03_ud4a1x6_scaled_dual_carrier.py
```
