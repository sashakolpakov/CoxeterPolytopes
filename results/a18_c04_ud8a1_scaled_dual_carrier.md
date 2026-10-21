# A scaled-dual `U+D8+A1` K3 carrier for `A18-C04 H10/28`

## Positive theorem

The finite-volume noncompact 28-facet `A18-C04` Coxeter chamber in `H10`
has a second, arithmetically different K3 divisor-lattice realization.

Let `L` be its archived even fixed lattice of rank `11` and determinant
`1,024`.  Its 28 facet reflections preserve a unique maximal integral
index-two overclosure `M`.  This closure is odd of determinant `256`.  Its
doubled dual

```text
N = M^vee(2),             Gram(N) = 2 Gram(M)^(-1)
```

is even of determinant `8`, with Smith invariants `1^8 2^3`.  The exact
certificate freezes an integral determinant-one matrix `X` satisfying

```text
X^T Gram(U+D8(-1)+A1(-1)) X = Gram(N).
```

The target is not an abstract genus label.  It is Roulleau's archived actual
complex-K3 Picard stock

```text
NS(Y) = U + D8(-1) + A1(-1),
```

in the published curve-class basis `A1,...,A11`.  Its own complete nef
chamber has 14 curve walls.  The theorem imports the different 28-wall mixed
reflection chamber into the full divisor lattice of such a K3 surface.

## Primitive covectors and all 28 halfspaces

Primitive reduction in the odd closure is essential.  The mirror exchange is

| primitive mirror in `M` | count | primitive mirror in `N` and `NS(Y)` |
|---|---:|---|
| `(-4, div 2)` | 16 | `(-2, div 1)` |
| `(-1, div 1)` | 8 | `(-2, div 2)` |
| `(-2, div 1)` | 4 | `(-4, div 2)` |

For an oriented primitive closure root `q`, put

```text
d = div_M(q),
p = Gram(M) q / d,
Phi(x) = Gram(M) x.
```

Then every `p` is primitive and

```text
(Phi(x),p)_N = (2/d) (x,q)_M.
```

The multiplier is positive.  Nineteen archived primitive presentations have
scale `-1` and nine have scale `+1`; the certificate applies those recorded
orientations before primitive reduction and checks the archived pairing
column exactly.  Hence this is transport of all 28 original halfspaces, not
only of their unoriented hyperplanes.

For every facet, the reflection on `N` is the transpose of the reflection on
`M`.  Under `X`, it becomes an integral, unimodular, involutive isometry of
the actual K3 Picard lattice.  Pairwise Gram entries differ only by the
positive row-and-column covector factors, so every normalized Coxeter
relation is unchanged.

The Hermite normal forms of both 11-by-28 primitive-root matrices are the
identity.  The transported divisors therefore generate the whole target
Picard lattice with index one.  Their oriented sum is a strict timelike
interior class of square `508`, with minimum wall pairing `3`.

## Literal Vinberg finite-volume gate

The inherited literal Vinberg audit gives:

| invariant | value |
|---|---:|
| dimension | `10` |
| facets | `28` |
| finite vertices | `1,087` |
| cusps | `74` |
| elliptic rank-9 edges | `6,956` |
| connected affine components | `204` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Thus every elliptic edge has exactly two finite-or-ideal endpoints and every
irreducible affine component has the required full parabolic
zero-neighborhood.  The chamber is finite-volume and noncompact.

## What is new, and what is not

The abstract `H10/28` chamber is not new.  It already had a project-level
geometric realization: a symplectic Nikulin involution on a special
Picard-rank-19 K3 surface has fixed Néron--Severi lattice equal to the
original determinant-`1,024` lattice `L`.

The new positive enrichment is a genuinely different full-rank carrier:

- target Picard determinant `8`, rather than `1,024`;
- an actual very-general Roulleau K3 divisor stock;
- an explicit stable-closure/scaled-dual isometry;
- exact transport of all oriented facets and reflections;
- index-one generation of the complete target Picard lattice.

This 28-wall mixed arrangement is not the K3's complete nef chamber.  Four
primitive target mirrors have square `-4`, and additional `(-2)` roots may
refine the imported chamber.  The 24 square-minus-two divisor classes are not
asserted irreducible.  No surface automorphism preserving the transported
chamber is claimed.  There is no new abstract-polytope or worldwide-priority
claim.

## Reproduction

- [certificate generator](../scripts/certify_a18_c04_ud8a1_scaled_dual_carrier.py)
- [exact JSON certificate](../data/enumerations/a18_c04_ud8a1_scaled_dual_carrier.json)
- [focused regression test](../tests/test_a18_c04_ud8a1_scaled_dual_carrier.py)
- [original special-period K3 realization](a18_c03_c04_special_k3_symplectic_lifts.md)
- [Roulleau rank-11 stock theorem](roulleau_rank11_direct_graph_denominator.md)

Replay from the repository root:

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_a18_c04_ud8a1_scaled_dual_carrier --check-archive
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_a18_c04_ud8a1_scaled_dual_carrier.py
```

