# Hilbert-square scaled-dual carrier for `SS2-2C H14/98`

## Positive theorem

The certified finite-volume noncompact `SS2-2C` 98-wall chamber in
`H14` has a complete characteristic-zero hyperkähler divisor-lattice
realization.

There is one necessary intermediate step.  Let `L` be the archived even
rank-15 fixed lattice of determinant `256`.  Its 98 reflections preserve a
unique integral index-two overclosure `M`.  The lattice `M` is odd, has
determinant `64`, and is not itself the target divisor lattice.  Its doubled
dual

```text
N = M^vee(2),             Gram(N) = 2 Gram(M)^(-1)
```

is even of determinant `512`, with Smith invariants `1^6 2^9`.  The exact
certificate gives an integral determinant-one matrix `X` satisfying

```text
X^T Gram(Pic(Hilb^2(X))) X = Gram(N).
```

Here `X` is a K3 surface with `NS(X)=(14,8,1)`, and the target is Alexeev's
archived actual divisor lattice

```text
Pic(Hilb^2(X)) = (15,9,1).
```

Thus the characteristic-two fixed-lattice construction and the classical
characteristic-zero Hilbert-square construction are joined at the level of
the full integral lattice and every oriented wall, not merely by an equal
Coxeter edge histogram.

## Why the odd closure is essential

Taking `2 L^(-1)` directly gives the wrong determinant.  In the stable
closure `M`, primitive reduction changes the 98 source mirrors as follows:

| primitive type in `M` | count | primitive type in `N= M^vee(2)` |
|---|---:|---|
| `(-1, div 1)` | 70 | `(-2, div 2)` |
| `(-4, div 2)` | 14 | `(-2, div 1)` |
| `(-2, div 1)` | 14 | `(-4, div 2)` |

The 70 odd roots arise because 70 primitive roots of `L` become twice a
primitive vector in `M`.  Their hyperplanes and reflections do not change.
This accounts exactly for the short/long exchange that was visible in the
earlier 98-by-98 Gram comparison.

## Exact positive halfspace transport

For a primitive closure root `q`, set

```text
d = div_M(q),
p = Gram(M) q / d,
Phi(x) = Gram(M) x.
```

Then `p` is primitive in `N`, and the certificate checks for every facet

```text
(Phi(x),p)_N = (2/d) (x,q)_M.
```

The factor is positive.  The original `L`-root is itself a positive multiple
of `q`, so all 98 original inequality signs survive as well.  Under the
explicit isometry `X`, each `p` equals the corresponding wall in Alexeev's
archived 98-wall list under the frozen facet permutation.

The primitive target split is

```text
14 * (-2,div1) + 70 * (-2,div2) + 14 * (-4,div2).
```

The Hermite normal form of the 15-by-98 target-root matrix is the identity.
Consequently these 98 divisors generate all of `Pic(Hilb^2(X))`, with index
one.  Every target reflection is integral, determinant minus one,
involutive, and preserves the Beauville--Bogomolov Gram matrix.  Facet by
facet, the reflection on the doubled dual is the transpose of the closure
reflection and is intertwined by `X` with the Hilbert-square reflection.

The sum of the 98 primitive doubled-dual normals is a strict timelike
interior class of square `40,432`; its minimum wall pairing is `76`.

## Literal Vinberg finite-volume gate

The transported chamber inherits the already independent literal Vinberg
audit, because the positive linear similarity transports the entire cone:

| invariant | value |
|---|---:|
| dimension | `14` |
| facets | `98` |
| finite vertices | `75,898` |
| cusps | `2,114` |
| elliptic rank-13 edges | `646,786` |
| connected affine components | `5,838` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
irreducible affine component has the required full parabolic
zero-neighborhood.  Hence the chamber is finite-volume and noncompact; this
is not inferred merely from the existence of cusps.

## Exact novelty boundary

This is a positive carrier bridge, but it does **not** add a new abstract
Coxeter polytope or a new geometric-chamber slot.  Alexeev already published
the `(15,9,1)` 98-wall reflection chamber on the Hilbert-square Picard
lattice, and the earlier lineage theorem proved that `SS2-2C H14/98` is the
same abstract chamber after doubled duality.

The new project-level gain is the stronger coordinate statement:

- the unique stable closure is included explicitly;
- the doubled dual is explicitly isometric to the actual hyperkähler Picard
  lattice;
- all 98 oriented halfspaces, primitive roots, and reflections are carried
  across;
- the transported roots generate the full target divisor lattice.

No claim is made that all 98 roots are MBM or prime-exceptional divisors, or
that this chamber is the Kähler or movable cone.  No automorphism of
`Hilb^2(X)` preserving the selected chamber is asserted.  The characteristic-
two `2c` Frobenius class remains outside the geometric K3 subgroup.  There is
no worldwide priority claim.

## Reproduction

- [certificate generator](../scripts/certify_ss2_2c_h14_98_hilb2_scaled_dual_carrier.py)
- [exact JSON certificate](../data/enumerations/ss2_2c_h14_98_hilb2_scaled_dual_carrier.json)
- [focused regression test](../tests/test_ss2_2c_h14_98_hilb2_scaled_dual_carrier.py)
- [original `H13/H14` lineage theorem](borcherds_char2_h13_h14_lineage.md)
- [algebraic reflection baseline](higher_dimensional_algebraic_reflection_examples.md)

Replay from the repository root:

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_ss2_2c_h14_98_hilb2_scaled_dual_carrier --check-archive
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_ss2_2c_h14_98_hilb2_scaled_dual_carrier.py
```

