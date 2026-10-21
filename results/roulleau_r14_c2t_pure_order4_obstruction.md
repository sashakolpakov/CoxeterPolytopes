# R14 transposition: exact-fixed pure order-four K3 lift obstruction

## Outcome

The Roulleau `U+E8(-1)+A1(-1)^4` transposition defining the `H12/23`
section has **no purely non-symplectic order-four K3 automorphism lift whose
global fixed lattice has rank 13**, at any period.  Equivalently, it has no
such lift realizing `H12/23` as the full automorphism-fixed hyperbolic
section.

This retires the proposed pure order-four special-period route.  In
particular, the earlier candidate with square-invariant Picard lattice
`U+E8(-1)+D4(-1)+A1(-1)^2` was not a positive example: its algebraic
coinvariant contains a `(-2)`-vector, so it fails the invariant Kähler-class
gate.  That invalid candidate script has been removed.

The hyperbolic result itself is unchanged.  The archived section is a genuine
finite-volume, noncompact Coxeter polytope in `H^12`, with 23 facets, 320
finite vertices, 17 cusps, 2,308 elliptic edges, and no bad edge or bad
parabolic neighborhood.

## Determinant/classification gate

Let `F` be the fixed lattice of the archived divisor-lattice transposition.
Exact reconstruction gives

```text
rank(F) = 13
det(F)  = 16
A_F     = Z/2 + Z/2 + Z/4.
```

For a purely non-symplectic order-four K3 action, write `(r,l,m)` for the
multiplicities of `1`, `-1`, and the conjugate pair `i,-i` on `H^2`.  The
order-four Lefschetz relations

```text
r + l + 2m = 22,
4 alpha = r - l - 2
```

leave only

| `(r,l,m)` | `alpha` |
|---|---:|
| `(13,3,3)` | 2 |
| `(13,7,1)` | 1 |

The complete Brandhorst--Hofmann pure order-four classification contains
exactly one action with each character:

| Character | Database ID | Fixed determinant | Fixed Smith factors |
|---|---|---:|---|
| `(13,3,3)` | `0.4.2.9` | 64 | `1^10,4,4,4` |
| `(13,7,1)` | `0.4.0.2` | 256 | `1^6,2^6,4` |

Any exact-fixed lift fixes every vector of `F`.  Since `F` and the global
fixed lattice `S(g)` both have rank 13, this would be a finite-index containment
`F <= S(g)`, hence

```text
det(F) = [S(g):F]^2 det(S(g)).
```

The two required squared indices would be `16/64=1/4` and `16/256=1/16`.
Neither is an integer square.  Both classified actions are therefore
impossible extensions of the archived transposition.

## Reproducibility and provenance

The certifier checks the two pinned fixed Gram matrices offline.  It can also
replay all 79 rows of the external classification directly from
`K3Groups.zip`:

```bash
PYTHONPATH=.:src python3 scripts/certify_roulleau_r14_c2t_pure_order4_obstruction.py \
  --brandhorst-hofmann-zip /path/to/K3Groups.zip --check-archive
```

External source: Brandhorst--Hofmann, *Finite subgroups of automorphisms of
K3 surfaces*, [classification dataset](https://doi.org/10.5281/zenodo.10410278),
version 2.1, CC-BY-4.0.  The pinned ZIP MD5 is
`143047b9cb76f4ad4cb2aadd5959d879`; its SHA-256 is
`c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c`.
The internal `results/purely_ns/order04.txt` SHA-256 is
`d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8`.

## Scope

This is a complete obstruction only for a **purely non-symplectic order-four
lift with global fixed rank 13**.  It does not rule out a pure order-four lift
whose global fixed lattice strictly enlarges the archived one (which would no
longer make `H12/23` the full fixed section), an order-four action with
symplectic square, or a higher-order lift.  It does not invalidate the
abstract `H12/23` Coxeter polytope.  Those possibilities belong to different
search lanes.

Artifacts:

- `scripts/certify_roulleau_r14_c2t_pure_order4_obstruction.py`
- `data/enumerations/roulleau_r14_c2t_pure_order4_obstruction.json`
- `tests/test_roulleau_r14_c2t_pure_order4_obstruction.py`
