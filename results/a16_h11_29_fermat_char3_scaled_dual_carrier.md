# A characteristic-three Fermat carrier for `A16_6 H11/29`

## Positive theorem

The finite-volume noncompact 29-facet `A16_6` Coxeter chamber in `H11` has
an additional algebraic carrier on the supersingular Fermat quartic in
characteristic three.

Let `L` be its archived even rank-12 fixed lattice of determinant `-256`.
The 29 facet reflections preserve a unique maximal integral index-eight
overclosure `M`.  This closure is odd of determinant `-4`.  Its doubled dual

```text
N = M^vee(2),             Gram(N) = 2 Gram(M)^(-1)
```

is even of determinant `-1,024`, with Smith invariants `1^2 2^10`.  The exact
certificate constructs an integral determinant-one matrix `A` satisfying

```text
A^T Gram(NS(X)^C2) A = Gram(N),
NS(X)^C2 = U + E8(-2) + <-2>^2.
```

Here `NS(X)^C2` is the actual saturated fixed-divisor lattice in the
published rank-22 Néron--Severi basis of the characteristic-three
supersingular Fermat quartic, not an abstract genus label.

## The exceptional primitive isotropic bridge

The archived Scharlau--Kirschmer normal form first gives

```text
N ~= U + A1(-1)^10.
```

The obvious copy of `U` has complement `A1(-1)^10`, so it does not display
the Fermat decomposition.  In that normal form, the certificate instead
uses the exceptional primitive isotropic pair

```text
u = (3,4,2,1,1,1,1,1,1,1,1,0),
v = (4,3,2,1,1,1,1,1,1,1,1,0).
```

They satisfy `u^2=v^2=0` and `(u,v)=1`, and the full split basis has
determinant `-1`.  Their positive rank-ten complement has exactly the short
vector behaviour of `E8(2)+<2>^2`.  A frozen determinant-minus-one
PARI/GP `qfisom` matrix is replayed by the permanent certificate through the
entrywise identity

```text
X^T K X = E8(2) + <2>^2.
```

Composing this with the archived standard-to-Fermat basis gives the final
determinant-one isometry.  No classification-only inference is used.

## All 29 oriented halfspaces

Primitive reduction in the odd closure gives the mirror exchange

| primitive mirror in `M` | count | primitive mirror in `N` and `NS(X)^C2` |
|---|---:|---|
| `(-2, div 1)` | 17 | `(-4, div 2)` |
| `(-1, div 1)` | 7 | `(-2, div 2)` |
| `(-4, div 2)` | 5 | `(-2, div 1)` |

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

The multiplier is positive.  All 29 archived presentation orientations are
`+1`; the certificate also checks each archived pairing column.  Hence this
transports every original halfspace, not only the unoriented mirrors.

Every reflection on `N` is the contragredient of its reflection on `M`.
After the explicit Fermat isometry it is integral, unimodular, involutive,
and preserves the actual fixed-divisor Gram.  The primitive 12-by-29 root
matrices have identity Hermite normal form in both bases, so these divisors
generate the whole fixed lattice with index one.  Their oriented sum is a
strict timelike interior vector of square `436` and minimum wall pairing `4`.

## Literal Vinberg finite-volume gate

The inherited literal Vinberg certificate gives:

| invariant | value |
|---|---:|
| dimension | `11` |
| facets | `29` |
| finite vertices | `3,876` |
| cusps | `253` |
| elliptic rank-ten edges | `28,798` |
| connected affine components | `668` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Thus every elliptic edge has exactly two finite-or-ideal endpoints and every
irreducible affine component has the required full parabolic
zero-neighborhood.  The chamber is finite-volume and noncompact.

## What is new, and what is not

The known `H11/29` chamber already has a full-rank complex-K3 divisor carrier
with very-general Néron--Severi lattice
`U+D8(-1)+A1(-1)^2`, of determinant `-16`.  The result here is an
arithmetically different, additional algebraic carrier of determinant
`-1,024` on an actual characteristic-three surface.

The target lattice is itself known: it already carries the Fermat 206-wall
section and Alexeev's 14-wall full reflection chamber.  Consequently this is
not a new lattice or maximal reflection group, and there is no new
abstract-polytope or worldwide-priority claim.

The imported 29-wall arrangement is not the complete nef chamber of the
Fermat quartic.  Seventeen primitive target mirrors have square `-4`, extra
roots may refine the arrangement, and the twelve square-minus-two divisor
classes are not asserted irreducible curves.  No surface automorphism
preserving this transported chamber is claimed; in particular this is not
an automorphism-fixed-section theorem.

## Reproduction

- [certificate generator](../scripts/certify_a16_h11_29_fermat_char3_scaled_dual_carrier.py)
- [exact JSON certificate](../data/enumerations/a16_h11_29_fermat_char3_scaled_dual_carrier.json)
- [focused regression test](../tests/test_a16_h11_29_fermat_char3_scaled_dual_carrier.py)
- [original H11/29 complex-K3 carrier](alexeev_h11_29_ud8_divisor_overlattice_carrier.md)
- [Fermat 206-wall lineage](ss_fermat_char3_c2_h11_lineage.md)

Replay from the repository root:

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_a16_h11_29_fermat_char3_scaled_dual_carrier --check-archive
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_a16_h11_29_fermat_char3_scaled_dual_carrier.py
```
