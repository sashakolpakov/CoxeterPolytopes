# Complex-K3 carrier for the `SS2-2A H13/50` chamber

## Positive theorem

The certified finite-volume noncompact 50-facet Coxeter chamber in `H13`
obtained from the geometric `2a` class of the characteristic-two
Borcherds--Dolgachev--Kondō stock has a second, genuinely complex-K3 divisor
realization.

Let `M` be its archived fixed lattice.  It is the even 2-elementary lattice
`(14,10,0)`, of determinant `-1024`.  Pass to the scaled dual

```text
N = M^vee(2),             Gram(N) = 2 Gram(M)^(-1).
```

Then `N` is even of determinant `-16`, with Smith invariants
`1^10 2^4`.  The certificate gives an explicit integral unimodular matrix
`X`, of determinant `-1`, satisfying

```text
X^T Gram(U+D8(-1)+D4(-1)) X = 2 Gram(M)^(-1).
```

Thus

```text
N = U + D8(-1) + D4(-1),
```

the actual Néron--Severi lattice of the archived Brandhorst--Hofmann/Roulleau
complex K3 carriers.  Two independently replayed Brandhorst--Hofmann rows
already certify genuine complex K3 surfaces with this exact Picard lattice.

This is a scaled-dual, contragredient construction, not a same-form
overlattice of `M`.

## All 50 mirrors survive

For a primitive source normal `r`, put

```text
p = Gram(M) r / gcd(Gram(M) r).
```

Writing `d=div_M(r)` and `Phi(x)=Gram(M)x`, the certificate also checks the
halfspace identity

```text
(Phi(x),p)_N = (2/d) (x,r)_M.
```

The factor `2/d` is positive.  Therefore all inequality signs and the entire
50-wall chamber are transported, not merely the abstract reflection
involutions.

These are its primitive coordinates in `N=M^vee(2)`.  The mirror exchange is
exact:

```text
40 source (-4, divisibility 2) mirrors -> 40 primitive (-2), divisibility 1;
10 source (-2, divisibility 1) mirrors -> 10 primitive (-4), divisibility 2.
```

The 14-by-50 matrix of the `p` has Hermite normal form equal to the identity.
Consequently the 50 primitive mirrors generate the whole K3 Picard lattice,
not merely a proper finite-index sublattice.

For every facet the certificate checks the source type, primitive reduction,
the scaled-dual square law, and the contragredient identity

```text
R_N = R_M^T.
```

After applying `X`, every target reflection is integral, unimodular,
involutive, Gram-preserving, and satisfies `R_K3 X = X R_N`.  Hence all 50
hyperplanes are honest algebraic divisor mirrors on the complex K3 carrier.

The sum of the 50 primitive normals is a strict timelike carrier interior of
square `2012`; its minimum wall pairing is `6`.  A generic nearby class avoids
all extra root hyperplanes.  Weyl-transporting the marking to the ample cone
and applying K3 Riemann--Roch chooses the effective signs of the 40 primitive
`(-2)` classes simultaneously.  These effective divisors may be reducible
and are not asserted to be nef-cone facets.  The 10 primitive `(-4)` divisibility-two
mirrors remain algebraic divisor reflections, not asserted curve walls.

## Two independent finite-volume audits

The inherited literal Vinberg audit gives:

| invariant | value |
|---|---:|
| dimension | `13` |
| facets | `50` |
| finite vertices | `45,954` |
| cusps | `2,598` |
| elliptic rank-12 edges | `391,614` |
| connected affine components | `6,210` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Thus every elliptic edge has exactly two finite-or-ideal endpoints and every
affine component has the required full zero-neighborhood.  The chamber is
finite-volume and noncompact.

An independent GMP-cdd exact double-description audit starts from the pinned
50-by-14 inequality matrix and obtains exactly `48,552` extreme rays:

```text
45,954 timelike + 2,598 isotropic + 0 hyperideal.
```

The cusp active-facet incidence split is

```text
1,656 on 14 facets, 872 on 15, 68 on 16, and 2 on 17.
```

This exactly reproduces the Vinberg finite/ideal counts by a separate cone
calculation.

## Exact boundary

This is not the full nef cone: the 10 primitive `(-4)` mirrors are not K3
curve walls, and additional `(-2)` roots of the Picard lattice may refine the
50-wall chamber.  It is also not asserted that either archived
Brandhorst--Hofmann automorphism preserves this transported chamber.

This is not a new abstract polytope.  The literal `H13/50` chamber and its
characteristic-two geometric construction were already certified.  The new
positive gain is a complete reflection-by-reflection realization on a
different, characteristic-zero K3 Picard lattice.  No worldwide priority
claim is made.

## Reproduction

- [certificate generator](../scripts/certify_ss2_2a_h13_50_ud8d4_divisor_carrier.py)
- [exact JSON certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_divisor_carrier.json)
- [focused regression test](../tests/test_ss2_2a_h13_50_ud8d4_divisor_carrier.py)
- [original `H13/50` lineage theorem](borcherds_char2_h13_h14_lineage.md)
- [complex `U+D8+D4` K3 theorem](roulleau_ud8d4_bh_k3_lifts.md)

Replay from the repository root:

```bash
python3 -m scripts.certify_ss2_2a_h13_50_ud8d4_divisor_carrier --check-archive
python3 -m pytest -q tests/test_ss2_2a_h13_50_ud8d4_divisor_carrier.py
```
