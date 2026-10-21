# Polytope-first OG10 carrier theorem for `BOR-D4-H21-210`

## Result and directional status

Starting only from the frozen `210 x 210` Coxeter Gram matrix of the
classical Borcherds `D4` chamber, one can reconstruct a unique primitive
embedding of its rank-22 wall lattice

```text
L = U + D20(-1)
```

in the OG10 Beauville lattice.  The corresponding zero-dimensional marked
OG10 period locus with exact Picard lattice `L` is nonempty, and all 210 input
facets become genuine OG10 wall-divisor classes.

The stronger cone promotion is impossible.  In the reconstructed basis,

```text
v = (-1,1,0,...,0)
```

is an additional primitive OG10 wall divisor of square `-4` and ambient
divisibility one.  The positive vector

```text
x = G^(-1) (1,...,1)^t,       x^2 = 581
```

pairs strictly positively with all 210 input walls—minimum `1`, maximum
`421`—but satisfies `(x,v)=0`.  Thus `v^perp` cuts the strict interior of the
210-wall Coxeter chamber.  No OG10 Kähler/nef chamber can equal that chamber.

This is a genuine

```text
POLYTOPE -> ALGEBRAIC VARIETY -> CONE CRITERION
```

calculation.  It is **not** a nonexistence theorem for the hyperbolic
polytope, which exists classically, and it does not claim that OG10 carriers
exhaust all possible geometric carrier species.  Under the strengthened
project mandate it is calibration for the future uniform transfer theorem,
not completion of that mandate.

## 1. Raw Coxeter reconstruction

The input is Alexeev's pinned `(22,2,0)` matrix, SHA-256

```text
4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457
```

The first 22 facets form a nonsingular block `G` with

```text
det G       = -4
signature   = (1,21)
Smith form  = 1^20,2,2.
```

If `M` is the full Coxeter Gram matrix, the matrix

```text
C = G^(-1) M[1..22,*]
```

is integral and satisfies `C^t G C=M` exactly.  Hence no prior identification
with a K3 surface, OG10 manifold, or named lattice is needed to recover the
rank-22 wall lattice and all 210 wall coordinates.  Their intrinsic types in
`L` are

```text
42  x (square -2, divisibility 1),
168 x (square -4, divisibility 2).
```

The identification `L = U+D20(-1)` agrees with the classical Borcherds and
Dolgachev--Kondō lineage; it is not a new polytope.

## 2. Exhaustive primitive embedding in the OG10 lattice

The discriminant group of `L` is `(Z/2)^2`.  Its three nonzero classes all
have discriminant value `-1 mod 2Z`.  An OG10 lattice has determinant three,
so a primitive embedding must kill the entire 2-primary discriminant.  The
gluing subgroup therefore has order four, and the positive rank-two
orthogonal complement `K` must satisfy

```text
det K = 3 * 4^2 / 4 = 12.
```

Gauss reduction gives exactly two positive even binary forms of determinant
12:

```text
<2> + <6>                 q on 2-torsion = {1/2,3/2,0},
A2(2) = [[4,2],[2,4]]     q on 2-torsion = {1,1,1}.
```

Only `A2(2)` can anti-glue to `L`.  There are 36 ordered choices of generators
and six distinct graph subgroups.  The integral isometry group of `A2(2)` has
order 12 and induces the full `S3` of order six on its nonzero 2-torsion.
Consequently all six anti-isometries lie in one orbit: the primitive embedding
orbit is unique.

Using either sign convention for the off-diagonal entry, two displayed glue
vectors give an even integral overlattice with

```text
signature       (3,21)
determinant     -3
Smith form      1^23,3
discriminant q  4/3 mod 2Z.
```

This is

```text
Lambda_OG10 = 3U + 2E8(-1) + A2(-1).
```

The embedded copy of `L` is primitive.  The 42 short walls retain type
`(-2,div 1)`, while all 168 long walls change from `(-4,div 2)` in `L` to
`(-4,div 1)` in `Lambda_OG10`.

[Rapagnetta](https://arxiv.org/abs/math/0606409) determines the OG10 second
Betti number and Beauville lattice.  Period-map surjectivity for marked
hyperkähler manifolds, in the form proved by
[Huybrechts](https://arxiv.org/abs/alg-geom/9705025), supplies a marked OG10
manifold with this Hodge lattice.  Since the orthogonal complement is positive
rank two, the exact-Picard period locus is zero-dimensional; since `L`
contains positive classes, it is projective.

## 3. Necessary-and-sufficient OG10 cone gate

[Mongardi--Onorati](https://arxiv.org/abs/2010.12511), in their wall-divisor
classification and Kähler-cone theorem, prove that primitive OG10 wall divisors
have precisely the following numerical types:

```text
(-2,div 1), (-4,div 1), (-6,div 3), (-24,div 3),
```

and that a Kähler cone is a connected component of the positive cone after
removing their orthogonal hyperplanes.  Thus, for an exact-Picard OG10
carrier, the complete promotion criterion is intrinsic:

> the proposed Coxeter cone is a Kähler/nef chamber if and only if it is a
> component of the complement of every OG10 wall-divisor hyperplane.

No primitive vector of `L` can have ambient divisibility three.  Indeed, if
`c` did, then `Gc=0 mod 3`; but `det G=-4` is a unit modulo three, so
`c=0 mod 3`, contrary to primitivity.  Therefore only the `-2` and `-4`
branches need be searched inside this `L`.

All 210 input facets pass the wall-divisor gate.  Completeness fails by the
single exact witness `v` above.  Because the primitive embedding orbit is
unique, changing the OG10 marking cannot remove the obstruction.

The precise moduli conclusions are therefore

```text
weak exact-Picard OG10 carrier locus                 nonempty;
exact-Picard OG10 locus with the 210-wall full cone  empty.
```

## 4. What this does—and does not—transfer back

The theorem gives a reusable **promotion test**: enumerate all carrier
embedding orbits, use the deformation-type wall classification, and reject a
proposed complete cone as soon as one wall hyperplane meets its strict
interior.  Here both the embedding enumeration and the obstruction are exact.

It still concerns one polytope and one carrier species.  It therefore does
not imply a new bound for arbitrary arithmetic Coxeter polytopes.  To obtain
the project-level inverse theorem demanded in `GOAL.md`, the same compiler
must be made total for a non-singleton surface-free class, every allowed
carrier species must be exhausted, and the resulting algebraic
classification or emptiness statement must be transferred back to that
whole class.

## 5. Replay

With Alexeev's v4 ancillary file at the default pinned path:

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_borcherds_d4_h21_og10_carrier_obstruction.py \
  --check-archive

PYTHONPATH=src:. pytest -q \
  tests/test_borcherds_d4_h21_og10_carrier_obstruction.py
```

Frozen files:

- certificate:
  `data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json`;
- certifier:
  `scripts/certify_borcherds_d4_h21_og10_carrier_obstruction.py`;
- tests:
  `tests/test_borcherds_d4_h21_og10_carrier_obstruction.py`.
