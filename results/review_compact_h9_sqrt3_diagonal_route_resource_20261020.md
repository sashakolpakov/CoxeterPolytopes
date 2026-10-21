# Independent audit: square-root-of-three diagonal-edge construction route

20 October 2026. Source/paper review only. No execution, import of the
scientific worker, enumeration, replay, build or Lean run.

Reviewed source:
`results/compact_h9_sqrt3_diagonal_construction_route_20261020.md`,
SHA256 `1ac6604b8760f5406bfc679fe362b8274c1542165c6e32136cf4f84d85b6b96b`.

**Verdict: PASS as a complete necessary endpoint inventory and bounded
first-edge search reduction, with one harmless wording correction below.**
No compactness, finite volume, actual first simple mirror or completion
of a chamber has yet been proved by this source.

## Integral host and primitive normalization

I independently checked the norm-Euclidean argument: simultaneous
rounding in Z[sqrt(3)] bounds the absolute norm of a remainder by 3/4.
Thus primitive coordinates have a Bezout relation. Reflection
integrality first implies that q divides each 2(G r)_i, by multiplying
the integral rank-one reflection entries by that Bezout relation.
Multiplying the relation again by 2sqrt(3) gives q dividing 2sqrt(3).
This justifies the primitive divisor bound rather than assuming it.

The only ramified prime elements needed are eta=1+sqrt(3), of norm -2,
and s=sqrt(3), of norm -3, with eta^2=2epsilon. At the second embedding
both change sign; all positive units have positive conjugate since a
norm-minus-one unit would give a square equal to -1 modulo 3. Total
positivity therefore leaves valuation pairs (0,0), (2,0), (1,1).
Unit rescaling of a primitive vector changes its norm by a unit square.
Consequently the six displayed representatives are a complete safe
necessary pool. They are not asserted all to occur.

Wording correction, preserving the pinned source: for a positive unit
1<=u<epsilon, its s coefficient is in [0,1), hence is the INTEGER zero,
and u=1. For 1<u<epsilon it is strictly between zero and one, a
contradiction. The source's strict interval at the endpoint and its
antecedent "it is rational and equals1" should be read with this
correction. The claimed unit group and all subsequent algebra remain
valid.

The two omitted examples really are integral reflectors. Their squares
3-s and 2-s normalize respectively to (3+s)epsilon and epsilon by unit
rescaling. Hence the old norms-one/two screen cannot exclude this host.
The nonphysical embedding is positive definite, so anisotropy and the
conditional assertion that a finite-volume integral-reflection chamber
would be compact are also valid.

## Zero chamber and all finite A8 endpoint profiles

The B9 zero-height chamber is complete. One can even bypass the cited
finite-reflection enlargement argument: an integral spatial isometry
has each entry bounded in absolute value by one at both embeddings.
For an entry A+B s, this implies |B|<=1/s<1, so B=0. Orthogonality
then makes the matrix a signed permutation matrix.

At a first finite diagonal endpoint the eight actual A8 simple walls
remain present. The finite rank-nine classification gives exactly
A8+A1, A9, B9 or D9. A nonzero contact is one gap: a3 at path end for
A9, a3 at the next position for D9, or a4 at an end for B9. Other
finite exceptional types have insufficient rank and cannot contain
this eight-node induced path.

For a3 the required integral gap square is q/2. Of the pool, only q=2
works: epsilon/2 has square root (1+s)/2 in the field but not in the
integer ring. For a4, only q=1 and q=2epsilon have integral square
roots, respectively 1 and 1+s. Nonsquare rational norms dispose of
the last two pool classes, and the unit-square classification handles
epsilon. Thus no hidden long-root or extra finite label is omitted.

Writing z=9b+kg and X=3a gives exactly

    z^2-s X^2=9q-k(9-k)g^2.

The isolated profiles give six right sides 9q; the four nonisolated
rows give 10, 4, 1, 2epsilon, with their two marked placements. The
congruences X in 3O and z-kg in 9O recover all coordinates, not merely
the scalar norm. All these right sides are positive at both embeddings.

## Explicit finite upper reflector

I expanded U^3 directly and obtained the displayed coefficients
(98+57s)+(75+43s)alpha. Cubing that expression modulo 9 gives
Z=8 and X=3s modulo 9, so a=X/3 and b=(Z+1)/9 are integral.
From Z^2-sX^2=1,

    9b^2-s a^2=2b,

and the reflection in rho has precisely the two displayed integral
coordinate formulas. Its square 2b is positive. Primitivity of this
particular presentation is unnecessary for an integral reflection.
Its diagonal collision satisfies lambda_*^2<s/9, so it is genuinely
finite rather than an ideal or spacelike cut.

The entire integral isometry group is discrete at the physical
embedding: compact bounds there, together with bounded matrices in
the positive-definite conjugate group, give finitely many integral
matrix entries. Its reflecting mirrors are locally finite. Thus the
known finite collision supplies a true upper endpoint for the first
collision, without proving that this particular upper mirror is simple.

## Orbit representative boxes and the five powers

For any physical positive solution, w=z+alpha X and R/w=z-alpha X
are positive. Multiplication by U or its inverse preserves O-pairs
and R. Reducing w/sqrt(R) to [1,U) makes z positive and X nonnegative,
with the exact upper bounds

    z<(2+s)sqrt(R),        X<(1+s)sqrt(R).

At the conjugate embedding z'^2+sX'^2=R', giving both finite conjugate
bounds. Together these bound all four integer coefficients of z and X.
No assertion that U generates the full relative unit group is needed:
every solution belongs to some orbit of this explicitly integral unit,
and each such orbit has a representative in the strip being enumerated.

For the explicit upper mirror, the normalized w ratio is exactly
U^(9/2). The collision X/z, and hence lambda=sX/(3z), increases strictly
with w/sqrt(R). A reduced representative has ratio in [1,U); powers
j<0 have X<0, while every j>=5 is beyond the upper mirror. Thus testing
j=0,...,4 is exhaustive in the declared comparison range. Congruences,
primitive normalization and integral reflection must still be checked
after each power; they are not assumed preserved by U.

The exact comparison X(Z+1)<=X_9 z follows by cross multiplication of
positive denominators, so no algebraic square-root approximation or
assumed fundamental relative unit enters the cutoff.

## Boundaries of this PASS

The finite Coxeter classification and elementary arithmetic of the
integer ring are explicit mathematical inputs. The note does not
claim that its future implementation has enumerated the boxes, that
an arbitrary mirror at the least displayed collision is the actual
simple wall, or that closing this edge proves a finite-volume chamber.
Its warning about simultaneous incident mirrors and the actual finite
sector is necessary and correctly retained. Full cone/vertex coverage
would still be required for any compact H9 example.

No scientific worker source is included in this verdict; root owns its
separate review and serialized execution.
