# TOOT: a three-parameter `D_infinity^3` family in H12

Status: **all exact checks pass** (2026-08-08).

## Positive theorem

The finite-volume noncompact 172-facet Jacobian–Kummer seed polytope in
`H^12` contains three mutually commuting infinite-dihedral doubling factors.
For every triple of positive integers `(N1,N2,N3)` and every choice of the
three path orientations, their finite chamber box is a finite-volume
noncompact hyperbolic Coxeter polytope containing exactly `N1*N2*N3` seed
chambers.  Hence

```text
index  = N1*N2*N3,
volume = N1*N2*N3 * Vol(P_seed).
```

Its exact facet count is

```text
f = [177 N1 N2 N3
     + 29(N1 N2 + N1 N3 + N2 N3)
     + 13(N1 + N2 + N3)
     + 40 + epsilon] / 2.
```

The phase convention is essential.  Let `sigma` be the parity of the three
endpoint-orientation flips, and let

```text
E = 0  if N1,N2,N3 are all even,
E = 1  otherwise.
```

Then

```text
epsilon = E      if sigma is even,
epsilon = 2 - E  if sigma is odd.
```

Thus all phases have the same count whenever at least one `Ni` is odd.  When
all three are even, the two phase-parity classes differ by exactly one facet.
For example,

```text
(N1,N2,N3)=(2,2,2):  941 facets in even phase,
                      942 facets in odd phase.
```

This explicitly corrects the preliminary phase-blind formula.  The test
suite contains a regression assertion that `941 != 942`.

Every boundary normal is a primitive integral reflection root in the same
rank-13 lattice.  Every member retains the weak full-Picard K3 carrier

```text
U(2) + D5(-1)^2 + <-4>.
```

No full-nef or surface-automorphism realization is claimed.

## The two exact configurations

The independent 124-wall doubling-stock search finds precisely two direct
`D_infinity^3` configurations.  Wall numbers are one based:

```text
A: (W48,W42) x (W51,W39) x (W56,W36),
B: (W58,W33) x (W59,W32) x (W60,W31).
```

For either row the selected six-wall Gram matrix is

```text
[-4  4  0  0  0  0]
[ 4 -4  0  0  0  0]
[ 0  0 -4  4  0  0]
[ 0  0  4 -4  0  0]
[ 0  0  0  0 -4  4]
[ 0  0  0  0  4 -4].
```

All six reflections are integral involutive lattice isometries.  The three
parallel factors commute pairwise.  Their primitive isotropic half-sums are
the same vector, so the three parabolic translations fix a common cusp:

```text
u=(3,3,3,3,2,4,-6,-7,-6,-7,-5,-6,-8)
```

up to sign.  If `Ti=I+Ki`, exact arithmetic gives

```text
rank(Ki)=2, rank(Ki^2)=1, Ki^3=0,
Ki Kj=0 for i != j,
image(Ki^2)=Q u.
```

This proves that the factors are genuinely parallel, not merely disjoint or
ultraparallel, and supplies the algebra used in the universal collision
proof.

## Exhaustive all-`Z^3` collision proof

Remove the six selected endpoint walls.  The remaining 166 seed facets are
classified independently relative to each parallel pair:

```text
B = orthogonal to both endpoints,
L = orthogonal only to the left endpoint,
R = orthogonal only to the right endpoint,
N = orthogonal to neither endpoint.
```

There are 12 `BBB` facets.  No other pattern contains `B`.  The remaining
154 facets occupy 27 exact `L/R/N` patterns; the complete histogram is frozen
in the JSON certificate.

The universal proof considers

```text
T1^z1 T2^z2 T3^z3 S_R1^b1 S_R2^b2 S_R3^b3(r_i) = +/- r_j
```

for every `(z1,z2,z3) in Z^3`, all eight branch bits, both signs, and every
ordered pair of the 166 side facets.  This is 440,896 signed equations.  No
finite parameter box is imposed.

For each of the 1,232 moving branch/source pairs, the four vectors

```text
[K1 v, K2 v, K3 v, u]
```

have rank four.  Modulo the common null line this gives the unique rational
candidate `(z1,z2,z3)`.  The verifier tests its integrality and then all 13
coordinates exactly.  The complete result is

```text
462 finite solutions, all positive self-stabilizers;
 96 all-Z^3 identities, exactly 12 BBB facets in each branch;
  0 negative-root collisions;
  0 cross-facet collisions.
```

The six endpoint orbits are audited separately.  Pairing with `u` separates
them from every moving side orbit; all different-factor endpoint classes are
distinct modulo `u`; and the 12 fixed `BBB` roots are distinct from all six
endpoint roots.  The one-factor endpoint equations leave exactly two outer
walls per factor, hence six endpoints in every finite box.

## Derivation of the formula

For one factor and a path of length `N`, a seed facet contributes

```text
w_B(N)=1,
w_L(N)=ceil(N/2),
w_R(N)=floor(N/2)+1,
w_N(N)=N.
```

Reversing that factor exchanges `w_L` and `w_R`.  The all-`Z^3` theorem says
that there are no additional identifications between different seed facets,
so the side-facet count is the exact sum of the products of these weights
over the frozen 28-pattern histogram.  Adding the six outer endpoints and
simplifying in all 64 combinations of parameter parity and phase bits gives
the displayed closed formula and epsilon convention.  Every one of the 64
symbolic differences is identically zero.

As a redundant finite replay, the certifier constructs all eight phases for
all 27 triples in `{1,2,3}^3`, plus `(4,4,4)`, `(2,3,4)`, and `(5,2,3)`: 240
exact boundaries in total.  All match the phase-refined formula, and every
displayed root remains primitive and integrally reflective.  Sample counts
are

```text
f(1,1,1)=172,
f(2,1,1)=296,
f(2,2,1)=523,
f_even(2,2,2)=941,
f_odd (2,2,2)=942,
f(3,3,3)=2860,
f_even(4,4,4)=6458.
```

## Symmetry and the phase obstruction

The exact published `D16` chamber group has order 23,040.  The centralizer of
the trace-nine involution has order 768.  It carries configuration A to each
of A and B in exactly 384 ways, so the two six-wall configurations form one
exact lattice-isometry orbit.

The stabilizer of A has order 384 and realizes all six permutations of the
three factors, but only 24 of the 48 possible labelled factor/endpoint
actions: precisely the even endpoint-flip subgroup, with 16 lattice lifts per
labelled action.  This is the exact symmetry reason the phase parity cannot
be discarded.  When all `Ni` are even, the two phase classes are certainly
nonisometric because their facet counts differ by one.  When some `Ni` is
odd, their facet counts agree; no wider isometry identification is claimed
here.

## Scoped maximum-product theorem

The exhaustive relation graph on the 124 doubling walls has

```text
 264 orthogonal pairs,
1138 parallel pairs,
6224 ultraparallel pairs.
```

The graph whose vertices are parallel pairs and whose edges mean complete
cross-orthogonality has only six edges.  Its maximal factor cliques are
exactly the two triples above.  Consequently, inside this fixed doubling
stock:

- the maximum number of mutually orthogonal `D_infinity` factors is three;
- there is no `D_infinity^4` product;
- no additional `C2` factor is orthogonal to either `D_infinity^3` triple;
- the maximum `D_infinity^2 x C2^k` has `k=1`;
- the maximum `D_infinity x C2^k` has `k=3`, attained by exactly the eight
  already frozen path-cube configurations;
- the maximum pure commuting cube is `C2^4`, attained by 16 labelled
  four-wall sets.

This theorem is deliberately scoped to direct products assembled from
parallel doubling-wall pairs and pairwise orthogonal factors in this exact
124-wall stock.  It is not a classification of arbitrary non-product
redoubling trees.

## Finite volume and the full-Picard carrier

Finite volume is inherited from the independently certified seed because a
box is the union of exactly `N1*N2*N3` interiors-disjoint seed chambers.  It
remains noncompact: it contains a closed noncompact seed chamber and retains
the common ideal direction.  This uses the literal Vinberg
edge/endpoint/parabolic-neighbourhood certificate and the independent exact
GMP extreme-ray certificate, not the K3 lattice alone.

Every box boundary contains the 166 nonendpoint seed roots and one outer
endpoint root from each factor.  All eight resulting 169-root phase cores
have Smith diagonal `1^13`; hence every descendant spans the full rank-13
Picard carrier.  Integral isometries transport primitive seed reflection
roots, so primitivity and integral reflectivity persist for all parameters.

The carrier is a seven-dimensional very-general projective K3 family with
full Picard lattice `U(2)+D5(-1)^2+<-4>`.  This remains a weak divisor-
reflection carrier, not a nef-cone or surface-automorphism realization.

## Lineage and claim boundary

The pinned CoxIter ordinary and index-two lanes have H12 facet counts at most
35.  Neither the 172-facet seed nor any nonseed descendant, whose first count
is 296, occurs by literal H12 header.  This is bounded same-header absence
only: those lanes are not a publication-complete denominator, and worldwide
publication priority remains open.

The project-level positive advance is the exact `D_infinity^3`
three-parameter family, its phase-refined formula, the all-`Z^3` proof, and
the scoped maximum-product theorem.

## Reproduction

```bash
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h12_triple_parallel_family.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_triple_parallel_family.py
```

Artifacts:

- [`jacobian_kummer_h12_triple_parallel_family.json`](../data/enumerations/jacobian_kummer_h12_triple_parallel_family.json)
- [`certify_jacobian_kummer_h12_triple_parallel_family.py`](../scripts/certify_jacobian_kummer_h12_triple_parallel_family.py)
- [`test_jacobian_kummer_h12_triple_parallel_family.py`](../tests/test_jacobian_kummer_h12_triple_parallel_family.py)

Fresh verification: `5 passed in 94.63s`.

Frozen SHA-256 values:

```text
script       cd27d2201e2ccbe5cb665637fe225bfbe40c62554a41f4290f09cbd1978f135e
certificate  bf6fa667f224df5117b97a4590a1f673b0504d3144dfd573219e177cec487c5a
test         1ee981a2898eb30067ba7a71c39070704582861ae5bd9f28455dd1614b2a6a04
```
