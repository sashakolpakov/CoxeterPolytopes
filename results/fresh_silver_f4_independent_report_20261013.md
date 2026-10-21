# Independent silver F4 audit

The saved silver F4 constructions pass an independent exact-arithmetic audit.
The H4 calibration is a compact five-wall simplex; the three H9 finite cones
all have feasible spacelike extreme rays. No new finite-volume example is
certified. Two H9 presentations have the same full integral reflection group
as the earlier silver K3 carrier, as distinct from being integrally isometric
lattices.

The [standalone verifier](../scripts/fresh_silver_f4_independent_verify.py) uses
only Python's `Fraction` and standard-library arithmetic. It imports no project
arithmetic or search code. The [audit output](fresh_silver_f4_independent_audit_20261013.json)
records:

| Host coefficient | Dimension | Walls | Dominant zero roots | Bounded controls | Geometry |
| --- | ---: | ---: | ---: | ---: | --- |
| `lambda=1+sqrt2` | 4 | 5 | 2 | 9 | Five finite vertices; compact simplex |
| `lambda` | 9 | 92 | 6 | 17 | Escape norm `109-77sqrt2>0` |
| `sqrt2` | 9 | 40 | 6 | 19 | Escape norm `92-65sqrt2>0` |
| `2lambda` | 9 | 92 | 7 | 21 | Escape norm `218-154sqrt2>0` |

## Enumeration source review

The positive block is the equal-norm F4 Gram with diagonal 2 and chain
pairings `-1,-sqrt2,-1`. Its determinant is 1. Sylvester's criterion is checked
exactly in both real embeddings. The producer searches dominant covectors
`d=F4*v` with nonnegative primary entries. Cauchy--Schwarz gives
`d_i^2<=2*target` in both embeddings; its integer alphabet encloses both
bounds and filters exactly. The independent checker uses a larger, separately
derived coefficient box.

For a covector prefix of length k, the true minimum quadratic energy after
relaxing the remaining coordinates is
`d_prefix^T * inverse(F4_leading_k) * d_prefix`.
This is the producer's `SCHUR[k-1]`, by the block inverse formula. Rejecting
when that minimum exceeds the target in either embedding is valid. Its second
lower bound uses the leading part of `inverse(F4)` only in the primary
embedding, where all entries and every dominant coordinate are nonnegative.
It is therefore a valid monotone partial-sum bound. No conjugate nonnegativity
of dominant coordinates is assumed.

At full length, unimodularity makes the reconstructed spatial coordinates
integral. The producer enumerates all block combinations and solves the last
scalar square exactly. The time-control bounds are precisely the primary
height bound and the positive-definite conjugate bound. Reflection integrality
is tested in every coordinate before chronological selection. Equal-height
packets are checked mutually acute, and oriented parallel duplicates are
identified before appending.

The checker independently repeats the complete zero-height dominant census
by a direct fourfold product, without the source Schur recursion, and repeats
every time-control census. It does **not** independently repeat every
positive-height spatial recursion; that portion was reviewed at source level.

## Primitive norm classes and zero-height completeness

Write `s=sqrt2`, `eta=2+s`. The ring `Z[s]` is norm-Euclidean, hence a PID.
For a primitive reflective integral normal, its totally positive norm q
satisfies `q | 2 det(B)`: reflection integrality gives `q | 2Br`, and the
adjugate identity plus primitivity gives the determinant divisibility.
The only relevant prime is the ramified prime s. Totally positive units are
squares, and `eta=s*lambda`, `eta^2=2*lambda^2`. Thus representatives for the
possible valuations are `1,eta,2,2eta,4`; the three host coefficients require
respectively the first 3, 4, and 5 classes. This includes all primitive
reflection directions, although the raw candidate census may also contain
nonprimitive vectors.

The independent finite reflection closure has exactly 48 F4 norm-two roots,
with two dominant roots. The direct dominant census finds precisely those
roots in either F4 block, together with scalar roots of norm 1, 2, and (when
allowed) 4. The latter scalar vectors are positive multiples of one another.
Finite F4 reflection reduction takes every zero-height root to a dominant
one, preserving integral reflection conditions. Hence the initial negative
simple-coordinate walls give the complete zero-height arrangement; no extra
mixed-block wall is missing.

## Geometry and Coxeter gates

For every saved wall, the checker recomputes the norm, lattice integrality,
reflection coefficients, chronological height, and strict interior pairing.
For every pair it requires a nonpositive pairing and, below the disjoint or
parallel threshold, the exact condition

`4*(r,t)^2/(r^2*t^2) in {0,1,2,3,2+sqrt2}`.

These are the field's possible Coxeter values with labels 2, 3, 4, 6, and 8.
In particular, the other positive integral-reflection value `2-sqrt2`,
corresponding to angle `3pi/8`, is rejected.

Each H9 escape has positive norm and all wall pairings nonpositive, with
nine independent active facets. Each cone also has a strict timelike
interior. Thus the escape proves infinite volume for the saved finite cone,
not nonreflectivity of the entire lattice.

The H4 calibration is checked by computing all five vertices exactly.
All are timelike and feasible, and the strict interior is a positive
barycentric combination. This certifies the bounded compact simplex. Its
five-cycle diagram has one label-four edge and four label-three edges.

## Audit of the group identifications

The checker reconstructs the unsaturated one-four C5 Gram and its displayed
saturation independently. It replays the saved 5-by-5 integral isometry and
its integral inverse. It then reconstructs the old carrier Gram as saturated
C5 plus F4 plus a norm-two scalar, and verifies its full transport into
`M=<-lambda>+F4^2+<1>`.

The carrier image is exactly the sublattice `L_scalar={z divisible by s}`;
its module index has cardinality 2. It is **not** an integral isometry onto M.
Likewise `L_time={x divisible by s}` is the `alpha=2lambda` host after the
displayed coordinate inclusion.

For an integral vector of M, write a_x,a_z for the constant coefficients of
its time and scalar coordinates modulo two. F4 contributes zero modulo 2,
and `(a+bs)^2=a mod2`, so

`q mod2 = (a_x+a_z)+a_x*s`.

Thus L_time and L_scalar are intrinsic subsets of M determined by the norm:
the former has even s-coefficient of q, and the latter has even sum of the
two coefficients. Every integral isometry of M preserves each sublattice.
The independent checker repeats the residue identity on all time/scalar
residue classes.

Conversely each smaller lattice has cyclic discriminant module `O/(2)`:
all its other blocks are unimodular. This module has exactly one nonzero
proper O-submodule, of order 2, independently checked by enumerating the
four-element residue ring's ideals. The associated integral overclosure is
M. The whole dual is not integral, since its exceptional dual generator has
norm `-1/(2lambda)` or `1/2`. Every isometry of the smaller lattice must
preserve the unique overclosure M. Consequently the full integral orthogonal
groups agree under the embeddings, and hence so do the full reflection
groups. This proof is valid and does not require an integral isometry of the
three lattices.

The two 92-wall lists also agree exactly as oriented hyperplanes after the
time scaling. The previous K3 code searched only its norm-two subgroup;
this identification concerns the full groups, not equality with that
previously searched subgroup.

Finally the displayed F4-to-diagonal field isometry is independently
replayed. Its nonintegral denominator prevents an integral lattice or
full-group identification with the diagonal lattice. The `alpha=sqrt2`
case therefore is not established to be a new group, nor proved to duplicate
the full group of that diagonal lattice.

## Validation

[Focused regressions](../tests/test_fresh_silver_f4_independent.py) compare the
producer's embedding bounds and zero-block recursion against independent
brute-force calculations, reject a forged escape norm, reject angle `3pi/8`,
and check the compact calibration.

```sh
PYTHONPATH=. python3 scripts/fresh_silver_f4_independent_verify.py \
  --output results/fresh_silver_f4_independent_audit_20261013.json
PYTHONPATH=src:. python3 -m pytest -q tests/test_fresh_silver_f4_independent.py \
  tests/test_fresh_rank3_h9.py
```

## Why the earlier norm-two stabilizer proof does not transfer

The old certificate's involution `H` swaps walls 4 and 18 and is reflection
in their difference d, of norm `4+2sqrt2`. Transport d to M using the
verified carrier embedding. It is divisible by sqrt2 in M. After orientation,
the primitive normal `r=-T*d/sqrt2` is

```
(3+2s, 6+4s, 9+6s, 8+6s, 4+3s,
 2+2s, 3+3s, 4+2s, 2+s, 1+s),  s=sqrt2.
```

Its norm is eta=`2+s`, and every reflection coefficient is integral. The
transported H is exactly reflection in this normal. The former H-fixed
regular norm-two witness is on this additional mirror, so it is not regular
for the full reflective arrangement.

This is **already wall 46** (zero-based) of the 92-wall lambda cone, at
height `12+17s/2<30`. Its pairing with the saved escape is `-7+2s<0`, so
it supplies no new cut. Its pairing with the strict interior is
`-17/8-3s/2<0`; applying H reverses that pairing. Hence H itself exchanges
the two halfspaces of this full-alphabet mirror, rather than preserving a
full arrangement chamber. No separator search is needed to expose this
specific failure of transferring the norm-two proof.

This says nothing by itself about whether the old **infinite-order product**
could stabilize some full chamber. Full-host reflectivity remains open.

The [standalone check](../scripts/fresh_silver_f4_stabilizer_check.py) replays
the transported normal, reflection matrix, fixed witness, wall-list identity,
and the strict side exchange. The [certificate](fresh_silver_f4_stabilizer_check_20261013.json)
contains the literal normal and matrix.

```sh
PYTHONPATH=. python3 scripts/fresh_silver_f4_stabilizer_check.py \
  --output results/fresh_silver_f4_stabilizer_check_20261013.json
```
