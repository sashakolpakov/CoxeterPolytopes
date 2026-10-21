# Jacobian--Kummer class_0007: an exact 84-facet polytope in H10

## Result

The class_0007 involution in the full published lattice/chamber stabilizer
`O(S16,D16)` produces an exact finite-volume, noncompact Coxeter polytope in
hyperbolic dimension 10.  Its chamber has 84 facets, 15,698 finite vertices,
and 862 ideal vertices.

This construction uses the divisor lattice and complete induced chamber of a
general Jacobian--Kummer K3 surface.  It does **not** impose the extra condition
that the involution come from a surface automorphism.  The relevant finite
symmetry is a symmetry of the weighted lattice chamber.

The existence and finite-volume statements below are proved twice: first by
Vinberg's complete edge/endpoint and full parabolic-neighborhood criterion,
and independently by exact GMP rational double description of the chamber
cone.  No literature-priority claim is made here; novelty in print requires a
separate bounded literature audit.

## Algebraic and lattice input

The source is the minimal resolution `Y16` of a general Jacobian Kummer
quartic, with divisor lattice `S16 = NS(Y16)`.  Dolgachev--Shimada construct
the induced chamber `D16` from Kondo's embedding into `II_(1,25)`, and
Shimada's exact computational data give:

- the rank-17 Gram matrix `GramS16`, of signature `(1,16)`;
- the interior ample vector `ample16`;
- all 316 walls `D16WallRecs`; and
- the full chamber stabilizer `O(S16,D16)`, of order 23,040 and structure
  `(C2)^5 : S6`.

Primary references are the
[paper](https://doi.org/10.1007/s12215-019-00464-7),
[arXiv version](https://arxiv.org/abs/1908.05390), and
[computational data](https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/15NodalQuart/15nodal/compdata/15nodalcompdata.txt).

The exhaustive GAP subgroup census has 2,455 conjugacy classes.  Class_0007
has:

| datum | exact value |
|---|---:|
| subgroup | `C2` |
| conjugacy-class size | 60 |
| normalizer order | 384 |
| trace on `S16` | 5 |
| fixed rank | 11 |
| hyperbolic dimension | 10 |

The representative matrix is frozen in the screen certificate and replayed
against `GramS16`, `ample16`, and all 316 published walls.

## Exact restricted chamber

The involution has 174 wall orbits: 32 singleton orbits and 142 two-element
orbits.  Summing the wall normals over each orbit spans the complete rank-11
fixed Lorentzian space, with exact inertia `(1,10,163)` for the redundant
174-by-174 orbit Gram matrix.

Exact conic reduction gives:

| restricted-wall outcome | count |
|---|---:|
| genuine facets | 84 |
| redundant inequalities | 90 |

Every removed inequality has an archived positive rational conic dependence.
Every retained inequality has a strict timelike relative-interior witness:
orthogonal projection of `ample16` to that wall is timelike, lies on the
target wall, and is strictly inside all other 173 inequalities.  Thus this is
an exact irredundancy proof, not a floating-point hull heuristic.

Facet-normal squares are:

| square | facets |
|---:|---:|
| `-1` | 64 |
| `-2` | 8 |
| `-4` | 12 |

All 3,486 facet pairs pass the exact crystallographic Coxeter gate:

| relation | pairs |
|---|---:|
| orthogonal (`m=2`) | 1,166 |
| `m=3` | 72 |
| `m=4` | 112 |
| parallel | 860 |
| ultraparallel | 1,276 |
| forbidden normalized product | **0** |

## Intrinsic fixed lattice and mirrors

Returning to the integral published `S16` lattice, a Smith decomposition of
the involution equations gives the saturated fixed lattice

\[
L\cong U\oplus E_8(-2)\oplus\langle-4\rangle .
\]

It has signature `(1,10)`, determinant `1024`, Smith invariants
`1^2,2^8,4`, and discriminant group
`(Z/2)^8 + Z/4`.  The discriminant form splits as the plus-type
eight-dimensional even quadratic space over `F2` orthogonally summed with
`<-1/4>`.  The same split is exhibited for the displayed standard model;
Nikulin's indefinite uniqueness criterion applies at equality
`rank(L)=length(A_L)+2=11`.

The 84 primitive mirrors have intrinsic types

| square and divisibility | mirrors |
|---|---:|
| `(-2,1)` | 8 |
| `(-4,2)` | 33 |
| `(-4,4)` | 43 |

Every one of the 84 corresponding reflections is an integral involutive
isometry of `L`.  Thus literature comparison can use the intrinsic lattice
and mirror data, rather than only the parent Jacobian--Kummer construction.

Exact rational diagonalization and local Hasse comparison, archived in
[`higher_dimensional_rational_lineage_certificate.json`](../data/enumerations/higher_dimensional_rational_lineage_certificate.json),
prove that `L_Q` is similar to the odd-unimodular `I_(1,10)` quadratic space
underlying Vinberg's classical `H10` example.  It is not similar to McLeod's
`f_3^10` space, with local obstructions at 2 and 3.  Hence this 84-wall chamber
is not a new arithmetic commensurability class relative to Vinberg's example;
possible novelty remains only at the finer integral-lattice,
reflection-subgroup, chamber, or algebraic-extraction levels.

## Maximal integral-reflection chamber and classical lineage

An independent all-reflective-root Vinberg search on `L` exhausts the only
possible primitive reflecting norms `-2`, `-4`, and `-8`.  It terminates at
13 simple roots: `1 x (-2,div 1)`, `11 x (-4,div 2)`, and
`1 x (-4,div 4)`.  Exact cone enumeration gives 29 finite vertices, three
cusps, and no hyperideal rays, so this is the complete fundamental chamber of
the full integral reflection group `W(L)`.

After the certified facet permutation, its normalized Coxeter relation
matrix agrees entry by entry with the pinned CoxIter record for Vinberg's
classical

\[
f_2^{10}=-2x_0^2+x_1^2+\cdots+x_{10}^2.
\]

Therefore the 13-wall maximal abstract chamber is known, and the 84 facet
reflections generate a reflection subgroup of exact index 16,450,560 in it.
Solomon's finite-special-subgroup sum gives
`chi(W_84)=-527/11520` and covolume `pi^5*527/340200`; an independent CoxIter
1.1.8 run returns the same invariants and full face vector.  The only remaining
chamber-priority question is whether this particular 84-wall subgroup chamber
was printed before.  See the
[`maximal-reflection lineage audit`](jacobian_kummer_h10_maximal_reflection_lineage.md)
and its
[`exact maximal-chamber certificate`](../data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json)
and
[`index certificate`](../data/enumerations/jacobian_kummer_h10_reflection_index.json).

## Literal Vinberg finite-volume proof

The optimized exact component search found all connected elliptic and affine
subdiagrams.  Its order histograms are:

| order | connected elliptic | connected affine |
|---:|---:|---:|
| 1 | 84 | 0 |
| 2 | 184 | 860 |
| 3 | 456 | 712 |
| 4 | 528 | 606 |
| 5 | 480 | 348 |
| 6 | 384 | 280 |
| 7 | 0 | 144 |

The 2,950 connected affine diagrams yield 862 distinct full zero
neighborhoods.  Every one is componentwise affine of total Coxeter rank 9;
there are zero invalid neighborhoods.

Vinberg's criterion was then applied literally:

| object | exact count |
|---|---:|
| finite vertices (elliptic rank 10) | 15,698 |
| ideal vertices (full parabolic rank 9 neighborhoods) | 862 |
| elliptic rank-9 edges | 107,233 |

Every edge has exactly two endpoints.  Endpoint types are:

| endpoints | edges |
|---|---:|
| finite--finite | 56,560 |
| finite--ideal | 43,860 |
| ideal--ideal | 6,813 |

There are zero bad edges, zero unrecognized endpoint incidences, and zero
invalid parabolic neighborhoods.  Consequently Vinberg's exact decision is
`finite_volume_noncompact`.

## Independent exact-cone proof

The independent certificate reconstructs the fixed cone from the primary
316-wall data and replays all 90 redundancy certificates and all 84 strict
facet witnesses.  It then sends an 84-inequality, rank-11 rational
H-representation to `cddexec_gmp --rep`.

The GMP double-description result has exactly 16,560 primitive extreme rays:

| Lorentz type | rays |
|---|---:|
| timelike (finite vertices) | 15,698 |
| lightlike (ideal vertices) | 862 |
| spacelike/hyperideal | **0** |

All finite rays have exactly 10 active facets.  Ideal rays have 12, 13, or 14
active facets; their 862 incidence sets agree exactly with the 862 exhaustive
affine full neighborhoods.  Reconstructing every rank-9 edge incidence from
the extreme rays reproduces all 107,233 Vinberg edges and the same endpoint
type counts.  This proves independently that the cone has no infinite-volume
hyperideal escape.

## Frozen artifacts and replay

- [`jacobian_kummer_h10_class0007_screen.json`](../data/enumerations/jacobian_kummer_h10_class0007_screen.json)
- [`jacobian_kummer_h10_class0007_component_audit.json`](../data/enumerations/jacobian_kummer_h10_class0007_component_audit.json)
- [`jacobian_kummer_h10_class0007_vinberg_certificate.json`](../data/enumerations/jacobian_kummer_h10_class0007_vinberg_certificate.json)
- [`jacobian_kummer_h10_class0007_exact_cone_certificate.json`](../data/enumerations/jacobian_kummer_h10_class0007_exact_cone_certificate.json)
- [`jacobian_kummer_h10_class0007_exact_cone_rays.jsonl.gz`](../data/enumerations/jacobian_kummer_h10_class0007_exact_cone_rays.jsonl.gz)
- [`jacobian_kummer_h10_class0007_fixed_lattice_certificate.json`](../data/enumerations/jacobian_kummer_h10_class0007_fixed_lattice_certificate.json)
- [`jacobian_kummer_h10_maximal_reflection_lineage.json`](../data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json)
- [`jacobian_kummer_h10_reflection_index.json`](../data/enumerations/jacobian_kummer_h10_reflection_index.json)
- [`screen_jacobian_kummer_h10_class0007.py`](../scripts/screen_jacobian_kummer_h10_class0007.py)
- [`enumerate_jacobian_kummer_h10_class0007_components.py`](../scripts/enumerate_jacobian_kummer_h10_class0007_components.py)
- [`verify_jacobian_kummer_h10_class0007_vinberg.py`](../scripts/verify_jacobian_kummer_h10_class0007_vinberg.py)
- [`independent_certify_jacobian_kummer_h10_class0007.py`](../scripts/independent_certify_jacobian_kummer_h10_class0007.py)
- [`certify_jacobian_kummer_h10_class0007_fixed_lattice.py`](../scripts/certify_jacobian_kummer_h10_class0007_fixed_lattice.py)
- [`certify_jacobian_kummer_h10_maximal_reflection_lineage.py`](../scripts/certify_jacobian_kummer_h10_maximal_reflection_lineage.py)
- [`certify_jacobian_kummer_h10_reflection_index.py`](../scripts/certify_jacobian_kummer_h10_reflection_index.py)

Deterministic archive replay:

```sh
PYTHONPATH=src:. python3 scripts/screen_jacobian_kummer_h10_class0007.py --check-archive
PYTHONPATH=src:. python3 scripts/enumerate_jacobian_kummer_h10_class0007_components.py --check-archive
PYTHONPATH=src:. python3 scripts/verify_jacobian_kummer_h10_class0007_vinberg.py --check-archive
PYTHONPATH=src:. python3 scripts/independent_certify_jacobian_kummer_h10_class0007.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h10_class0007_fixed_lattice.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h10_maximal_reflection_lineage.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_jacobian_kummer_h10_reflection_index.py --check-archive
python3 -m pytest -q tests/test_jacobian_kummer_h10_class0007.py tests/test_jacobian_kummer_h10_class0007_fixed_lattice.py tests/test_jacobian_kummer_h10_maximal_reflection_lineage.py tests/test_jacobian_kummer_h10_reflection_index.py
```
