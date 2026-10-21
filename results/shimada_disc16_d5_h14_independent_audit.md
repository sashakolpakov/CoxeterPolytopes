# Independent audit of the Shimada disc16 D5 fold in H14

## Decision

An intrinsic involution in Shimada's chamber `D5` for the discriminant-16
singular K3 surface produces an exact 60-facet Coxeter polytope in
hyperbolic dimension 14.  The polytope is **finite-volume and noncompact**.

This conclusion does not come from pair compatibility alone.  It is proved by
an independent GMP rational double-description of the complete pointed
chamber cone: all 28,872 extreme rays are timelike or lightlike, and none is
spacelike/hyperideal.

The audit was reconstructed directly from Shimada's official `disc16.txt`
data and did not read the candidate Gram matrix produced by the discovery
search.  The numerical counts agree with that search, but the selection,
walls, irredundancy proof, pair gate, and volume computation were rerun in a
separate script and artifact chain.

No global literature-priority claim is made here.  Whether this 60-facet
chamber, or an isometric copy or finite-index reflection subgroup, has
appeared in print is a separate bounded literature and identification audit.

## Algebraic source

The surface is Shimada's No. 10 singular complex K3 surface, whose oriented
transcendental lattice has Gram matrix

```text
[2 0]
[0 8]
```

and discriminant 16.  Its Néron--Severi lattice therefore has rank 20,
determinant `-16`, and signature `(1,19)`.  Shimada's paper is
[An algorithm to compute automorphism groups of K3 surfaces and an application to singular K3 surfaces](https://arxiv.org/abs/1304.7427);
the exact atlas is the official
[`threediscs.zip`](https://home.hiroshima-u.ac.jp/ichiro-shimada/preprints/AlgoAutK3/threediscs.zip)
supplement.

The audited member is `threediscs/disc16.txt`, with SHA-256
`a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6`.
The reconstruction uses only:

- `GramMatS` and `prMatSdual`;
- `chamberdata[5][weyl]`;
- the archived stabilizer order and generators of `D5`; and
- all `chamberdata[5][wallorbit][*]` records.

## Intrinsic representative and complete restricted walls

The chamber `D5` has an archived stabilizer of order 20 and 165 primitive
dual-basis walls, arranged in 16 published stabilizer orbits.  Enumerating the
entire order-20 stabilizer gives exactly one element with fixed rank 15 on the
Néron--Severi lattice.  It has trace 10 and order two.

The dual action is reconstructed as

\[
G^{-1} M G,
\]

where `G = GramMatS` and `M` is Shimada's row action.  Exact replay verifies
that it preserves the dual Gram matrix, fixes the projected Weyl interior
vector, and permutes all 165 walls.

The wall permutation has 100 orbits:

| orbit size | count |
|---:|---:|
| 1 | 35 |
| 2 | 65 |

Their orbit-sum normals span the entire fixed rank-15 space.  The redundant
100-by-100 Gram matrix has exact inertia `(1,14,85)`.

## Exact irredundancy and Coxeter gate

The projected Weyl vector has square `35/2` and lies strictly inside all 100
restricted inequalities.  Exact reduction gives:

| outcome | count |
|---|---:|
| genuine facets | 60 |
| redundant inequalities | 40 |

Every retained inequality has a strict timelike relative-interior witness:
the Weyl vector projected orthogonally to that wall lies on the wall and is
strictly inside the other 99 inequalities.  Every removed inequality has a
frozen positive rational conic dependence on other orbit sums.

Facet-normal squares are:

| square | facets |
|---:|---:|
| `-1/2` | 20 |
| `-2` | 30 |
| `-4` | 10 |

All 1,770 facet pairs pass the exact crystallographic Coxeter gate:

| relation | pairs |
|---|---:|
| orthogonal (`m=2`) | 865 |
| `m=3` | 35 |
| `m=4` | 70 |
| parallel | 310 |
| ultraparallel | 490 |
| forbidden normalized product | **0** |

This proves that the complete restricted chamber is a hyperbolic Coxeter
polytope candidate.  It does not yet decide its volume.

## Independent GMP volume proof

The independent cone script re-runs the primary-data screen, replays all 40
conic dependencies and all 60 strict facet witnesses, selects 15 independent
facet normals, and constructs the exact 60-by-15 homogeneous
H-representation.  It has rank 15, zero lineality, and the projected Weyl
vector is a strict rational interior point.

`cddexec_gmp --rep` returns exactly 28,872 primitive extreme rays.  Exact
integer replay of every ray against all 60 inequalities gives:

| Lorentz type | extreme rays |
|---|---:|
| timelike / finite vertices | 27,845 |
| lightlike / ideal vertices | 1,027 |
| spacelike / hyperideal | **0** |

All finite vertices have exactly 14 active facets.  The cusp incidences are:

| active facets | cusps |
|---:|---:|
| 15 | 432 |
| 16 | 540 |
| 17 | 55 |

Because the cone is full-dimensional and pointed and every projective
extreme ray lies in the timelike cone or on the light cone, its closure in
the projective hyperbolic ball has no hyperideal escape.  Hence the polytope
has finite volume.  The 1,027 lightlike rays show that it is noncompact.

## Frozen independent artifacts

- [`shimada_disc16_d5_h14_independent_screen.json`](../data/enumerations/shimada_disc16_d5_h14_independent_screen.json)
- [`shimada_disc16_d5_h14_independent_exact_cone_certificate.json`](../data/enumerations/shimada_disc16_d5_h14_independent_exact_cone_certificate.json)
- [`shimada_disc16_d5_h14_independent_exact_cone.ine`](../data/enumerations/shimada_disc16_d5_h14_independent_exact_cone.ine)
- [`shimada_disc16_d5_h14_independent_exact_cone_rays.jsonl.gz`](../data/enumerations/shimada_disc16_d5_h14_independent_exact_cone_rays.jsonl.gz)
- [`independent_screen_shimada_disc16_d5_h14.py`](../scripts/independent_screen_shimada_disc16_d5_h14.py)
- [`independent_certify_shimada_disc16_d5_h14.py`](../scripts/independent_certify_shimada_disc16_d5_h14.py)
- [`test_shimada_disc16_d5_h14_independent.py`](../tests/test_shimada_disc16_d5_h14_independent.py)

With `disc16.txt` at the archived local path, deterministic replay is:

```sh
PYTHONPATH=src:. python3 scripts/independent_screen_shimada_disc16_d5_h14.py --check-archive
PYTHONPATH=src:. python3 scripts/independent_certify_shimada_disc16_d5_h14.py --check-archive
python3 -m pytest -q tests/test_shimada_disc16_d5_h14_independent.py
```

