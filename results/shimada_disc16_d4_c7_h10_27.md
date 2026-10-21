# The archived `D0004/C007` 27-wall chamber in `H10`

## Outcome

An order-four subgroup of the archived stabilizer of Shimada's
discriminant-16 chamber `D0004` produces a **27-facet, finite-volume,
noncompact hyperbolic Coxeter polytope in `H10`**.

This has passed two complete exact volume audits:

- GMP rational double description gives 895 finite vertices, 67 ideal
  vertices, and **zero** hyperideal rays;
- a literal face-driven Vinberg audit gives 6,221 elliptic rank-nine edges,
  every one with exactly two endpoints, and validates the full rank-nine
  parabolic zero-neighbourhood of every cusp.

This is a new positive construction in the project ledger.  Its exact
lattice, reflection-group, subgroup-index, and bounded-corpus lineage is now
closed in the
[lineage theorem](shimada_disc16_d4_c7_h10_27_lineage.md).  It is not a new
stable lattice or maximal reflection group, and no worldwide priority claim
is made for the literal 27-wall subgroup chamber.

## Algebraic and search provenance

The source is the rank-20 Néron--Severi lattice and complete published wall
stock of Shimada's singular K3 chamber `D0004`.  The fixed section comes from
one conjugacy class of order-four subgroups of the chamber's archived
stabilizer `Aut(D)`:

| datum | value |
|---|---:|
| archived stabilizer order | 48 |
| subgroup order | 4 |
| subgroup conjugacy-class size | 3 |
| saturated fixed-lattice rank | 11 |
| fixed-lattice inertia | `(1,10,0)` |
| fixed-lattice determinant | 512 |

The construction uses the divisor lattice as the source of hyperbolicity:
the Hodge-index signature on the complete saturated fixed lattice is exactly
the Lorentz signature needed for `H10`.  No surface-automorphism condition is
used as an additional gate.

The construction certificate starts from Shimada's **archived** chamber
stabilizer.  The subsequently completed atlas-wide weighted-wall audit proves
that for chamber `D0004` the archived group and the full exact weighted-wall
group are equal, both of order 48.  Thus no additional chamber symmetry is
missing from this source stock.

## Exact wall reduction and Coxeter gate

The 158 published walls fold to 62 wall orbits.  Exact conic reduction in the
saturated fixed lattice gives:

| result | count |
|---|---:|
| strict hyperbolic facets | 27 |
| exact conic redundancies | 35 |
| unresolved inequalities | 0 |

Every removed inequality is replayed as a nonnegative rational sum of the
final 27 facets.  Every retained facet has an exact timelike relative-
interior witness which is zero on that facet and strictly positive on all
other 61 orbit inequalities.

All 351 final facet pairs are checked exactly:

| Coxeter relation | pairs |
|---|---:|
| orthogonal, `m=2` | 235 |
| `m=3` | 15 |
| `m=4` | 14 |
| parallel | 48 |
| ultraparallel | 39 |
| forbidden | **0** |

This proves that the complete restricted chamber is an acute Coxeter
chamber.  The pair gate alone does not establish finite volume.

## Complete GMP cone audit

The 27 exact inequalities form a rank-11, full-dimensional pointed
homogeneous cone.  The projected ample class has minimum strict wall pairing
one.  `cddexec_gmp --rep` enumerates every projective extreme ray over
`Q`:

| Lorentz type | rays |
|---|---:|
| timelike / finite | 895 |
| lightlike / ideal | 67 |
| spacelike / hyperideal | **0** |

All `962 x 27 = 25,974` exact ray--inequality incidences are replayed.  The
active-facet histogram is:

| ray type | active facets | rays |
|---|---:|---:|
| finite | 10 | 895 |
| ideal | 11 | 6 |
| ideal | 12 | 58 |
| ideal | 13 | 3 |

The absence of a hyperideal extreme ray proves finite volume for the complete
pointed projective cone.  The 67 lightlike rays prove noncompactness.

## Literal Vinberg audit

The diagram verifier independently enumerates every exact pair-compatible
clique through order 12.  Among 55,451 compatible rank-nine subsets, exact
finite Coxeter classification leaves 6,221 elliptic edge diagrams.  These
are entrywise equal to the 6,221 geometric edges reconstructed from the GMP
face incidences, and each edge has exactly two endpoints:

| endpoint types | edges |
|---|---:|
| finite--finite | 3,161 |
| finite--ideal | 2,628 |
| ideal--ideal | 432 |

At rank ten, the 895 elliptic diagrams are entrywise equal to the 895 finite
vertex incidence sets.

For the cusp half of Vinberg's criterion, all connected affine subdiagrams
are enumerated using exact deletion, determinant, and rank tests.  There are
198.  The full global zero-neighbourhood of every component is one of the 67
certified ideal incidence sets.  Conversely, every ideal incidence set is
componentwise affine of total Coxeter rank nine:

| affine components at a cusp | cusps |
|---:|---:|
| 2 | 6 |
| 3 | 58 |
| 4 | 3 |

Hence the literal Vinberg audit has zero bad elliptic edges and zero bad
parabolic neighbourhoods, independently returning finite-volume
noncompact.

## Closed lineage and remaining priority question

The 27 primitive roots span the determinant-512 fixed lattice.  Their
reflections preserve a unique index-four integral overlattice, explicitly
`U + D6(-1) + A1(-1)^3`, whose full reflection chamber is the known 14-wall
Scharlau--Kirschmer `List10/Gram10` entry 20.  Complete Solomon enumeration
gives exact subgroup index `2016` and covolume
`pi^5 * 527/48,988,800`.  Neither pinned CoxIter corpus contains any
27-facet `H10` record.  The only remaining lineage issue is worldwide
priority for the literal subgroup chamber or this divisor-fixed realization;
see the [exact lineage report](shimada_disc16_d4_c7_h10_27_lineage.md).

## Frozen artifacts and replay

- [archived fixed-section case](../data/enumerations/shimada_disc16_atlas_archived_high_sections/D0004_C007.json)
- [GMP cone certificate](../data/enumerations/shimada_disc16_d4_c7_h10_27_exact_cone_certificate.json)
- [exact H-representation](../data/enumerations/shimada_disc16_d4_c7_h10_27_exact_cone.ine)
- [compressed exact rays](../data/enumerations/shimada_disc16_d4_c7_h10_27_exact_cone_rays.jsonl.gz)
- [literal Vinberg certificate](../data/enumerations/shimada_disc16_d4_c7_h10_27_vinberg_face_certificate.json)
- [compressed Vinberg edge audit](../data/enumerations/shimada_disc16_d4_c7_h10_27_vinberg_edges.jsonl.gz)
- [exact-cone certifier](../scripts/certify_shimada_disc16_d4_c7_h10_27_exact_cone.py)
- [Vinberg certifier](../scripts/certify_shimada_disc16_d4_c7_h10_27_vinberg_faces.py)
- [focused regression tests](../tests/test_shimada_disc16_d4_c7_h10_27.py)
- [lineage certificate](../data/enumerations/shimada_disc16_d4_c7_h10_27_lineage.json)
- [lineage certifier](../scripts/certify_shimada_disc16_d4_c7_h10_lineage.py)
- [lineage report](shimada_disc16_d4_c7_h10_27_lineage.md)

Deterministic replay, with GMP cddlib installed, is:

```sh
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d4_c7_h10_27_exact_cone.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d4_c7_h10_27_vinberg_faces.py --check-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_shimada_disc16_d4_c7_h10_27.py
```
