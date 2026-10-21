# Complete subgroup fishing in Shimada's discriminant-16 chamber `D5`

## Outcome

The complete subgroup-fixed-section census of Shimada's 165-wall chamber
`D5` contains a new positive construction candidate:

> A noncentral involution class in the full weighted-wall stabilizer produces
> a **30-facet finite-volume noncompact Coxeter polytope in `H10`**.

This is no longer a pair-compatible candidate.  Two complete exact volume
audits agree:

- GMP rational double description gives `1,333` finite vertices, `207` ideal
  vertices, and **zero** hyperideal rays;
- a face-driven literal Vinberg audit gives `10,030` elliptic rank-nine
  edges, each with exactly two endpoints, and validates the full rank-nine
  parabolic zero-neighbourhood of every cusp.

The result is stock-bounded.  It exhausts every subgroup of the exact
weighted-wall stabilizer of this one complete `D5` chamber.  It neither
exhausts other Shimada chambers nor makes a global literature-priority claim.
The fixed lattice, maximal reflection group, and weighted-diagram lineage are
being audited separately.

## Exact source and search denominator

The input is the official `threediscs/disc16.txt` supplement to Ichiro
Shimada's [algorithm for automorphism groups of singular K3
surfaces](https://arxiv.org/abs/1304.7427).  The frozen member has SHA-256

```text
a19b797d6c314bc85e1ac4991fe3b190fe29a18a4746103bc6fa9b05e26a7db6
```

The reconstruction uses the rank-20 Néron--Severi Gram matrix, the complete
165-wall `D5` record, its Weyl interior vector, and its archived lattice
stabilizer.  No additional surface-automorphism interpretation is imposed as
a search gate.

An independent exact weighted-graph automorphism enumeration agrees with the
archived lattice stabilizer.  The resulting group is the dihedral group of
order 20, with element orders

| order | elements |
|---:|---:|
| 1 | 1 |
| 2 | 11 |
| 5 | 4 |
| 10 | 4 |

All 22 subgroups, in ten conjugacy classes, are enumerated exactly.  Their
complete order/dimension histogram is:

| subgroup order | fixed hyperbolic dimension | subgroups |
|---:|---:|---:|
| 1 | 19 | 1 |
| 2 | 14 | 1 |
| 2 | 11 | 5 |
| 2 | 10 | 5 |
| 4 | 8 | 5 |
| 5 | 3 | 1 |
| 10 | 3 | 1 |
| 10 | 2 | 2 |
| 20 | 2 | 1 |

Thus the list below really is the full above-`H9` denominator.  Every larger
nontrivial subgroup has fixed dimension at most eight.

## Complete above-`H9` screen

| section | conjugacy-class size | restricted inequalities | facets | forbidden final pairs | result |
|---|---:|---:|---:|---:|---|
| identity `H19` | 1 | 165 | 165 | 2,625 | not Coxeter |
| central involution `H14` | 1 | 100 | 60 | 0 | known Alexeev chamber |
| noncentral involution `H10` | 5 | 86 | 30 | 0 | finite-volume noncompact |
| noncentral involution `H11` | 5 | 91 | 31 | 109 | not Coxeter |

For every nontrivial section, orbit-sum normals span the entire fixed space.
An inequality is removed only after an exact rational nonnegative expression
in the final retained facets is found.  Every retained facet has a strict
timelike relative-interior witness.  In the `H11` class, two orthogonal Weyl
projections initially lie on two additional walls; exact rational tangent
perturbations make both witnesses strict before the 109-pair obstruction is
asserted.

The central `H14` replay gives the previously certified `60 + 40` reduction
and zero bad angles.  That chamber is already identified entrywise with
Alexeev's published `(r,a,delta)=(15,7,1)` chamber, so it is a calibration,
not a new example.

## The `H10` 30-wall chamber

The noncentral fixed-rank-11 class has 86 wall orbits.  Exact global reduction
gives:

| outcome | count |
|---|---:|
| strict facets | 30 |
| redundant inequalities | 56 |

Every redundancy certificate uses only the final 30 retained facets.  The
complete pair gate is:

| Coxeter relation | pairs |
|---|---:|
| orthogonal, `m=2` | 276 |
| `m=3` | 15 |
| `m=4` | 30 |
| parallel | 62 |
| ultraparallel | 52 |
| forbidden | **0** |

This proves that the complete restricted chamber is an acute hyperbolic
Coxeter chamber.  It does not, by itself, prove finite volume.

## Independent GMP finite-volume certificate

Thirty exact inequalities are written in an eleven-dimensional facet-normal
basis.  The H-representation has rank 11, zero lineality, and a strict
rational Weyl interior point.  `cddexec_gmp --rep` returns every projective
extreme ray:

| Lorentz type | rays |
|---|---:|
| timelike / finite | 1,333 |
| lightlike / ideal | 207 |
| spacelike / hyperideal | **0** |

All `1,540 x 30 = 46,200` ray--inequality incidences are replayed exactly.
Every finite vertex has ten active facets.  The ideal active-facet histogram
is:

| active facets | ideal vertices |
|---:|---:|
| 10 | 48 |
| 11 | 98 |
| 12 | 52 |
| 13 | 9 |

Because the full-dimensional pointed projective cone has no hyperideal
extreme ray, the chamber has finite volume.  The 207 lightlike rays prove it
is noncompact.

## Literal Vinberg edge and cusp criterion

The Vinberg cross-check does not infer volume from the preceding counts.  It
reconstructs the relevant Coxeter subdiagrams and exact face incidences.

First, every pair-compatible clique through order 12 is enumerated.  There
are 133,463 pair-compatible rank-nine subsets; exact finite Coxeter
classification leaves exactly 10,030 elliptic diagrams.  Independently, the
complete GMP incidence lattice has exactly the same 10,030 geometric edges.
Every one has two endpoints:

| endpoint types | edges |
|---|---:|
| finite--finite | 4,423 |
| finite--ideal | 4,484 |
| ideal--ideal | 1,123 |

Similarly, exact classification of all 123,508 pair-compatible rank-ten
subsets leaves exactly 1,333 elliptic diagrams, entrywise equal to the finite
vertex incidence sets.

For cusps, the verifier enumerates all connected affine subdiagrams rather
than selecting convenient parabolic subsets.  The exact deletion,
determinant, and rank gates leave 436 affine components.  For every one, its
*full* global zero-neighbourhood is one of the 207 ideal incidence sets.  In
the reverse direction, every ideal incidence set is componentwise affine of
total Coxeter rank nine.  The component-count histogram is:

| affine components at cusp | cusps |
|---:|---:|
| 1 | 48 |
| 2 | 98 |
| 3 | 52 |
| 4 | 9 |

Thus there are zero bad elliptic edges and zero bad or unrecognized
parabolic neighbourhoods.  Vinberg's criterion independently returns
finite-volume noncompact.

## What is and is not new

The proven mathematical output is a new positive construction in this
project's ledger: a complete algebraic wall stock and an exact subgroup
fixed-section mechanism produce a 30-facet finite-volume noncompact Coxeter
polytope in `H10`.

That statement is deliberately weaker than a global novelty claim.  The
remaining lineage questions are:

1. identify the saturated fixed lattice and every primitive reflecting root;
2. determine the maximal reflection chamber and the index, if this is a
   reflection-subgroup chamber;
3. compare the exact weighted diagram against the pinned ordinary and
   reflection-doubling CoxIter corpora and the wider literature; and
4. decide whether this particular divisor-lattice realization has appeared
   before even if the abstract arithmetic group is known.

## Frozen artifacts and replay

- [complete subgroup census](../data/enumerations/shimada_disc16_d5_full_subgroup_fishing.json)
- [GMP cone certificate](../data/enumerations/shimada_disc16_d5_h10_30_exact_cone_certificate.json)
- [exact H-representation](../data/enumerations/shimada_disc16_d5_h10_30_exact_cone.ine)
- [compressed exact rays](../data/enumerations/shimada_disc16_d5_h10_30_exact_cone_rays.jsonl.gz)
- [literal Vinberg certificate](../data/enumerations/shimada_disc16_d5_h10_30_vinberg_face_certificate.json)
- [compressed Vinberg edge audits](../data/enumerations/shimada_disc16_d5_h10_30_vinberg_edges.jsonl.gz)
- [subgroup certifier](../scripts/certify_shimada_disc16_d5_full_subgroup_fishing.py)
- [exact-cone certifier](../scripts/certify_shimada_disc16_d5_h10_30_exact_cone.py)
- [Vinberg certifier](../scripts/certify_shimada_disc16_d5_h10_30_vinberg_faces.py)
- [focused regression tests](../tests/test_shimada_disc16_d5_full_subgroup_fishing.py)

Deterministic replay, with the pinned `disc16.txt` member at its archived
local path and GMP cddlib installed, is:

```sh
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d5_full_subgroup_fishing.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d5_h10_30_exact_cone.py --check-archive
PYTHONPATH=src:. python3 scripts/certify_shimada_disc16_d5_h10_30_vinberg_faces.py --check-archive
python3 -m pytest -q tests/test_shimada_disc16_d5_full_subgroup_fishing.py
```
