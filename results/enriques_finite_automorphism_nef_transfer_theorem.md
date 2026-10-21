# The finite-automorphism Enriques seven-image transfer theorem

**Frozen scope:** complex Enriques surfaces with finite automorphism group,
their complete sets of smooth rational curves, and the resulting full nef
chambers in `Num(Y)`.  This is a strong surface realization of known
Nikulin--Kondō configurations, not a worldwide novelty claim for the surfaces
or reflection groups.

## The theorem

Let `M_Enr^fin` be the moduli problem of numerically marked complex Enriques
surfaces with finite automorphism group.  A marking identifies

\[
 \operatorname{Num}(Y)\cong L_{10}=U\oplus E_8(-1),
 \qquad \operatorname{sig}(L_{10})=(1,9).
\]

For such a surface let

\[
 P(Y)=\bigl\{[x]\in \mathbb H^9(L_{10}):
          x\cdot C\geq 0\text{ for every smooth rational curve }C\subset Y
       \bigr\}.
\]

Then the forgetful map

\[
 \Phi_{\rm Enr}:M_{\rm Enr}^{\rm fin}\longrightarrow
 \{\text{hyperbolic Coxeter polytopes up to congruence}\},
 \qquad (Y,\eta)\longmapsto P(Y),
\]

has exactly seven images.  They are the Nikulin--Kondō types I--VII, with
Shimada records

```text
I=12A, II=12B, III=20B, IV=20F,
V=20A, VI=20E, VII=20C=20D.
```

Every image is a genuine finite-volume noncompact Coxeter polytope in
`H9`.  Types I and II have 12 facets; Types III--VII have 20 facets.  In
particular, the full-nef finite-automorphism Enriques class contains no
`H9/14` polytope.  This last statement is scoped to this carrier class and is
not a global `H9/14` nonexistence theorem.

There is also a surface-free converse formulation.  Let `C_Enr` be the class
of Coxeter polytopes whose primitive integral normalized wall Gram is one of
the seven explicit matrices in the machine certificate.  Every member of
`C_Enr` is congruent to one of the seven chambers above and therefore has a
full-nef complex-Enriques carrier.  Conversely, every full-nef chamber in the
stated surface class lies in `C_Enr`.  Thus

\[
 P\in C_{\rm Enr}
 \quad\Longleftrightarrow\quad
 \Phi_{\rm Enr}^{-1}(P)\ne\varnothing .
\]

This is the first exact finite calibration of the proposed
polytope-to-surface transfer principle.  It is a lookup criterion for seven
explicit integral Gram types, not yet the desired broad universality theorem
for an infinite class of Coxeter polytopes.

The finite-automorphism assumption is subsequently shown to follow from
finite hyperbolic volume of the complete nef chamber.  Thus these same seven
types exhaust the full-nef finite-volume outputs of **all** complex Enriques
surfaces; see the
[full-nef classification theorem](enriques_full_nef_finite_volume_classification.md).

## Exact polytope census

For roots of square `-2`, an intersection `0` gives Coxeter label `m=2`, an
intersection `1` gives `m=3`, and an intersection `2` gives parallel walls.
No selected configuration has an ultraparallel pair.  `F` and `I` below are
finite and ideal vertices; “affine audits” counts every connected affine
component whose full zero-neighborhood was tested.

| type | Shimada | facets | `m=3` pairs | parallel pairs | `F` | cusps `I` | elliptic rank-8 edges | affine audits |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| I | 12A | 12 | 10 | 3 | 14 | 9 | 114 | 12 |
| II | 12B | 12 | 15 | 0 | 38 | 17 | 270 | 20 |
| III | 20B | 20 | 16 | 36 | 200 | 44 | 1,301 | 84 |
| IV | 20F | 20 | 40 | 10 | 720 | 87 | 4,405 | 184 |
| V | 20A | 20 | 24 | 31 | 192 | 29 | 1,182 | 68 |
| VI | 20E | 20 | 45 | 10 | 600 | 57 | 3,420 | 124 |
| VII | 20C | 20 | 30 | 25 | 450 | 51 | 2,490 | 92 |

The two one-dimensional surface families, Types I and II, therefore map to
single polytope types.  Types III--VII are isolated complex surfaces in the
classification.  Markings still make the forgetful map many-to-one.

## Why these are full nef chambers

There are four separate gates.

1. **Free hyperbolicity.**  The numerical divisor lattice of an Enriques
   surface is the even unimodular lattice `L10` of signature `(1,9)`, so its
   positive rays form `H9`.
2. **Complete geometric walls.**  Kondō's classification gives the dual graph
   of *all* smooth rational curves: 12 in Types I--II and 20 in Types
   III--VII.  Martin's Main Theorem restates the classification as an iff for
   the complete graph.  A negative irreducible curve spans an extremal ray of
   the effective cone, hence its dual nef inequality is a real facet.  There
   is no substitution of a convenient curve subconfiguration for the full
   wall set.
3. **Coxeter gate.**  Shimada's exact vectors replay as
   `configuration = W Gram(L10) W^t`; every vector is primitive of square
   `-2`, the vectors span rank ten, and every off-diagonal entry is `0`, `1`,
   or `2`.  The positive-normal Gram has inertia `(9,1,f-10)`.
4. **Volume gate.**  The exact convex-cone and literal Vinberg audits below
   both pass.  Finite automorphism or finite reflection index is not used as
   an unexamined synonym for finite hyperbolic volume.

Kondō's Proposition 1.10 supplies an additional no-hidden-curve implication:
once the displayed reflection group has finite index, the displayed nodal
curves are the only smooth rational curves.  Here this agrees with the
classification's already complete curve lists.

## Two independent finite-volume certificates

### Exact double description

For each type, ten independent wall vectors give a rational homogeneous
`H`-representation.  `cddexec_gmp` enumerates every extreme ray.  Across the
seven cones it finds exactly 2,508 rays:

```text
finite = 14+38+200+720+192+600+450 = 2,214
ideal  =  9+17+ 44+ 87+ 29+ 57+ 51 =   294
hyperideal = 0
```

Every ray satisfies every inequality exactly, every causal ray has positive
Lorentz pairing with the strict ample witness, and every proposed wall has
active-ray span nine in the ten-dimensional homogeneous cone.  Thus the cone
is pointed and full-dimensional, all walls are irredundant, all projective
vertices lie in the closed hyperbolic ball, and every chamber is
finite-volume.  The 294 null rays prove noncompactness.

### Literal Vinberg criterion

Independently, the verifier enumerates every elliptic rank-eight edge diagram,
every elliptic rank-nine finite endpoint, and every connected affine diagram.
For each affine component it forms the **full** zero-neighborhood and checks
that it is componentwise parabolic of total Coxeter rank eight.  It then tests
the endpoint incidence of every elliptic edge.

All 13,182 elliptic edges have exactly two finite/ideal endpoints.  All 584
connected affine zero-neighborhood audits have rank eight.  There are zero
bad edges and zero bad affine neighborhoods.  This is precisely the
finite-covolume, rather than compact-only, Vinberg test.  The cusp counts
agree exactly with the independent double-description census.

## Seven genuinely different Gram types

No graph-isomorphism heuristic is required to keep the seven images apart.
The tuple

```text
(facet count, number of intersection-0 pairs,
 number of intersection-1 pairs, number of intersection-2 pairs)
```

is respectively

```text
I    (12,  53, 10,  3)
II   (12,  51, 15,  0)
III  (20, 138, 16, 36)
IV   (20, 140, 40, 10)
V    (20, 135, 24, 31)
VI   (20, 135, 45, 10)
VII  (20, 135, 30, 25).
```

These exact weighted-pair histograms are all different, so no facet
permutation can identify two normalized Gram matrices.  Shimada's records
20C and 20D, on the other hand, have literally equal configuration matrices
in the pinned ordering and give the single Type-VII image.

## Primary-source lineage

- Shigeyuki Kondō, *Enriques surfaces with finite automorphism groups*, Japan.
  J. Math. (N.S.) 12 (1986), 191--282,
  [doi:10.4099/math1924.12.191](https://doi.org/10.4099/math1924.12.191).
  The local open-access dissertation copy and page locators are recorded in
  [`data/external/kondo_enriques_1986/README.md`](../data/external/kondo_enriques_1986/README.md).
- Gebhard Martin, *Enriques surfaces with finite automorphism group in
  positive characteristic*, Algebraic Geometry 6 (2019), 592--649,
  [doi:10.14231/AG-2019-027](https://doi.org/10.14231/AG-2019-027).
  Its Main Theorem and Table 1 give the modern iff statement and moduli;
  local provenance is in
  [`data/external/martin_enriques_2019/README.md`](../data/external/martin_enriques_2019/README.md).
- Ichiro Shimada, *A note on configurations of (-2)-vectors on Enriques
  surfaces* and its official GAP supplement.  The exact type dictionary,
  files, and byte checksums are recorded in
  [`data/external/shimada_l10l26note/README.md`](../data/external/shimada_l10l26note/README.md).

## Replayable artifacts

- production verifier:
  [`scripts/certify_enriques_finite_aut_nef_chambers.py`](../scripts/certify_enriques_finite_aut_nef_chambers.py);
- exact matrices and summary:
  [`data/known_examples/enriques_finite_automorphism_nef_chambers.json`](../data/known_examples/enriques_finite_automorphism_nef_chambers.json);
- all 2,508 exact rays:
  [`data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz`](../data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz);
- regression tests:
  [`tests/test_enriques_finite_aut_nef_chambers.py`](../tests/test_enriques_finite_aut_nef_chambers.py).
- no-import raw-data hostile replay:
  [`results/enriques_finite_aut_nef_chambers_hostile_audit.md`](enriques_finite_aut_nef_chambers_hostile_audit.md), with its
  [exact certificate](../data/enumerations/enriques_finite_aut_nef_chambers_hostile_audit.json).

The full raw-to-cdd-to-Vinberg replay passes (`3 passed` in about 126 seconds
on the freezing machine).  A second certifier, whose sole chamber input is
Shimada's raw GAP file and which imports neither the production certifier nor
either production archive, also passes (`3 passed in 123.76s`).  It confirms
the precise `20C=20D` statement: the configuration matrices are literally
equal, while the two raw wall-coordinate arrays are not.  SHA-256 values are
byte-integrity pins, not a substitute for the exact mathematical checks:

```text
raw Shimada data  218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c
summary           245715b7352219a7ed1bfc34a43b5b37921e289d39829df51015c347b97072bc
rays              b72a2a5479a773c8151ed92cd8b698e71fe4b8cb0b4bf527245847df39cb46c5
hostile summary   224ded8c22c223cf379461f54650d2f7fd560292e1de79bbefe62e5d51d55b1b
```

## What this does and does not settle

This is a positive result: seven actual full-nef hyperbolic Coxeter polytopes
are obtained from an algebraic classification, with exact finite-volume
certificates and a converse carrier criterion for their seven Gram types.  It
shows concretely how a softer surface moduli problem can have fewer outputs
than surfaces and can classify the outputs completely.

Combined with Kondō's reflection fundamental-domain/Torelli bridge and the
Artebani--Hausen--Laface Cox-ring theorem, it now also proves that finite
hyperbolic volume of the complete Enriques nef chamber is equivalent to
finite generation of the surface's Cox ring.  The strengthened statement and
its exact scope are in
[`enriques_full_nef_finite_volume_classification.md`](enriques_full_nef_finite_volume_classification.md).

It does **not** prove the broad conjecture that every crystallographic Coxeter
polytope satisfying some presently known list of lattice conditions has a
surface carrier.  Nor does it lower the global compact `H9/14` boundary: every
Enriques image here is noncompact and has 12 or 20 facets.  The next useful
step is to abstract which parts of the Enriques proof survive beyond a finite
lookup list—most likely saturated root data, effectivity/extremality,
no-hidden-walls, and a period/moduli existence theorem—while the direct
higher-dimensional fishing continues independently.
