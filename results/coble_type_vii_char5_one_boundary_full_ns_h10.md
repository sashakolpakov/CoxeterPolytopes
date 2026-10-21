# The characteristic-five Type-VII Coble full-nef chamber

**Status (2026-08-09): proved exact positive `SURFACE -> POLYTOPE`
theorem.**  This result has inverse-carrier score **N/A/L0**: it starts from
Kondō's surface and makes no `POLYTOPE -> SURFACE` implication.

## The theorem

Let `S` be the unique strong Coble surface of Type VII over an algebraically
closed field of characteristic `5`.  Then the projectivized full
Néron--Severi nef cone of `S` is a finite-volume, noncompact arithmetic
Coxeter polytope

```text
COBLE-VII-p5-b1-H10-72  in H10
```

with the following exact data.

| facets | primitive wall norms `-1,-2` | finite vertices | ideal vertices | hyperideal rays | edges |
|---:|---:|---:|---:|---:|---:|
| 72 | `52,20` | 18,761 | 1,731 | 0 | 143,565 |

All 72 walls are irredundant.  Their reflections are integral involutive
isometries, all wall pairs have Coxeter angles or are parallel/ultraparallel,
and their primitive normals span the whole lattice

\[
              \operatorname{NS}(S)\cong I_{1,10}.
\]

This supplies a third exact one-boundary Coble full-nef chamber in `H10`,
after Types I and II, and the first such chamber in positive characteristic
completed by this project.

The reproducible artifacts are:

- [summary certificate](../data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json);
- [complete exact ray/incidence archive](../data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10_rays.jsonl.gz);
- [complete elliptic-edge/endpoint archive](../data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10_edges.jsonl.gz);
- [certifier](../scripts/certify_coble_type_vii_char5_one_boundary_full_ns_h10.py) and its [independent finite-diagram enumerator](../scripts/certify_coble_type_vii_char5_one_boundary_diagrams.cpp);
- [tests](../tests/test_coble_type_vii_char5_one_boundary_full_ns_h10.py).

## Why the 72 walls are the complete nef wall set

Kondō's Theorem 1.1 and Table 2 give a unique Type-VII Coble surface in
characteristic `5`, with one boundary component, twenty effective
irreducible roots, and automorphism group `S5`.  His Example 4.15 constructs
the surface and says that its twenty curves have the same graph
`Gamma_VII` as the Type-VII Enriques configuration; Proposition 4.16 records
the invariant `(A9+A1,Z/2)`.

Write

\[
             B\in|-2K_S|,\qquad b=B/2=-K_S.
\]

Since `B` is a rational `(-4)`-curve, `b^2=-1`.  Kondō's Proposition 3.1
gives

\[
 \operatorname{CM}(S)=K_S^\perp\cong E_{10},\qquad
 \operatorname{NS}(S)=E_{10}\mathbin\perp\langle-1\rangle
                     \cong I_{1,10}.
\]

The integral transport to Shimada's exact `20C` embedding follows from a
small intrinsic glue lemma, not merely from recognizing the graph.  The
abstract lattice `L` spanned by the twenty roots has

```text
det(L)=-4,       A_L=(Z/2)^2.
```

In an exact HNF basis, the three nonzero discriminant classes have quadratic
values `-5/2`, `-2`, and `-7/2`.  Exactly one is even/isotropic modulo `2Z`.
Consequently `L` has a unique even index-two overlattice; it is unimodular of
signature `(1,9)`, hence `E10`.  The entrywise `Gamma_VII` Gram isometry
therefore extends across this unique overclosure to an `E10` isometry.  This
proves the existence of a marking carrying the Coble configuration and its
cusp cone to `20C`; no uniqueness of the marking itself is asserted.  The
executable freezes all four discriminant representatives, their q-values,
the unique glue vector, and the resulting unimodular overbasis.  It also
checks the Smith diagonal `1,1,1,1,1,1,1,1,1,2` of the root span in the
stored `E10` marking.  This is more than an abstract unweighted-graph match.

Kondō's Proposition 3.7 makes those twenty roots the **complete** effective
irreducible-root list.  With one boundary, every such root is the class of a
`(-2)`-curve.

The remaining negative curves are recovered from the cusps, not guessed.
For every primitive isotropic ray `f` of the complete Coble--Mukai chamber,
Kondō's Lemma 3.11 gives a genus-one fibration.  His Halphen description in
Proposition 3.8 says that a one-boundary surface comes from blowing up the
node or cusp of one `I1` or `II` fibre on an index-two Halphen surface.  Its
unique exceptional curve `E` satisfies

\[
                         B+2E=2f,
                 \qquad E=f-b.                 \tag{1}
\]

Conversely every `(-1)`-curve satisfies, by adjunction,

\[
 f=E+b\in b^\perp,\qquad f^2=0,
\]

and `f.R=E.R>=0` for every effective CM root `R`; hence `f` is a primitive
nef CM cusp.  Thus (1) is a bijection

```text
all primitive CM cusp rays  <-->  all (-1)-curves.
```

Exact double description of the `20C` CM chamber gives 51 cusps.  The full
primitive nef-wall list is therefore

\[
 \{b\}\ \cup\ \{R_1,\ldots,R_{20}\}\ \cup
 \ \{f-b:f\text{ is one of the 51 CM cusps}\},               \tag{2}
\]

so it has `1+20+51=72` members.  Physically these are one boundary
`(-4)`-curve, twenty `(-2)`-curves, and fifty-one `(-1)`-curves; the boundary
hyperplane uses its primitive norm-`-1` normal `b`.  Kondō's negative-curve
trichotomy says that boundary, `(-2)`, and `(-1)` curves are the only
negative curves on a Coble surface.  Equation (2) exhausts all three kinds,
so this is the complete full nef cone rather than a selected reflection
subcone.

## Coxeter and arithmetic checks

In the exact `E10+<-1>` marking the complete primitive normal set has

```text
norm -1: 52       norm -2: 20       wall rank: 11
Smith span: 1^11  (the walls generate all of I_(1,10))
```

The pair profile is

```text
(-1,-1): pairing 0/1/2/3 = 405/561/330/30
(-2,-1): pairing 0/1/2/3 = 520/300/180/40
(-2,-2): pairing 0/1/2   = 135/30/25.
```

Consequently the angle census is

```text
m2=1,060; m3=30; m4=300; parallel=586; ultraparallel=580;
non-Coxeter=0.
```

Every one of the 72 reflection matrices is checked entrywise to be integral,
involutive, and lattice-preserving.  Once finite covolume is established
below, the group is therefore an arithmetic reflection subgroup of the full
`I_(1,10)` reflection group.

## Literal Vinberg finite-volume certificate

The finite-volume conclusion does **not** come from saying that the chamber
is noncompact, nor merely from finding no negative-square ray in one partial
enumeration.

First, `cddexec_gmp --repall` performs complete rational double description
and returns

```text
20,492 extreme rays = 18,761 finite + 1,731 ideal + 0 hyperideal.
```

Every facet's incident rays span a rank-ten face.  All finite vertices are
simple, while the ideal vertices have the following complete incidence
census:

| incident walls | finite vertices | ideal vertices |
|---:|---:|---:|
| 10 | 18,761 | 120 |
| 11 | 0 | 500 |
| 12 | 0 | 721 |
| 13 | 0 | 295 |
| 14 | 0 | 90 |
| 15 | 0 | 5 |

The sum of all exact oriented rays gives a strict positive-square interior
point with strictly positive slack at every wall.  The certifier also checks
its exact Lorentz pairing with **every** causal extreme ray and requires one
strict common sign.  This rules out accidentally mixing the two components
of the positive cone.

Second, the finite Coxeter subdiagram pass is exhaustive, including every
disconnected diagram.  A deterministic compiled classifier screens every
pair-compatible clique against the complete finite crystallographic Coxeter
list:

| order | compatible cliques screened | elliptic diagrams |
|---:|---:|---:|
| 9 | 15,619,250 | 143,565 |
| 10 | 24,466,316 | 18,761 |

The 143,565 elliptic rank-nine masks equal, set-for-set, the supports of the
143,565 exact geometric edges from cdd adjacency.  The 18,761 elliptic
rank-ten masks equal, set-for-set, the finite endpoints.  Containment in the
complete finite and ideal active sets gives exactly two endpoints for every
edge:

```text
FF=60,250; FI=67,110; II=16,205; bad edges=0.
```

Third, the parabolic test also includes disconnected cusps.  Starting from
all 4,842 connected finite masks, the verifier attaches every compatible
nonorthogonal node through order 12 and uses an exact determinant-zero test
with a positive-definite principal deletion.  It finds all 4,943 connected
affine components:

```text
order 2/3/4/5/6/7/8/9/10
count 586/1330/1590/672/250/240/135/20/120.
```

For each component `J`, its complete zero-neighbourhood `N(J)` is exactly
one of the 1,731 certified ideal vertices.  Conversely every ideal active
set decomposes into affine connected components of total Coxeter rank nine,
and every component has that whole ideal set as its full zero-neighbourhood.
The numbers of connected affine components per cusp are

```text
1/2/3/4/5/6 components: 120/500/721/295/90/5 cusps.
```

Thus Vinberg's rank-nine edge-extension criterion and his full parabolic
zero-neighbourhood criterion both pass literally, with zero failures.  The
ideal vertices prove noncompactness; the two Vinberg tests prove finite
volume.

## Euler characteristic, index, and lineage

The independent connected-elliptic/Solomon calculation gives

\[
             \chi(W_{\rm VII,full})=-\frac{527}{7560}.
\]

The even-dimensional hyperbolic Gauss--Bonnet formula consequently gives the
exact chamber volume

\[
 \operatorname{vol}(P_{72})
   =-\frac{32\pi^5}{945}\,\chi(W_{\rm VII,full})
   =\frac{2108\pi^5}{893025}.
\]

Comparison with `chi(W(P10))=-1/11147673600` gives index

\[
                    [W(P10):W_{\rm VII,full}]=777093120.
\]

This integer is **not new lineage evidence**.  It equals Shimada's published
index for the known Type-VII `20C` induced chamber in `H9`; the analogous
index equality already occurs for the Type-I and Type-II one-boundary
full-NS lifts.  The equality is instead a striking structural invariant of
the cusp-extension construction and deserves a separate general theorem.

Nor is the bare number `51` a novelty claim.  Mukai's 2018 Oberwolfach
report already records 51 genus-one fibrations, with orbit/type counts
`A8:20`, `A7+A1:15`, `A4+A4:6`, and `A5+A2+A1:10`, for a characteristic-two
Coble surface `R_{-10}` with automorphism group `S5`.  That is a different
carrier—ten boundary components and five `(-2)`-curves—not the present
characteristic-five one-boundary Type-VII surface.  Nevertheless this is an
**exact abstract seed identification**, not a numerical resemblance:
decomposing all 51 archived `20C` ideal active sets into nonorthogonal affine
components gives

```text
A8:20,  A7+A1:15,  A4+A4:6,  A5+A2+A1:10,
```

entry-for-entry the same type census as Mukai's report.  Thus the Type-VII
`H9` cusp census was already present in that different characteristic-two
geometry.  The new project content here is
the **bijection of those 51 cusps with all 51 `(-1)`-curves on this
one-boundary surface**, and the resulting complete `H10/72` full-nef cone.

The status is therefore deliberately split:

- the Type-VII surface, its twenty-root configuration, and the `20C=20D`
  `H9` Coble--Mukai/Enriques chamber are known;
- the 51 cusp lifts, complete 72-wall full-NS Gram, exact face lattice, and
  literal Vinberg certificate are new in this project ledger;
- worldwide priority for the literal abstract `H10/72` chamber remains
  open, so no worldwide-new-polytope claim is made.

As a bounded check, the pinned CoxIter corpus at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe` contains no `H10/72` header among
1,214 ordinary records or 2,959 archived doubling records.  This is only a
bounded exclusion.

## Positive surface-list completion

The subsequent positive search closed Type VI in characteristic five, Type
MI in characteristic three, and finally the characteristic-three Type-V
two-boundary row.  Type V has a unique full-`NS` glue; seven of its 29 CM
cusps split over both boundaries and give fourteen exceptional curves.  The
result is the complete finite-volume noncompact `H11/36` chamber with
`6,356+432+0` exact rays and 47,929 Vinberg-certified edges.  Thus all nine
finite-automorphism Coble templates now have exact full-Néron--Severi
polytope certificates.  See the [`Type-MI theorem`](coble_type_mi_full_ns_h11.md)
and [`Type-V theorem`](coble_type_v_two_boundary_full_ns_h11.md).

## Sources

- Shigeyuki Kondō, [*Coble surfaces with finite automorphism group*](https://arxiv.org/abs/2102.00591), especially Theorem 1.1/Table 2, Propositions 3.1, 3.4, 3.7, 3.8 and 4.16, Lemma 3.11, and Example 4.15.
- E. B. Vinberg, [*Hyperbolic reflection groups*](https://www.mathnet.ru/eng/rm2140), Russian Mathematical Surveys 40:1 (1985), especially the finite-volume/parabolic-neighbourhood criterion.
- Shigeru Mukai, [*Automorphism of K3 surfaces and decomposition groups of rational sextics*](https://ems.press/content/serial-article-files/46758?nt=1), Oberwolfach Reports 33/2018, pp. 2024--2027, for the distinct characteristic-two 51-fibration census.
