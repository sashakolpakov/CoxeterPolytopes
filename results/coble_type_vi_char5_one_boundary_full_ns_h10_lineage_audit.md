# Hostile lineage audit: Type-VI characteristic-five full-nef chamber

Status date: 2026-08-09.

## Verdict

The construction is a certified positive hyperbolic example:

```text
COBLE-VI-p5-b1-H10-78
dimension                         10
facets                            78
finite vertices                   21,661
ideal vertices                    1,797
hyperideal extreme rays           0
edges                             166,320
volume status                     finite, noncompact
ambient lattice                   I_(1,10)
index in W(P10)                   906,608,640
volume                            1054*pi^5/382725.
```

The exact novelty claim is narrower than “new arithmetic Coxeter group.”
Kondō's Type-VI surface is known, Shimada's twenty-wall `20E` chamber in
`H9` is known, and the ambient odd unimodular arithmetic group is
classical.  The project-new positive content is the complete one-boundary
lift of all 57 CM cusps to all 57 `(-1)`-curves, the resulting full
78-wall Néron--Severi chamber, and its exact face/Vinberg certificate.
Worldwide priority for the literal abstract subgroup chamber remains open,
so it is not claimed.

The direction is strictly `SURFACE -> POLYTOPE`, inverse score `N/A/L0`.
Nothing in the argument starts from an independently defined raw-polytope
class or forces a surface carrier for every member of such a class.

This verdict applies to the argument audited here.  A later independent
companion starts from the raw 78 walls, reconstructs the Desargues/cut plane
design and its integral sextic, and returns the same chamber; it is recorded
separately as structural P-first `L3`, not retroactively merged into this
surface-first proof.  See
[`coble_type_vi_char5_desargues_inverse.md`](coble_type_vi_char5_desargues_inverse.md).

## Primary-source payload

Kondō's classification supplies:

- the unique characteristic-five, one-boundary Type-VI isomorphism class;
- `CM(S)=E10` and the nef-slice interpretation;
- the negative-curve trichotomy;
- the finite-index argument exhausting effective irreducible roots;
- the Halphen description and primitive cusp/fibration correspondence;
- the twenty `Gamma_VI` curves, `Aut(S)=S5`, and
  `R=(E6+A4,{0})`.

It does not print “57 `(-1)`-curves,” “78 full-nef walls,” the full `H10`
Gram, or the vertex/edge/cusp census.  Those statements are derived and
certified here.

## Integral lineage bridge

The imported CM data are Shimada's exact `20E` vectors, but the transport is
intrinsic.  The labelled `Gamma_VI` and `20E` root Gram matrices agree
entrywise.  The twenty roots have coordinate-span Smith diagonal `1^10`,
an HNF basis of determinant one, and unimodular Gram determinant `-1`.
They therefore span `E10` itself.  The root-Gram isometry extends to an
integral `E10` marking without a discriminant-form or glue choice.

This closes the gap between an unweighted graph identification and a valid
transport of all 57 exact ideal rays.

## Completeness of the 78 walls

For the unique boundary `B`, put `b=B/2=-K`.  Then

\[
 NS(S)=E10\perp\langle-1\rangle,qquad b^2=-1.
\]

The twenty roots are all `(-2)`-curves.  Every primitive CM cusp `f` gives a
one-boundary index-two Halphen fibration and its unique exceptional curve

\[
                         E=f-b,qquad B+2E=2f.
\]

Conversely every `(-1)`-curve returns the primitive nef cusp `f=E+b`.
Thus the 57 exact cusps of `20E` give all 57 `(-1)`-curves.  Together with
the boundary, Kondō's negative-curve trichotomy proves the complete stock

```text
1 boundary + 20 (-2)-curves + 57 (-1)-curves = 78 walls.
```

Independent cdd irredundancy checks show that every one is a genuine facet.

## Hostile finite-volume audit

GMP rational double description enumerates all 23,458 extreme rays:

```text
21,661 finite; 1,797 ideal; 0 hyperideal.
```

All finite vertices have ten active facets.  The ideal incidence histogram
is `11:480, 12:772, 13:460, 14:75, 15:10`.  A strict timelike interior
vector pairs with every causal ray with one sign; the exact minimum is
`362134>0`.

The independent exhaustive finite-diagram pass screens

```text
26,769,085 rank-9 masks  -> 166,320 elliptic edge diagrams;
40,619,993 rank-10 masks ->  21,661 finite endpoints.
```

These sets agree with cdd edge and finite-vertex supports set-for-set.  Every
edge has exactly two endpoints:

```text
FF 69,210; FI 78,190; II 18,920; failures 0.
```

The affine pass finds 5,551 connected components, and this example is a
particularly strong disconnected-cusp test: **none** of the 1,797 maximal
cusps is connected.  They contain two through six affine components, with
histogram

```text
2:480, 3:772, 4:460, 5:75, 6:10.
```

Every component's complete zero-neighbourhood is a rank-nine parabolic
cusp, and every ideal active set decomposes componentwise to total rank
nine.  There are zero bad affine neighbourhoods.  Hence finite volume is
proved, rather than inferred from the absence of hyperideal rays alone.

## Arithmetic and catalogue lineage

All walls span `I_(1,10)` and all reflections are integral.  Solomon's
formula gives

\[
 \chi(W_{78})=-527/6480,
 \quad [W(P10):W_{78}]=906608640,
 \quad \operatorname{vol}(P_{78})=1054\pi^5/382725.
\]

The index is already encoded by Shimada's known `20E` row: its
`volindex=2^7*3^4*5=51840`, and division of the full mod-two orthogonal
order `46998591897600` gives `906608640`.  Thus neither the arithmetic
class nor this index is new.

The pinned CoxIter scan at commit
`b800be48d8240fbabd2d1c9542e61054584ad7fe` covers 1,214 ordinary records
and 2,959 doubling records, with no `H10/78` header.  Exact searches for the
facet count, Euler characteristic, and index find the known Type-VI/`20E`
lineage but no source printing this full chamber.  This is bounded B1-style
discovery evidence only.  CoxIter is not a worldwide census, and worldwide
literal priority stays unresolved.

The defensible novelty boundary is therefore:

> A new exact full-divisor realization and positive chamber in this project;
> not a new dimension, lattice, maximal reflection group, arithmetic class,
> or proved worldwide-new abstract polytope; and no inverse-carrier theorem.

## Reproducibility

- [main theorem](coble_type_vi_char5_one_boundary_full_ns_h10.md)
- [summary certificate](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json)
- [complete ray archive](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_rays.jsonl.gz)
- [complete edge archive](../data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10_edges.jsonl.gz)
- [certifier](../scripts/certify_coble_type_vi_char5_one_boundary_full_ns_h10.py)
- [regression and full-replay tests](../tests/test_coble_type_vi_char5_one_boundary_full_ns_h10.py)

The three archive/construction tests pass independently of a fresh cdd
replay.  A complete live replay passed on 2026-08-09 in 255.66 seconds; it
regenerated the rational cone, every finite mask, every affine
zero-neighbourhood, and the archives from the current certifier.  A separate
root-agent check-only replay then exited zero against the frozen bytes.

Frozen artifact SHA-256 values are

```text
summary  2f58973a6e7305032588681cdec0f2186d9ce640207a72eeecadb6397a0bfffc
rays     72ec876fcc319d5366f90303fed2ce7c6d1b1d6f37d3f2a601898ab91de1c438
edges    d841d5d747041824be1d1b3c839cb64f5665130e75bd7cc94f07e5d7b342d6db
```
