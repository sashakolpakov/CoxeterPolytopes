# A rational-involution `H10/37` chamber on Shimada's singular K3 surface `X_2`

## Result

There is an exact finite-volume noncompact hyperbolic Coxeter polytope in
\(\mathbb H^{10}\) with **37 facets** obtained from a genuine rational
involution of Shimada's singular complex K3 surface \(X_2\).  It has

- 3,236 finite vertices;
- 261 cusps; and
- 24,350 elliptic rank-9 edge diagrams.

Every one of those 24,350 edges has exactly two finite-or-ideal endpoints,
and every irreducible affine component has a full zero-neighbourhood of
parabolic rank 9.  Thus this is a literal application of Vinberg's
finite-volume criterion, not an inference from pairwise Coxeter angles.

The fixed divisor lattice is

\[
F=\operatorname{NS}(X_2)^\iota
  \cong U\oplus E_8(-2)\oplus\langle-2\rangle,
\]

the unique even 2-elementary lattice with
\((r,a,\delta)=(11,9,1)\).  Its determinant is 512 and its Smith form is
\(1^2 2^9\).  The 37 facet mirrors generate all of \(F\), and exhaustive
discriminant-module closure proves that \(F\) has no proper integral
overlattice preserved by all 37 reflections.

This is a newly extracted K3/divisor realization candidate and a new positive
entry in the maintained ledger.  It is **not** a new arithmetic lattice or a
new full reflection group.  Its reflection group has exact index

\[
[W_r(F):W_{37}]=241{,}920
\]

in the known 13-wall full reflection group.  Worldwide priority of the
literal 37-facet abstract polytope remains open.

## The complete algebraic denominator

The input is not a sampled wall list.  Shimada proves for \(X_2\) that:

1. \(D^{(0)}\) is an induced chamber contained in the nef cone;
2. its translates tessellate the nef cone;
3. its complete wall set consists of 478 walls in eight orbits of sizes
   \(36,12,40,90,30,30,120,120\); and
4. the ample stabilizer is
   \(\operatorname{Aut}(X_2,a_2)\cong\operatorname{PGL}_2(\mathbb F_9)\),
   of order 720.

All 478 published wall vectors and all 720 genuine surface-automorphism
matrices are frozen in the source stock.  Shimada's proposition on
involutions says that this stabilizer has 45 symplectic involutions and 36
rational involutions.

The exact fixed-rank histogram of all 720 matrices is

| fixed rank | 3 | 4 | 6 | 8 | 11 | 12 | 20 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| element count | 324 | 144 | 90 | 80 | 36 | 45 | 1 |

Every nonidentity element of fixed rank at least 11 is an involution.  The
high elements contain exactly 120 closed \(V_4\)'s:

- 90 have element-rank pattern \((11,11,12)\) and common fixed rank 7;
- 30 have pattern \((12,12,12)\) and common fixed rank 8.

This closes the proper-subgroup denominator above \(H^9\).  Indeed, if a
subgroup \(H\) had common fixed rank at least 11, every element of \(H\)
would be one of the 81 high involutions.  Two distinct elements would have a
high involutive product, hence commute and generate one of the enumerated
\(V_4\)'s—but every such \(V_4\) drops below rank 11.  Therefore the only
proper high possibilities are the two involution classes themselves.

## Restricting the complete chamber

For an involution \(\iota\), every one- or two-element orbit of walls is
restricted to \(\operatorname{Fix}(\iota)\) by its orbit sum.  Positive-ray
deduplication is performed only after all 478 walls are included.

For the rational class, the 478 walls form 242 involution orbits.  They give
207 distinct oriented restrictions:

| restricted square | negative | zero | positive |
|---|---:|---:|---:|
| count | 37 | 45 | 125 |

The nonnegative vectors cannot cut the interior of the selected component of
the hyperbolic positive cone.  Every one of the 37 negative vectors has an
exact strict timelike relative-interior witness against all other 206
restrictions, so none is redundant.  Their norms are
\(27\times(-2)+10\times(-4)\), and their Coxeter-pair census is

| relation | \(m=2\) | \(m=3\) | \(m=4\) | parallel | ultraparallel | forbidden |
|---|---:|---:|---:|---:|---:|---:|
| pairs | 405 | 15 | 70 | 121 | 55 | 0 |

For comparison, the symplectic class has fixed rank 12 and gives 90 strict
negative facets, but 352 of their pairs have forbidden normalized Gram
ratios.  It is decisively non-Coxeter.  Together with the \(V_4\) argument,
this proves that the rational class is the unique positive fixed-section
type above \(H^9\) in this complete 720-element denominator.

## Literal Vinberg finite-volume proof

The exact integer Coxeter Gram was passed to the same independent C++17
Vinberg enumerator used for the other high-dimensional certificates.  It
enumerates all finite special subsets and all irreducible affine components;
signed 128-bit Bareiss determinants are overflow-guarded.

The finite-special counts by rank are

```text
37, 490, 3415, 14490, 39842, 72352, 85585, 62400, 24350, 3236
```

Thus there are 24,350 elliptic rank-9 edges and 3,236 elliptic rank-10 finite
endpoints.  The endpoint types are

| finite endpoints | ideal endpoints | edge count |
|---:|---:|---:|
| 0 | 2 | 2,360 |
| 1 | 1 | 11,620 |
| 2 | 0 | 10,370 |

There are 813 irreducible affine components, whose full zero-neighbourhoods
give 261 distinct rank-9 cusps.  The bad-edge count and bad-parabolic-
neighbourhood count are both zero.  The polytope therefore has finite volume;
the 261 cusps prove that it is noncompact.

## Divisor lattice and reflection closure

The saturated Smith kernel of \(\iota^T-I\) has rank 11.  In its archived
integral basis, the fixed Gram has determinant 512, inertia \((1,10,0)\), and
Smith invariants \(1^2 2^9\).  Its discriminant quadratic form has nonintegral
values, so \(\delta=1\).  Nikulin's uniqueness theorem for indefinite even
2-elementary lattices gives

\[
F\cong U+E_8(-2)+\langle-2\rangle.
\]

The primitive mirror types are

- 21 roots of type \((-2,\operatorname{div}=2)\);
- 6 roots of type \((-2,\operatorname{div}=1)\); and
- 10 roots of type \((-4,\operatorname{div}=2)\).

All 37 reflection matrices are integral and unimodular on \(F\).  The
37-by-11 coordinate matrix has Smith diagonal \(1^{11}\), proving that the
mirrors span the full saturated divisor lattice.

The discriminant group has order 512.  Exhausting every self-integral
discriminant element, its orbit under all 37 reflections, every additive
orbit closure, and all sums of invariant closures leaves stable integral
subgroup orders `[1]` only.  Hence there is no hidden proper integral
reflection-stable overclosure.  The scaled dual \(F^\vee(2)\), which is not an
overlattice, is the odd determinant-four class in Scharlau--Kirschmer
`List10`, entry 6.

## Exact full-group lineage

Alexeev--Engel's finite-covolume theorem and Coxeter figure give the
13-wall full root-reflection chamber \(W_r\) for the same lattice
\((11,9,1)\).  The already frozen geometric `K3-S3-C2-H10-13` chamber has the
exact relation graph of that figure, 34 finite vertices, and four cusps.  In
Scharlau--Kirschmer notation this is `List10` entry 6, genus
\(I_{10,1}(2^2)\), after passing to the scaled dual/sign convention.

Solomon's formula gives

\[
\chi(W_{37})=-\frac{527}{138240},\qquad
\chi(W_r)=-\frac{527}{33443020800}.
\]

Every displayed reflection is a root reflection of \(F\), so
\(W_{37}\subset W_r(F)\).  Both groups have finite covolume and their
nonzero Coxeter orbifold characteristics therefore give

\[
[W_r(F):W_{37}]
=\frac{\chi(W_{37})}{\chi(W_r)}
=241920.
\]

This settles the arithmetic lineage completely: no new lattice class and no
new maximal group are involved.  What may be new is the literal 37-wall
finite-index chamber and, more concretely, its extraction from a genuine
rational involution on the singular K3 surface \(X_2\).

Before insertion into the ledger, no `H10/37` header occurred in the frozen
project ledger.  A strengthened [bounded lineage audit](shimada_aut3singk3_x2_h10_37_global_lineage.md)
scans every 37-node record in the pinned CoxIter ordinary corpus (1,214
files) and index-two doubling archive (2,959 files), regardless of declared
dimension, at commit `b800be48d8240fbabd2d1c9542e61054584ad7fe`.
All 2 ordinary and 43 doubling records are exact nonmatches by the complete
relation-label histogram.  This is still a bounded comparison only; it does
not prove worldwide priority.

## Algebraic meaning

This construction does use the surface, rather than merely borrowing a
Lorentz lattice.  The involution is one of Shimada's actual rational surface
automorphisms.  Over \(\mathbb Q\), invariant Néron--Severi classes are the
pullbacks of divisor classes on the normal rational quotient
\(X_2/\langle\iota\rangle\).  The polytope is the intersection of Shimada's
complete induced nef chamber with this invariant divisor space.  It is a
finite-index subchamber of the full root chamber, not a claim that the 37
walls form the entire nef cone of the quotient.

This is therefore the requested rigid algebro-geometric mechanism in a
successful instance:

\[
(X_2,\iota,D^{(0)})
\longmapsto
\bigl(\operatorname{NS}(X_2)^\iota,
      D^{(0)}\cap\operatorname{Fix}(\iota)\bigr)
\longmapsto
P^{10}_{37}.
\]

## Reproducibility

- [complete 478-wall/720-automorphism stock](../data/enumerations/shimada_aut3singk3_x2_478_wall_stock.json)
- [composite `H10/37` certificate](../data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json)
- [bounded global-lineage certificate](../data/enumerations/shimada_aut3singk3_x2_h10_37_global_lineage.json)
- [stock ingester](../scripts/ingest_shimada_aut3singk3_x2_478_stock.py)
- [`H10/37` certifier](../scripts/certify_shimada_aut3singk3_x2_rational_involution_h10_37.py)
- [bounded global-lineage certifier](../scripts/certify_shimada_aut3singk3_x2_h10_37_global_lineage.py)
- [stock regression test](../tests/test_shimada_aut3singk3_x2_478_stock.py)
- [`H10/37` regression test](../tests/test_shimada_aut3singk3_x2_rational_involution_h10_37.py)
- [bounded global-lineage regression test](../tests/test_shimada_aut3singk3_x2_h10_37_global_lineage.py)

Replay with

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_shimada_aut3singk3_x2_rational_involution_h10_37 \
  --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_shimada_aut3singk3_x2_478_stock.py \
  tests/test_shimada_aut3singk3_x2_rational_involution_h10_37.py
```

Frozen SHA-256 values at certification time:

```text
complete stock JSON  10d1c767d21d7194c6ccb603c5dff29520201954b4a07fd3c5f214d452ea2f74
H10/37 certifier    689c27f35f8be21cb5b7642211a77ba0b47fb777f0b077b63acafc7ca9b579df
H10/37 certificate  12914344403b7d3d6b69a2866c4150f0a0420a64914e4bc8584fa372799599a7
Vinberg C++ source  45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045
```
