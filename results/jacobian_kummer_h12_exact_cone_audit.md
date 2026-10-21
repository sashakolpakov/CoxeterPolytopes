# A Jacobian--Kummer Coxeter polytope in \(\mathbb H^{12}\)

## Outcome

The trace-9 involution class in the full lattice/chamber stabilizer

\[
O(S_{16},D_{16})\cong (C_2)^5\rtimes S_6
\]

produces a **172-facet, finite-volume, noncompact hyperbolic Coxeter
polytope in \(\mathbb H^{12}\)**.  The construction starts with the
Néron--Severi lattice and the complete 316-wall chamber \(D_{16}\) of a
general Jacobian Kummer K3 surface, then intersects that chamber with the
fixed Lorentzian space of the involution.

The involution belongs to the full published **lattice/chamber** stabilizer.
It is not assumed to come from an automorphism of the K3 surface.  Thus the
construction uses exactly the more permissive algebro-geometric idea: the
surface supplies the hyperbolic divisor lattice and a complete chamber, while
lattice symmetries supply candidate sections.

Finite volume has two exact, logically separate certificates:

1. a pure Coxeter-diagram application of Vinberg's edge-endpoint and full
   parabolic-neighborhood criterion; and
2. a GMP rational double-description of the complete chamber cone, with no
   spacelike extreme ray.

This is a genuinely positive example produced by the search.  A bounded
literature/corpus check found no prior occurrence of this 172-facet chamber,
but global priority is not claimed.

## Construction and exact Coxeter gate

Dolgachev--Shimada's published lattice
\(S_{16}=\operatorname{NS}(Y_{16})\) has rank 17, determinant 64, and
signature \((1,16)\).  The selected involution is the unique involution class
in \(O(S_{16},D_{16})\) with lattice trace 9.  Its fixed space has rank 13 and
signature \((1,12)\), hence projectivizes to \(\mathbb H^{12}\).

The involution has 180 orbits on the 316 published walls: 44 singleton
orbits and 136 pairs.  Restriction to the fixed space turns these into 180
homogeneous inequalities.  Exact conic reduction gives:

- 172 irredundant facets, each with a strict timelike relative-interior
  witness;
- 8 redundant inequalities, each with an explicit nonnegative rational
  dependence; and
- orbit-sum normal squares \((-1)^{124},(-2)^{40},(-4)^8\).

All \(\binom{172}{2}=14{,}706\) pairs pass the exact Coxeter gate:

| relation | count |
|---|---:|
| orthogonal, \(m=2\) | 2,820 |
| \(m=3\) | 144 |
| \(m=4\) | 800 |
| parallel | 1,750 |
| ultraparallel | 9,192 |

There is no forbidden angle.  The 172 inequalities are the complete section
of the complete \(D_{16}\) chamber, not a wall prefix.

## The fixed divisor lattice

The saturated integral fixed lattice \(L=S_{16}^{\langle g\rangle}\) has

\[
\operatorname{rk}L=13,\qquad
\operatorname{sig}L=(1,12),\qquad
|\det L|=256,
\]

and Smith invariants

\[
1^8,\;2^2,\;4^3.
\]

An explicit discriminant-group basis gives

\[
q_L\cong u(2)\oplus\langle 3/4\rangle
              \oplus\langle 3/4\rangle
              \oplus\langle 7/4\rangle.
\]

The standard model

\[
L\cong U(2)\oplus D_5(-1)^2\oplus\langle-4\rangle
\]

has the same signature and discriminant form.  Nikulin's indefinite-lattice
uniqueness criterion applies because \(13\geq l(A_L)+2=7\), giving the stated
isometry class.

After converting the 172 hyperplanes to primitive vectors of \(L\), their
types are:

| primitive square and divisibility | facets |
|---|---:|
| \((-2,1)\) | 16 |
| \((-4,2)\) | 36 |
| \((-4,4)\) | 96 |
| \((-8,4)\) | 24 |

Every corresponding reflection is integral on \(L\).  This identifies the
example intrinsically in divisor-lattice terms, independently of its
17-dimensional parent chamber.

## Literal Vinberg finite-volume proof

The diagram verifier first enumerates every connected elliptic component by
exact connected growth.  There are 16,908.  Every elliptic diagram is then
formed, without duplication, as an orthogonal union of these uniquely
determined connected components.  This gives

- 113,274 elliptic rank-12 finite vertices; and
- 910,316 elliptic rank-11 edge diagrams.

Separately, the verifier enumerates every connected affine component.  Their
order histogram is

\[
2^{1750},3^{2784},4^{3420},5^{3768},6^{1760},7^{544},
8^{1128},9^{1536},10^{384},
\]

for 17,074 components in total.  Taking the **full zero-neighborhood** of
each component gives exactly 6,214 distinct ideal vertices.  Every
neighborhood is componentwise affine of total Coxeter rank 11; the number of
invalid neighborhoods is zero.  Their active-facet histogram is

\[
13^{2304},\quad14^{3256},\quad15^{576},\quad16^{76},\quad18^2.
\]

The endpoint audit is exhaustive:

| endpoints of an elliptic rank-11 edge | edges |
|---|---:|
| finite--finite | 505,568 |
| finite--ideal | 348,152 |
| ideal--ideal | 56,596 |
| **total** | **910,316** |

Every edge has exactly two endpoints.  There are no bad edges and no
endpoint-generated nonedges.  The incidence checks are

\[
113274\cdot12=1{,}359{,}288=2(505568)+348152
\]

and

\[
461{,}344=348152+2(56596).
\]

Vinberg's criterion therefore proves **finite volume**.  The 6,214
parabolic vertices prove **noncompactness**.

## Independent exact cone certificate

In 13 independent normal coordinates, the complete chamber is a
full-dimensional pointed rational cone with 172 integral inequalities.  The
published ample class is a strict timelike interior point.  The GMP rational
build of `cddlib` gives all primitive extreme rays:

| extreme-ray kind | count |
|---|---:|
| timelike, hence a finite vertex | 113,274 |
| lightlike, hence an ideal vertex | 6,214 |
| spacelike/hyperideal | **0** |
| **total** | **119,488** |

Every inequality, active set, and Lorentz square was replayed with exact
integer arithmetic.  All finite vertices have exactly 12 active facets.  The
6,214 lightlike active sets agree exactly with the exhaustive affine
neighborhoods above.

Since the exact double-description lists every ray of the pointed cone, zero
spacelike rays rules out a hyperideal escape.  This supplies a direct
finite-volume/noncompact certificate independent of the pure-diagram
enumeration.

## Novelty audit

The following is evidence of absence, not a proof of global novelty.

- Dolgachev--Shimada publish \(D_{16}\), its 316 walls, and its stabilizer,
  but do not discuss this trace-9 fixed chamber or a 172-wall fold.
- A check of the official CoxIter example corpus at commit
  `b800be48d8240fbabd2d1c9542e61054584ad7fe` found 1,214 `.coxiter` files,
  none with more than 50 facets.  Its \(\mathbb H^{12}\) examples have only
  13, 14, 16, or 18 facets.
- Searches for the fixed-lattice isometry class, discriminant form, and the
  172-facet data found no match.

Accordingly, the safe status is: **newly constructed and exactly certified
here; apparently absent from the checked literature and corpus; global
priority unresolved**.

## Reproduction and sources

- Metric/facet certificate:
  [`jacobian_kummer_h12_interim_facet_gram.json`](../data/enumerations/jacobian_kummer_h12_interim_facet_gram.json)
- Complete affine-neighborhood audit:
  [`jacobian_kummer_h12_affine_neighborhood_audit.json`](../data/enumerations/jacobian_kummer_h12_affine_neighborhood_audit.json)
- Pure-diagram Vinberg certificate:
  [`jacobian_kummer_h12_vinberg_volume_certificate.json`](../data/enumerations/jacobian_kummer_h12_vinberg_volume_certificate.json)
- Independent exact cone certificate:
  [`jacobian_kummer_h12_exact_cone_certificate.json`](../data/enumerations/jacobian_kummer_h12_exact_cone_certificate.json)
- Primitive cone rays and active facets:
  [`jacobian_kummer_h12_exact_cone_rays.jsonl.gz`](../data/enumerations/jacobian_kummer_h12_exact_cone_rays.jsonl.gz)
- Pure-diagram verifier:
  [`verify_jacobian_kummer_h12_vinberg.py`](../scripts/verify_jacobian_kummer_h12_vinberg.py)
- Independent GMP verifier:
  [`independent_certify_jacobian_kummer_h12.py`](../scripts/independent_certify_jacobian_kummer_h12.py)
- Fixed-lattice and integral-reflection certificate:
  [`certify_jacobian_kummer_h12_fixed_lattice.py`](../scripts/certify_jacobian_kummer_h12_fixed_lattice.py),
  [`jacobian_kummer_h12_fixed_lattice_certificate.json`](../data/enumerations/jacobian_kummer_h12_fixed_lattice_certificate.json)
- Bounded primary-source and CoxIter-corpus comparison:
  [`jacobian_kummer_h12_literature_audit.md`](jacobian_kummer_h12_literature_audit.md)
- Regression tests:
  [`test_jacobian_kummer_h12_vinberg.py`](../tests/test_jacobian_kummer_h12_vinberg.py) and
  [`test_jacobian_kummer_h12_exact_cone.py`](../tests/test_jacobian_kummer_h12_exact_cone.py),
  [`test_jacobian_kummer_h12_fixed_lattice.py`](../tests/test_jacobian_kummer_h12_fixed_lattice.py)

Heavy exact replay commands are:

```bash
PYTHONPATH=src:. python3 scripts/verify_jacobian_kummer_h12_vinberg.py --check-archive
PYTHONPATH=src:. python3 scripts/independent_certify_jacobian_kummer_h12.py
```

Primary sources:

- I. Dolgachev and I. Shimada, *15-nodal quartic surfaces. Part II: the
  automorphism group*, [DOI 10.1007/s12215-019-00464-7](https://doi.org/10.1007/s12215-019-00464-7),
  [arXiv:1908.05390](https://arxiv.org/abs/1908.05390).
- V. V. Nikulin, *Integral symmetric bilinear forms and some of their
  applications*, [DOI 10.1070/IM1980v014n01ABEH001060](https://doi.org/10.1070/IM1980v014n01ABEH001060).
- R. Guglielmetti, *CoxIter -- Computing invariants of hyperbolic Coxeter
  groups*, [DOI 10.1112/S1461157015000273](https://doi.org/10.1112/S1461157015000273).
