# The characteristic-two involution divisor section: H10/47

Status: exact finite-volume noncompact Coxeter polytope.  Its ambient
arithmetic reflection group is classical; the explicit K3 fixed-divisor
realization is priority-open.

## Construction

Start with the complete 210-wall chamber in the Neron--Severi lattice

\[
  U\oplus D_{20}(-1)
\]

of the supersingular characteristic-two K3 surface of Artin invariant one.
The ambient wall Gram has rank 22 and inertia \((1,21,188)\).  We take an
involution in ATLAS wall-diagram class `2d` (class size 1008), which lies in
the index-two Dolgachev--Kondo geometric subgroup.

Its 210 walls split into six fixed walls and 102 pairs.  Restriction to the
fixed divisor space produces 108 orbit normals: 47 negative, 15 isotropic,
and 46 timelike.  The sum of all orbit normals is a strict timelike interior
vector of square 137214 and minimum wall pairing 198.  Orthogonal projection
to each negative mirror supplies a strict exact facet witness, proving that
all 47 walls are genuine and irredundant.

## Fixed divisor lattice and angles

The saturated integral fixed lattice has signature \((1,10)\), determinant
\(2^{12}=4096\), and Smith invariants

\[
  2^{10},4.
\]

After division of the form by two it is the classical determinant-two odd
Lorentzian lattice \(f_2^{10}\), equivalently Scharlau--Kirschmer `List10`
entry 2.  The primitive facets split as follows:

- six singleton roots of square \(-4\), divisibility four;
- five paired roots of square \(-2\), divisibility two;
- sixteen paired roots of square \(-4\), divisibility two;
- twenty paired roots of square \(-4\), divisibility four.

All 47 reflections preserve the saturated fixed lattice.  The complete pair
census is

\[
  m_2=490,\quad m_3=30,\quad m_4=20,\quad
  \text{parallel}=191,\quad \text{ultraparallel}=350,
\]

with no forbidden angle.

## Literal Vinberg volume certificate

Exact finite-type enumeration gives the elliptic face vector, by ranks one
through ten,

\[
[47,540,3230,12215,31212,54582,64050,47760,19910,3136].
\]

There are 216 full parabolic cusp neighborhoods.  Every one of the 19,910
elliptic rank-nine edge diagrams has exactly two endpoints:

- 12,050 finite--finite;
- 7,260 finite--ideal;
- 600 ideal--ideal.

All 568 connected affine components pass Vinberg's full zero-neighborhood
test.  Hence the chamber is finite volume and noncompact; this conclusion
does not rely merely on inheriting volume from the ambient section.

## Lineage and index

The underlying rational and integral arithmetic class is not new.  Its
maximal reflection chamber has 13 walls and three cusps.  Solomon's formula
gives

\[
  \chi(W_{47})=-\frac{527}{120960},\qquad
  \chi(W_{\max})=-\frac{31}{11147673600},
\]

so the 47-wall reflection group has exact index

\[
  [W_{\max}:W_{47}]=1,566,720
\]

and covolume \(\pi^5 527/3572100\).  Thus this is not a new lattice or
commensurability class.  No earlier source has yet been pinned for this exact
47-wall subgroup chamber; the defensible novelty claim is presently the
explicit geometric involution/fixed-divisor extraction, pending a complete
priority search.

## Replay

- `scripts/certify_borcherds_char2_2d_h10_47.py`
- `data/enumerations/borcherds_char2_2d_h10_47_certificate.json`
- `tests/test_borcherds_char2_2d_h10_47.py`

Run

```text
PYTHONPATH=src:. python3 -m scripts.certify_borcherds_char2_2d_h10_47 --check-archive
PYTHONPATH=src:. pytest -q tests/test_borcherds_char2_2d_h10_47.py
```
