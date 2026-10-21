# The characteristic-two order-3 divisor section: H11/22

Status: exact finite-volume noncompact Coxeter polytope; abstract chamber
known, algebraic extraction newly explicit in this project.

## Construction

Start with the complete 210-wall chamber of the supersingular
characteristic-two K3 surface of Artin invariant one.  Its Neron--Severi
lattice is

\[
  U\oplus D_{20}(-1),
\]

and the wall Gram has rank 22 and inertia \((1,21,188)\).  The full weighted
wall diagram has order 241920.  Its ATLAS class `3b` is the unique order-three
class with fixed rank 12 (class size 672).  Because the geometric
Dolgachev--Kondo subgroup has index two, every odd-order element, including
`3b`, lies in it.

The 210 walls split into 12 fixed walls and 66 three-cycles.  Restriction to
the fixed divisor space gives 78 wall-orbit normals.  Exactly 22 have negative
square and meet hyperbolic space; 50 are isotropic and six timelike, hence are
causally inactive.  The sum of all orbit normals is a strict timelike
interior vector of square 137214 and minimum wall pairing 99.  Orthogonal
projection of this vector to each of the 22 negative mirrors gives an exact
strict facet witness, so no wall is merely formal or redundant.

## Lattice and angles

The saturated fixed lattice has signature \((1,11)\), determinant
\(-3^5=-243\), and Smith invariants

\[
  1^7,3^5.
\]

The 22 primitive facets are

- 12 roots of square \(-2\), divisibility one;
- 10 roots of square \(-6\), divisibility three.

Every corresponding reflection is integral on the saturated fixed lattice.
The complete pair census is

\[
  m_2=186,\qquad m_3=35,\qquad m_6=10,
\]

with no parallel, ultraparallel, or forbidden pair.

## Literal Vinberg volume certificate

Exact finite-type enumeration gives the elliptic face vector, by ranks one
through eleven,

\[
[22,231,1480,6320,18704,38986,57182,57785,38140,14580,2304].
\]

There are 72 full parabolic cusp neighborhoods.  Every one of the 14,580
elliptic rank-ten edge diagrams has exactly two endpoints:

- 10,965 finite--finite;
- 3,414 finite--ideal;
- 201 ideal--ideal.

All affine zero-neighborhoods pass Vinberg's full parabolic test.  Therefore
the chamber is finite volume and noncompact.

## Lineage

The ambient H21/210 chamber is due to Borcherds, with the K3 realization due
to Dolgachev--Kondo; Alexeev's ancillary file supplies the pinned matrix used
here.  The resulting abstract H11/22 chamber is not new: its lattice, complete
elliptic face vector, and 72 cusps match Scharlau--Kirschmer `List11`, entry
27.  A separate height enumeration identifies the 22 walls with the complete
simple-root system, hence the containment index is one.  What remains
priority-open is the explicit order-three K3 fixed-divisor interpretation,
not the abstract polytope.

## Replay

- `scripts/certify_borcherds_char2_3b_h11_22.py`
- `data/enumerations/borcherds_char2_3b_h11_22_certificate.json`
- `tests/test_borcherds_char2_3b_h11_22.py`

Run

```text
PYTHONPATH=src:. python3 scripts/certify_borcherds_char2_3b_h11_22.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_borcherds_char2_3b_h11_22.py
```
