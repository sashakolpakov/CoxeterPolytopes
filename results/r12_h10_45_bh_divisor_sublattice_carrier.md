# TOOT — two order-four K3 carriers for the 45-wall `H10` polytope

## Positive result

The previously certified finite-volume, noncompact 45-facet Coxeter
polytope in \(\mathbb H^{10}\) now has two genuine algebro-geometric
carriers.  On a very general member of each of two explicit Weierstrass K3
families, the Néron--Severi lattice contains an automorphism-stable,
full-rank divisor sublattice

\[
M\cong U\oplus D_4(-1)^2\oplus A_1(-1)^2
\]

whose lattice-action pair is integrally conjugate to the archived `H10/45`
pair.  Transporting its 59 roots and taking the fixed section produces the
same 45-facet Coxeter polytope.

This is a new explicit algebraic realization in the project.  It is not a
claim that the abstract Coxeter polytope, its arithmetic class, or the K3
families are new to the literature.

| K3 carrier | very general \(\operatorname{NS}(X)\) | index of \(M\) | defining parity condition in the canonical BH basis |
|---|---|---:|---|
| Brandhorst--Hofmann `0.4.4.2` | \(U\oplus D_8(-1)\oplus A_1(-1)^2\) | 2 | coordinate 10 is even |
| Brandhorst--Hofmann `0.4.4.1` | \(U\oplus E_8(-1)\oplus A_1(-1)^2\) | 4 | coordinates 10 and 12 are even |

Thus the failed primitive special-period lane was too restrictive: the
right carrier is a nonprimitive divisor subsystem inside a slightly larger
Picard lattice.

## Explicit surfaces and order-four actions

The first family has equation

\[
y^2=x^3+P(t^2)x+tQ(t^2),
\qquad
\sigma(x,y,t)=(-x,iy,-t),
\]

with very general Néron--Severi lattice
\(U\oplus D_8(-1)\oplus A_1(-1)^2\).  The second has

\[
y^2=x^3+A(t^2)x+tC(t^2),
\qquad
\sigma(x,y,t)=(-x,iy,-t),
\]

with very general Néron--Severi lattice
\(U\oplus E_8(-1)\oplus A_1(-1)^2\).  The parameter formulae and nonempty
Zariski-open conditions are recorded in the separate family certificates.
In both cases \(\sigma\) has order four and multiplies the holomorphic
two-form by \(-i\).  Its image on the Picard lattice has order two and is the
pinned Brandhorst--Hofmann action.

The Brandhorst--Hofmann version-2.1 database rows are pinned byte for byte:
row 40 (`0.4.4.2`) and row 39 (`0.4.4.1`) of
`K3Groups/results/purely_ns/order04.txt`.  A replay from the external ZIP
reconstructs the full K3 lattice, the generic Néron--Severi lattice, and the
Picard action for both rows.  The source is
[*Existence and uniqueness of automorphisms of K3 surfaces*](https://doi.org/10.1017/fms.2023.50),
with the associated [database](https://doi.org/10.5281/zenodo.10410278).

## Exact lattice-action identification

For each carrier the parity sublattice has signature \((1,11)\), determinant
\(-64\), Smith diagonal

\[
(1^6,2^6),
\]

and the exact discriminant quadratic-value histogram of
\(U\oplus D_4(-1)^2\oplus A_1(-1)^2\).  Its Picard involution has
characteristic polynomial

\[
(\lambda-1)^{11}(\lambda+1).
\]

The fixed lattice has signature \((1,10)\), determinant 64, Smith diagonal
\((1^6,2^4,4)\), while the coinvariant lattice is exactly
\(\langle-4\rangle\).  More than these coarse invariants are checked.  On
the five nontrivial Smith coordinates \((2,2,2,2,4)\), the map

\[
\begin{pmatrix}
0&0&0&1&0\\
1&0&0&0&0\\
0&1&0&1&1\\
0&0&1&0&0\\
0&0&2&0&1
\end{pmatrix}
\]

is tested on all 64 discriminant classes and preserves the quadratic form
modulo \(2\mathbb Z\).  Together with the identity on the
\(\langle-4\rangle\) discriminant group it carries every primitive
fixed/coinvariant glue pair.

The fixed lattice has rank 11 and discriminant length 5.  Nikulin's
indefinite uniqueness and surjectivity theorem therefore lifts the displayed
discriminant isometry to an integral fixed-lattice isometry; compatibility
with the primitive eigen-glue extends it to the full lattice-action pair.
This proves integral conjugacy, rather than merely matching determinant or a
quadratic-value histogram.

## Why the 59 roots can be made effective

The archived integral vector

\[
h=(3,3,3,0,0,2,2,6,6,3,1,7)
\]

is fixed by the Picard involution, has \(h^2=40\), and pairs positively with
all 59 wall roots; the minimum pairing is 1.  In either ambient Picard
lattice the coinvariant lattice is the rootless lattice
\(\langle-4\rangle\).  Consequently the fixed positive cone is not trapped
in an extra root hyperplane.  Perturbing \(h\) in that cone gives an
invariant ambient Weyl chamber \(C\).

There is one Torelli subtlety worth making explicit.  The displayed
Weierstrass automorphism \(\sigma\) supplies the required full K3 Hodge
action, but its pre-existing ample chamber need not equal \(C\).  Choose a
Néron--Severi Weyl element \(w\) carrying \(C\) to the actual ample chamber.
A product of square-minus-two reflections acts trivially on the
Néron--Severi discriminant group, hence extends as the identity on the
transcendental lattice.  The conjugate

\[
w\sigma^*w^{-1}
\]

is an order-four Hodge isometry preserving the ample chamber.  Strong
Torelli realizes it as an automorphism \(\tau\) of the same explicit K3
surface.  Every \(w\)-transported wall root pairs positively with an ample
class, so K3 Riemann--Roch makes all 59 classes effective simultaneously.

The output is therefore an effective invariant **divisor-subsystem**
chamber.  Extra roots of the larger ambient Néron--Severi lattice may refine
it.  The certificate deliberately does not call the 59-root chamber the
full nef chamber of the carrier surface.

## Literal Vinberg conclusion

The fixed Coxeter section reuses the frozen exhaustive Vinberg audit, not a
heuristic vertex count.

| datum | exact value |
|---|---:|
| dimension | 10 |
| facets | 45 |
| elliptic rank-9 edges | 7,332 |
| finite vertices | 1,048 |
| ideal vertices/cusps | 93 |
| connected affine components audited | 302 |
| bad edge endpoints | 0 |
| bad affine zero-neighbourhoods | 0 |

Every elliptic codimension-one edge has exactly two finite or ideal
endpoints, and every maximal affine component has a valid full
zero-neighbourhood.  Hence the polytope has finite volume and is noncompact.

## A reusable positive criterion

The construction is not an isolated lucky basis change.  Let \(F\) be the
rank-11 fixed lattice of the `H10/45` action and let \(\gamma\in A_F\) be
the nonzero fixed component of the primitive eigen-glue.  A finite-index
Picard carrier is produced by an even isotropic subgroup

\[
H\subset A_F,\qquad H\perp\gamma.
\]

The ambient index is \(|H|\).  Exhausting the 64-element discriminant group
gives exactly:

- seven admissible subgroups of order 2; every extension is isometric to
  \(U\oplus D_8(-1)\oplus A_1(-1)^2\);
- three admissible subgroups of order 4; every extension is isometric to
  \(U\oplus E_8(-1)\oplus A_1(-1)^2\).

For every extension the code constructs the integral overlattice, verifies
evenness and hyperbolic signature, checks an exact discriminant-form
isometry to the reference lattice, and invokes indefinite uniqueness.  This
turns the positive example into a small finite search criterion for further
algebraic carriers.

## Reproducibility

- exact certificate:
  [`r12_h10_45_bh_divisor_sublattice_carrier.json`](../data/enumerations/r12_h10_45_bh_divisor_sublattice_carrier.json)
- certifier:
  [`certify_r12_h10_45_bh_divisor_sublattice_carrier.py`](../scripts/certify_r12_h10_45_bh_divisor_sublattice_carrier.py)
- regression tests:
  [`test_r12_h10_45_bh_divisor_sublattice_carrier.py`](../tests/test_r12_h10_45_bh_divisor_sublattice_carrier.py)
- explicit \(D_8\) family:
  [`ud8a1x2_order4_weierstrass_family.md`](./ud8a1x2_order4_weierstrass_family.md)
- explicit \(E_8\) family:
  [`ue8a1x2_order4_weierstrass_family.md`](./ue8a1x2_order4_weierstrass_family.md)

Replay the frozen certificate with:

```text
PYTHONPATH=.:src python3 scripts/certify_r12_h10_45_bh_divisor_sublattice_carrier.py --check-archive
```

To replay the two pinned Brandhorst--Hofmann rows from the original database
ZIP as well, add `--brandhorst-hofmann-zip K3Groups-v2.1.zip` and direct the
output to a temporary path.
