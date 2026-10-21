# Two coordinated walls close the r8 facet

The exact construction adds walls58 and59 to the frozen58-wall H9 chamber. The actual two-dimensional face with support

\[
W=\{0,1,2,3,5,6,55\}=S_6\cup\{r_8\}
\]

is now a **compact quadrilateral**, with cyclic sides7,30,59,58. Thus both open ends of the r8 facet inside the original H4+A2-supported H3 face are closed together. This is a local construction result; it does not establish finite volume of the ambient H9 chamber or compactness of the entire original H3 face.

The [script](../scripts/fresh_nonarithmetic_pi8_quad_completion_20261017.py) and [exact output](fresh_nonarithmetic_pi8_quad_completion_20261017.json) bind the [58-wall construction](fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.md) and the [complete continuous two-edge domains](fresh_nonarithmetic_pi8_two_side_domains_20261017.md). Every wall has square2. Rational outward intervals certify strict inequalities; prescribed equalities follow from the exact projection and quadratic equations below. The exact value2cos(pi/9) is enclosed using the frozen rational Machin/cosine routine, with its remainder bounds. No numerical equality is promoted to an exact contact.

## Exact recipes and native contacts

Each endpoint domain provides vectors p,A,B and a positive number rho, with

\[
A^2=B^2=0,\qquad A\cdot B=1,\qquad p\perp A,B,
\qquad p^2+\rho=2.
\]

Use the normal

\[
r(t)=p+\sqrt{\rho/2}\,(tA+B/t),\qquad t>0.
\]

For wall58 choose the [7,55] domain, profile1, labels(2,2,2,3) on(5,6,55,7), and kernel orientation−1. The equation r(t)·n49=0 has negative leading coefficient and positive constant after multiplication by t. Hence it has exactly one positive root, approximately5.8599637975.

For wall59 choose the [30,55] domain, profile3, labels(2,2,2,5) on(5,6,55,30), with the same kernel orientation. Solve r(t)·n58=−2cos(pi/9). Again the quadratic has exactly one positive root, approximately14.3730891074. The coefficients, discriminants, isolating intervals, projection vectors and null vectors are saved literally.

| New wall | Right contacts | Other finite contacts | Strictly dotted contacts |
|---|---|---|---:|
| 58 | W and49 | pi/3 to7 | 49 |
| 59 | W | pi/5 to30, pi/9 to58 | 50 |

All117 new native pairs pass. Every unspecified product is strictly below−2, so there is no unresolved finite-angle alphabet issue in this acceptance certificate. The bounded numerical event search was used only to discover this exact pair; no completeness claim is made for that search.

## Common interior and all60 facets

Both additions have strict negative product with the saved persistent old finite vertex V. The explicit point X=V+H, with H the old strict center, is future timelike and strictly inside all60 walls. For each wall n_i, the point

\[
X_i=X-\frac{X\cdot n_i}{2}n_i
\]

is future timelike, lies exactly on that wall and is strictly inside every other wall. The output contains all60 witnesses and their full product lists. Consequently all60 displayed inequalities are actual facets of a common full-dimensional Coxeter chamber.

## Complete quadrilateral face

The seven-root Gram on W is positive definite. Project the four side normals to its Lorentz3 orthogonal complement. Their squares, in cyclic order7,30,59,58, are

\[
4/3,\ 2,\ 2,\ 2.
\]

Three projected sides span this complement. A strict future point satisfies all four side inequalities, and the cone is pointed because the four covectors span the full dual space. Every extreme ray of this three-dimensional cone is therefore the intersection of at least two side planes.

The script checks all six possible side pairs. The adjacent pairs

\[
(7,30),\quad(30,59),\quad(59,58),\quad(58,7)
\]

give four future timelike rays. Each has a positive definite rank9 incident native Gram, and every nonincident native inequality is strict. The two opposite pairs(7,59) and(30,58) give spacelike intersection lines; each orientation violates one of the remaining side inequalities. Thus there are exactly four feasible extreme rays and all are finite. Every old inequality holds on their whole positive hull. This proves that the **entire actual W-face** is the compact quadrilateral, rather than just verifying four local corners.

Its intrinsic angle cosines, in that cyclic order, are

\[
\varphi\sqrt{3/8},\quad\cos(\pi/5),\quad
\cos(\pi/9),\quad\sqrt{3/8}.
\]

The intrinsic polygon need not itself be Coxeter. All native ambient contacts remain Coxeter.

## Persistent-cusp screen and scope

The complete60-wall finite-contact graph contains no irreducible affine circuit. An exact connected-growth scan over the complete Euclidean label alphabet visits670 positive definite and2 indefinite connected supports, with zero affine supports; there are no parallel pairs. Every connected affine diagram has a connected positive definite predecessor, so the growth procedure covers all irreducible affine circuits. This removes the previously encountered affine-plus-forbidden-angle obstruction for this retained graph; it is not a finite-volume criterion.

Generation and fresh saved replay pass, together with syntax, direct source hashes and report links. The output exports both computational and physical normals and a physical common center for a separate complete58-to60 cone computation. Independent verification and full ambient cone geometry are separate tasks.
