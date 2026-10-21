# A finite-volume H9 polytope with an all-ideal seven-simplex face

The explicit rational construction below is a finite-volume Coxeter
polytope in H9 with12 actual facets,14 finite vertices and9 proper ideal
vertices. The intersection of its original facets8 and9 is a
seven-dimensional simplex whose eight vertices are all ideal.

This realizes the all-ideal alternative retained in the
[seven-simplex theorem](fresh_finite_edge_cusp_tradeoff10_20261017.md).
It does not contradict that theorem's finite-edge coefficient10: this
face has no finite edge. It is not presented as a new dimension example,
as a new polytope, or as an identification with a named classification
entry. Its purpose is to certify a necessary scope distinction using a
fully specified example. The induced seven-simplex is not itself a
Coxeter simplex.

The [exact checker](../scripts/fresh_all_ideal_seven_simplex_realization_20261017.py)
and [saved certificate](fresh_all_ideal_seven_simplex_realization_20261017.json)
use rational SymPy matrices, with no project arithmetic imports. They
enumerate all220 possible rank-nine active supports and all4096 facet
subsets. The complete vertex and face incidences, all proper cusp
components, and strict interior witnesses for every facet are saved.

## Literal normals and a strict center

Take independent roots n0,...,n9 with squared norm2. Let G be their
10-by-10 Gram matrix. Its only nonzero offdiagonal entries are−1 on

\[
 0-1-2-3-4-5-6-7-0,\qquad7-8,\qquad3-9.
\]

Thus the first eight roots form a cycle, and8,9 are leaves at opposite
cycle nodes. The exact determinant and inertia are

\[
 \det G=-16,\qquad\operatorname{inertia}(G)=(9,1,0).
\]

Use coordinates in this original-root basis and write J=G^{-1}. Add
the two normals

\[
 n_{10}=-2Je_8,\qquad n_{11}=-2Je_9.             \tag{1}
\]

The inverse has J_jj=0 for j0,...,7 and

\[
 J_{\{8,9\}}=\begin{pmatrix}1/2&-1/2\\-1/2&1/2\end{pmatrix}.
\]

Consequently the two new normals have norm2. Their only nonzero products
with original normals are (n10,n8)=−2 and (n11,n9)=−2; also
(n10,n11)=−2. All other old/new products vanish. Every distinct pair
of the12 normals therefore has product0,−1 or−2, giving respectively
right, pi/3 or parallel Coxeter gates.

Set H=−J1. The exact identities

\[
 H^2=-70,\qquad(H,n_j)=-1\ (0\le j\le9),\qquad
 (H,n_{10})=(H,n_{11})=-8                         \tag{2}
\]

give a common strict timelike center and select the future sheet.
The normals span all ten dimensions. The polytope is the intersection
of their nonpositive halfspaces with that sheet.

Every wall is actual, with an explicit strict relative-interior point

\[
 H_j=H-\tfrac12(H,n_j)n_j.
\]

It has negative norm, vanishes on wall j, and has strictly negative
products with all other normals. The latter follows directly from
(2) and their nonpositive mutual products. Its product with H is
negative, so it lies on the same future sheet. The certificate checks
all12 witnesses literally.

## Complete cone and vertex enumeration

Put y=−Gx. The original ten inequalities are y>=0. The two new ones
are x8>=0 and x9>=0, or equivalently

\[
 (Jy)_8\le0,\qquad(Jy)_9\le0.                   \tag{3}
\]

The section sum(y)=1 is bounded, and y=1 is strictly feasible before
rescaling. Hence this is a pointed full-dimensional polyhedral cone.
Every extremal ray has nine independent active inequalities. Exhausting
all binomial(12,9)=220 supports gives exactly23 rays:50 feasible support
bases representing those rays and170 infeasible supports. No support
has rank below9.

An independent closed description of the resulting ray set is:

* Eight rays y=e_j, j0,...,7.
* For i8 or9 and each j0,...,7 with J_ij<0, the ray
  y=e_j−2J_ij e_i. There are seven for each i.
* The ray y=e8+e9.

The checker compares this formula with the support enumeration exactly,
including primitive rational normalization. For x=−Jy, the first eight
rays and the last ray have square zero. Each of the other14 has

\[
 x^2=y^tJy=-2J_{ij}^2<0
\]

in the displayed normalization. Every ray has (x,H)<0. There are no
spacelike extremal rays or untested exterior vertices. In particular
the whole projective polytope lies in the closed future hyperbolic ball.

The complete face-vector, including the full polytope as the last entry,
is

\[
 (f_0,...,f_9)=(23,114,288,462,504,378,192,63,12,1).
\]

The114 edges comprise42 finite–finite,42 finite–ideal and30 ideal–ideal
edges. These counts come from the complete closed facet-intersection
lattice, not from an assumed simple-vertex incidence rule.

## Every ideal vertex is a proper cusp

For each of the nine ideal rays the checker takes every incident wall,
splits its Gram into connected components, and verifies that each is PSD
with a one-dimensional strictly positive kernel. The positive combination
of its actual normals is a positive multiple of that same future ray.
The sum of the component ranks is8 in every case:

| Ideal ray in y coordinates | Affine components | Link dimension |
|---|---|---:|
| e0,e2,e4,e6 | E8 | 8 |
| e1,e5 | D8 | 8 |
| e3,e7 | E7+A1 | 7+1 |
| e8+e9 | A7+A1 | 7+1 |

Here the names denote affine diagrams; affine A1 is the parallel-root
pair with Gram [[2,−2],[−2,2]]. The code checks the literal graph arms or
cycle, the exact rank and the positive kernel, rather than inferring the
cusp type from determinant zero alone.

In the Euclidean quotient perpendicular to the cusp ray, these positive
kernel relations make the incident normals positively span all eight
directions. Their horospherical intersection is therefore bounded and
is a product of affine simplex alcoves. It is nonempty: an interior
sequence H+t x approaches each future ideal ray while retaining all
strict inequalities. Thus all nine cusps are proper. Their standard
horospherical volume integrals have the finite tail integral proportional
to integral(exp(−8t),t0,infinity). The remaining truncated portion is
compact because the complete cone has no other ideal or exterior region.
This proves finite volume of the entire H9 polytope, not merely of a face.

## All facets and the distinguished simplex

The certificate supplies the full vertex set and face lattice of every
facet. Their finite/ideal vertex counts and facet-face-vectors are:

| Facet indices | Finite | Ideal | Face-vector through dimension8 |
|---|---:|---:|---|
| 0,1,2,4,5,6 | 12 | 8 | (20,89,196,266,238,140,52,11,1) |
| 3,7 | 13 | 8 | (21,90,196,266,238,140,52,11,1) |
| 8,9 | 7 | 8 | (15,63,140,196,182,112,44,10,1) |
| 10,11 | 7 | 2 | (9,36,84,126,126,84,36,9,1) |

Grouping entries in this table asserts their displayed counts, not an
additional unproved isomorphism between facets. Facets10 and11 each
have exactly nine vertices in dimension8 and are eight-simplexes,
with seven finite and two ideal vertices.

The original roots8 and9 are orthogonal and span a positive-definite
two-plane. Their face has equations y8=y9=0. On it, (3) is redundant,
and its complete cone is the positive span of e0,...,e7. Thus it is an
actual dimension-seven simplex with exactly eight ideal vertices. Its
face-vector is

\[
 (8,28,56,70,56,28,8,1).
\]

The induced Gram of its black side normals is the Schur complement

\[
 G_{0:8,0:8}-G_{0:8,8:10}G_{8:10,8:10}^{-1}G_{8:10,0:8}.
\]

It has inertia(7,1,0). For example, its nodes2,3 have squares2 and3/2
and product−1. Their induced angle has cosine squared1/3, strictly
between the consecutive Coxeter values1/4 and1/2 for labels3 and4.
This gives a literal reason why the all-ideal seven-simplex is not a
Coxeter simplex, despite being a face of a Coxeter polytope. Applying
the dimension restrictions for Coxeter simplexes directly to this
induced face would be invalid.

## Verification and scope

**FINAL/FROZEN.** Generation and a fresh saved replay PASS. Syntax,
recursive source hashes and local links PASS. Root supplied the
construction; frontier_literature independently reconstructed the
literal matrices, all220 active supports, all4096 facet subsets, every
facet witness and all nine positive affine cusp decompositions. A separate
agent checked the finite-volume inference and the need for full-rank
bounded cusp links; it did not duplicate the enumeration.

The example uses exact rational coordinates and no sampled labels,
approximate root selection or omitted cone frontier. It establishes
realization of the all-ideal-seven-simplex alternative, with no claim
of novelty, classification identification, a new ambient dimension
record or sharpness of the finite-edge coefficient10. Earlier frozen
proofs and examples were preserved.
