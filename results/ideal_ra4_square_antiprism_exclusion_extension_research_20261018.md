# Excluding square antiprisms from the Oct/D facet class

18 October 2026. New own written proof; independent review pending.
No computation, build, Lean process, literature search or frozen-file edit.

**Precise extension.** Let P be an entirely ideal finite-volume
right-angled hyperbolic four-polytope. Suppose every facet is isometric to
the regular ideal octahedron O, its triangular-face double D, or an ideal
right-angled polyhedron combinatorially a square antiprism. Then no facet
has the third type. Consequently the existing internally reviewed Oct/D
decomposition and chain theorem applies unchanged to this enlarged class.
This proves no statement about an arbitrary ideal right-angled
four-polytope or universal volume divisibility.

The essential new ingredient is a triangular ridge whose three cusp
corners demand inconsistent types for its one neighboring facet. A formal
box pattern alone would not give the contradiction.

## 1. Elementary rigidity of the ideal square antiprism

This step avoids importing a general three-dimensional rigidity theorem.
Use Lorentz product of signature (-,+,+,+). Let b,c be outward unit
normals of the two square bases. Acuteness and orthogonal projection give
the distinct-normal gap0 or at most−1: if their product g lies in(−1,0),
their two-by-two Gram matrix is positive definite with nonnegative
inverse. Project a strict timelike interior point to the two walls.
The added wall coefficients are positive, so all other inequalities
remain strict; the resulting actual ridge would be nonright. Product0 gives an actual
ridge by positive-Gram projection. Product−1 makes b+c a future null
feasible point by acuteness, hence a shared actual cusp. Neither occurs
between the disjoint bases. Thus they are ultraparallel, and coordinates
may be chosen so

    b=(0,0,0,1),   c=(sinh L,0,0,−cosh L),   L>0.

The four upper cusps, in their cyclic order, are

    ξ_i=(1,cos θ_i,sin θ_i,0).

At ξ_i, the base b and the opposite lower lateral triangle normal ℓ_i
are opposite walls in its rectangular cusp link. Thus
ℓ_i=μ_i ξ_i−b with μ_i>0. This triangle is orthogonal to the lower
base c, so μ_i=coth L=:k>1, independent of i.

Consecutive ℓ_i and ℓ_(i+1) meet only at a lower cusp and are an opposite
pair there. Therefore

    −1=ℓ_i·ℓ_(i+1)=k²[cos(θ_(i+1)−θ_i)−1]+1.

Every cyclic angular gap is α or2π−α for one α∈(0,π), with
sin(α/2)=1/k. A gap2π−α would make the sum of all four positive gaps
at least2π+2α. Hence all four gaps equalα=π/2, so k=√2.
The upper base is regular; its four neighboring triangle planes are
determined by its edges and orthogonality to b. The lower triangle planes
ℓ_i are already fixed, as is c. All ten supporting planes are therefore
determined up to isometry.

This proves rigidity for the combinatorial square antiprism itself.
It does not infer any four-dimensional realization from that combinatorics.

## 2. An explicit model and its cusp-side signature

Put h=√2−1 and r²=2h=1−h². A representative of the unique antiprism
has upper vertices

    v_j=(1,r cos(jπ/2),r sin(jπ/2),h)

and lower vertices

    w_j=(1,r cos(jπ/2+π/4),r sin(jπ/2+π/4),−h),
    j=0,1,2,3.

The two base normals are B_±=(h/r,0,0,±1/r). The upper lateral triangle
T=(v_0,v_1,w_0) has normal

    N_0=(1/r,1,1,h/r).

At v_0 the four incident normals are

    B_+, N_+=(1/r,1,1,h/r), N_-=(1/r,1,−1,h/r),
    M=(1/r,√2,0,−h/r).

They are unit; neighboring pairs are orthogonal; the opposite pairs have
product−1. Their sums are

    N_++N_-=(2/r)v_0,
    B_++M=(√2/r)v_0.

At a cusp, coefficients of a common null vector in opposite unit-normal
sums are proportional to the corresponding Euclidean widths. This follows
directly from vertical upper-half-space walls: opposite unit normals have
sum equal to their signed separation times one fixed cusp null vector.
The rectangle therefore has aspect√2:1. The side in the square base
B_+ spans the larger width, so it is a **long side**. It is the only
quadrilateral side at this cusp. By vertex transitivity the same signature
holds at all eight cusps:

    A4: aspect√2:1; exactly one quadrilateral side, a long side.

The triangle T is a short side of this rectangle at v_0 and v_1. At w_0,
its normal N_0 is opposite the lower base normal B_−. Indeed

    B_−+N_0=(√2/r)w_0,

whereas the other opposite pair has sum(2/r)w_0. Thus T is a long
side at w_0. Every lateral triangle has the signature

    short at its two vertices on one square base;
    long at its third vertex on the other base.

These are actual intrinsic cusp rectangles and triangle planes of the
rigid hyperbolic antiprism, not arbitrary labels assigned to a box.

For comparison, the already explicit O/D models have signatures:

|Facet cusp|Aspect|Quadrilateral sides at that cusp|
|---|---|---|
|O|1|none|
|D, six square cusps|1|exactly one|
|D, three nonsquare cusps|2:1|both long sides|
|A4, every cusp|√2:1|exactly one long side|

## 3. Exhausting the ambient cusp boxes

An ideal right-angled four-dimensional cusp has a rectangular box link.
Write its side lengths x≥y≥z and normalize z=1. Every rectangular face
has aspect in{1,√2,2}. The complete ordered list is

    1:1:1, √2:1:1, √2:√2:1, 2:1:1, 2:√2:1, 2:2:1.

Each box edge is the horospherical corner of one actual triangular or
quadrilateral ridge. Its label must agree on both incident box faces.

* In2:2:1, four nonsquare D faces force all four edges on each square
  end to be quadrilateral. A square O cusp allows zero, a square D cusp
  only one. Impossible.
* In√2:√2:1, the four nonsquare A4 faces each require one long
  quadrilateral side. All four such sides lie on the two square ends,
  so the ends receive four quadrilateral incidences in total but allow
  at most two. Impossible.
* In2:√2:1, the two2-by-1 faces are D. Their long sides force both
  long sides of each2-by-√2 A4 face to be quadrilateral. A4 permits
  exactly one. Impossible.
* In2:1:1 the four nonsquare faces are D; their long sides are
  quadrilateral and their short sides triangular. The two square ends
  therefore have no quadrilateral side and are O.
* In√2:1:1 the four nonsquare faces are A4. Every short side is
  triangular, so both square ends are O. The long quadrilateral edges
  form one of the two perfect matchings of the four-face belt, because
  every A4 rectangle requires exactly one.

A cube has only square faces and hence contains no A4 cusp.
Consequently every cusp on a proposed A4 facet has shape√2:1:1,
four A4 belt facets, and two O ends. In this box a short side of an
A4 rectangle meets an O end; a long side meets another A4 belt facet.
This classification is a necessary condition for an actual polytope.
It does not assert that the formal four-A4/two-O box can extend to
whole facets; the next step proves that it cannot in this facet class.

## 4. The actual triangular-ridge contradiction

Choose a lateral triangular face T of a proposed A4 facet F. There is
one other ambient facet F' meeting F across the relative interior of T.
At either base-edge vertex of T its cusp side is short. The ambient
box calculation therefore forces F' to be an O end facet.

At the third vertex of the very same T its cusp side is long.
The ambient box calculation forces that same F' to be an A4 belt facet.
The types O and A4 are distinct. This contradiction excludes F.

The argument uses both the intrinsic triangle signature and the unique
ambient neighbor across an actual ridge. A box-face coloring alone would
miss this obstruction. We have not constructed a finite-volume
four-polytope with A4 facets; none is asserted.

## 5. What survives for further octahedral gluings

The A4 exclusion is more local than the old surface theorem. It remains
valid if D is replaced or augmented by any family X whose actual cusp
signatures are exactly those in the two D rows above, provided a square
cusp with zero quadrilateral sides occurs only in O. The same box list
and the same triangular-ridge contradiction apply. This is a signature
criterion, not a claim that a proposed gluing has those signatures.

The global octahedral decomposition does not follow from this signature
criterion alone. The old surface proof uses one actual internal ideal
triangle in each D, with its three sides lying in quadrilateral interiors.
A larger octahedral gluing can have several internal triangles, or cusp
rectangles with a new ratio, so the construction and the opposite-cusp
injection must be proved again for that specified geometry. Even if a
24-cell decomposition is obtained, its old path-rigidity proof uses the
fact that a boundary facet contains at most two octahedra. Permitting a
three-octahedron boundary facet removes precisely that contradiction.
It would be invalid to carry the chain conclusion over unchanged.

For the present theorem, no such step is needed: A4 is excluded, and the
facet class returns exactly to O/D. The existing internally reviewed
theorem then gives its chain, uniqueness within the restricted class and
volume formula. No universal decomposition or divisibility follows.

## 6. A consequence for the ideal5/6 carrier problem

The previously proved universal identity in
[the five-face carrier packet](right_angled_five_face_carrier_20261018.md) is

    Σ_(four-facets Q of an ideal RA5) κ(Q)=48,
    κ(Q)=3f_0(Q)−4f_3(Q)+24.

For an m-cell chain the proved counts give
f_0=18m+6, f_3=14m+10, so κ=2−2m≤0.
If every three-face of an ideal RA5 had the explicit O/D/A4 spectrum,
each of its four-facets would be a chain by the present extension.
Their kappa sum would be nonpositive, contradicting48. Thus every
ideal RA5 contains a three-face outside this spectrum. The same holds
for an ideal RA6 by applying the statement to a five-facet.

The earlier carrier constraints require a positive-kappa four-face,
beyond a24-cell or its double, and the new chain theorem already
rules out a positive carrier with only O/D three-faces. The extension
adds A4 to the excluded three-face spectrum. This consequence is
implicit in combining those identities with the new exclusion; it is
not an exclusion of ideal dimensions5 or6.

To translate “outside O/D/A4” into “at least ten vertices,” one needs
the separate complete classification of ideal right-angled
three-polyhedra with at most nine vertices. This report does not
establish or invoke that published classification, so it does not
silently make the numerical translation.

## 7. Translation and verification scope

The separate TeX appendix source is
[ideal_ra4_square_antiprism_extension_appendix.tex](../manuscripts/ideal_ra4_square_antiprism_extension_appendix.tex).
It is an input fragment with theorem/proof environments and a two-panel
cusp/triangle diagram; the existing Oct/D manuscript and PDF are unchanged.
All new calculations above are symbolic written mathematics. There is no
new exact checker, machine replay, Lean proof or external-priority claim.
Independent mathematical review of this new extension is still required.
