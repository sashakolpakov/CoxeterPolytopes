# Full independent review of the unified sharp-dimension manuscript

20 October 2026. Full written mathematical/source and reference review;
no build, checker, Lean, optimizer or scientific job was run in this lane.

Reviewed all 375 lines of
`manuscripts/right_angled_sharp_dimension_bounds.tex`, SHA256
`86e77ecdaf85d2d64560ad44a0226b9f5fddfa6b70f06617f93baa64d03170fa`.

**Verdict: PASS.** The manuscript gives complete human-readable arguments
for its two exclusions relative to the explicitly separated classical
geometric facts and published numerical inputs. I found no substantive
translation gap, incorrect equality branch, vertex-type restriction or
misleading mathematical feature in the TikZ figure. No source edit is
requested. This is a source review; no PDF was inspected here.

## 1. Statements, inputs and their separation

The general theorem concerns all finite-volume right-angled polytopes,
with no arithmetic or finite/ideal-vertex restriction. Its only numerical
geometric input is the six-dimensional **facet** minimum 27, attributed to
Dufour's Lemma 1 and Proposition 4 via Alexandrov's Section 4. I previously
checked that actual statement and its two underlying recurrences directly
in the [Alexandrov primary text](https://arxiv.org/html/2303.09533v1#S4).
The original Dufour full text was not retrieved in this review, and no
claim of independently reviewing its proof is made.

For sharpness I independently retrieved the official
[Potyagailo–Vinberg paper](https://ems.press/content/serial-article-files/43032).
Section 3, printed pages 66–67, constructs a finite-volume right-angled
family through dimension eight. Proposition 3.1 includes its dimension
eight member. This is the existence input used in the manuscript; the
paper's example with 27 facets in dimension six is not being mistaken for
a minimum theorem. Indeed that older paper separately poses the minimum
as a question. The manuscript correctly uses Dufour for the lower bound.

The ideal theorem additionally uses precisely Musin's cardinality bound
on S3 and Kolpakov's facet minimum with **metric** equality. The manuscript
has their correct scopes and locations: Musin's Main theorem, printed
page 4; Kolpakov's Theorem 4, printed page 1945. No equality classification
of a kissing configuration is assumed. Entire ideality, finite volume,
and right angles are explicitly inherited only by actual faces.

## 2. The common projection lemma is complete

The unit-normal convention is fixed before defining heights. A strictly
interior finite point has positive finite heights at all actual walls,
and the finite wall list has a maximum. The normalization in
`eq:projection` has square minus one, product zero with the selected
normal, and negative product with the original future point. Its future
orientation is therefore established, rather than assumed.

Equation `eq:allwall` is strict against **every** other actual wall.
It uses the stated all-wall acuteness theorem, including nonadjacent
supporting walls. Consequently the projected point belongs to the
relative interior of the actual facet. Actual ridge neighbors are
orthogonal to the selected normal and already unit in its perpendicular
Lorentz space. Their intrinsic heights are exactly the displayed old
height divided by the normalization factor. No omitted scaling enters.

The new point and the new maximum are intrinsic to the selected facet,
so the lemma can be reapplied. Monotonicity gives `eq:iterate`, including
strictness from the finite initial maximum. The proof never presumes a
unique maximizing wall or moves to a face where an ideal polygon angle
would incorrectly be called right-angled.

## 3. The general-eight proof and strict spherical endpoint

Three descents from dimension nine reach an actual finite-volume
right-angled six-face with all heights strictly below 1/sqrt3. Its
Lorentz space has dimension seven and its spatial sphere is S5 in R6.
The intrinsic normal decomposition, pairwise product estimate and
injectivity of the spatial directions are correct. The code counts
exactly its actual facets.

I reexpanded the tensor contractions in Section 3: the rank-two product
is x2−1/6; the two rank-three cross contractions are 3x and their
Kronecker-delta contraction is 24x, giving x3−3x/8. Their positive
combination is exactly the factored cubic in `eq:polynomial`.
The ordered-pair sum `eq:sos` is therefore a genuine sum of squared
Euclidean tensor norms with the stated constant N2/16.

The upper estimate uses the sign of the polynomial at every distinct
pair. Equality N=27 would force each such value to vanish; the strict
pairwise bound excludes the root 1/4, leaving only −1/2. The four-vector
squared-sum contradiction is exact. Thus the essential strengthening
from the weak upper bound 27 to the integer bound 26 is fully proved.
The lower bound and upper bound concern the same actual six-face.
The higher-dimensional tail uses an actual nine-face and makes no use
of the previous mixed-ten theorem or the ideal-four theorem.

## 4. The shortened ideal-four proof is valid

The new manuscript improves the presentation of the earlier reviewed
cusp-height proof. A single finite-point projection already gives a
four-facet with **strict** heights below one. Hence its directions are
a spherical S3 code with products below one half. Musin and Kolpakov
force the intrinsic polytope to be the regular ideal right-angled
24-cell.

Its displayed normals are outward unit vectors, and each spatial vector
has its negative among them. At a future unit point (a,w), opposite
heights at most one imply a<=1, while timelikeness gives a>=1. Thus the
point is the center and every height equals one, contradicting the
strict projected heights. The new proof needs no height propagation
around the facet graph and omits none of the necessary reasoning by
removing that longer route. Higher entirely ideal dimensions have actual
five-faces; the standard ideal four-cell supplies sharpness separately.

## 5. Figure, references and verification scope

The single figure labels the exact chain 9→8→7→6 and its successive
strict bounds correctly. Its caption explicitly distinguishes projection
arrows from polytope edges. The lower panels show the actual normal
decomposition and the resulting S5 code; they do not claim a hypothetical
wall set has been realized as a polytope. Figure dimensions and labels
match the proof. This review assesses its source semantics, not rendered
spacing or visual appearance.

The paper shares the projection foundation once, then proves the cubic
bound and the two short exclusions. It contains neither the obsolete
dimension-eleven argument nor an assumption that arbitrary ideal
four-polytopes decompose into 24-cells. The text makes no false assertion
that real geometry or the external theorems have been formalized in Lean.
Its displayed universal arguments suffice independently of the pending
exact tensor replay and conditional logic checks.
