# Rational-face obstruction and changed-incidence search checkpoint

No new compact or finite-volume Coxeter polytope, and no new dimension,
has been established. This checkpoint adds all-real obstructions and a
different compact-H9 incidence search, rather than extending the earlier
sampled parameter banks.

## Literature restriction retained

The [literature audit](literature_frontier_recheck_20261013.md) was checked
again on 14 October. Noncompact arithmetic reflection lattices can occur
only in dimensions at most 19 or exactly 21. Compact arithmetic examples
have the separate upper bound 29. Consequently the unrestricted H20/H22
construction must not be replaced by an excluded integral-lattice host
search. Allcock's published edge-walking algorithm itself works over Z;
using its traversal organization for nonarithmetic construction requires
separate realization and global finite-volume arguments.

## All-real H3 and H4 section obstruction

The [new proof](fresh_h20_h4_face_closure_20261014.md) excludes every
one-contact A4-to-H4 and A3-to-H3 section with the corresponding actual
face in the saved Borcherds source. This covers all real parameters in
3,320 A4 endpoint families and 4,230 A3 endpoint families in the saved
source-edge census.

The obstruction is geometric and arithmetic. An H4 face would be compact
of dimension 16. An exact 18-row Gram table forces it, or one of its
facets, to have only original integral roots as boundary normals. The
result would be a compact rational arithmetic polytope of dimension 16
or 15, contrary to isotropy of indefinite rational forms of rank at least
five. An H3 face either acquires precisely such an H4 subface or leads
to the same contradiction in dimensions 17 or 16. The hypothetical
section itself is not assumed arithmetic.

The [producer](../scripts/fresh_h20_h4_face_closure_20261014.py) verifies
the literal source Gram table and invariant integral lattice. A
[separate SymPy checker](../scripts/fresh_h20_h4_face_closure_20261014_review.py)
reconstructs the arithmetic, using characteristic-polynomial root counts
for the source signature and a larger parity-based adjacency table.
The proof replaces a numerical census of 152,368 candidates; those
numerical rejections are not premises of the theorem.

## Remaining single-factor cusp sections

The [single-factor report](fresh_h20_single_factor_links_20261013.md)
excludes all real shifts in 79 additional planes: 59 sums of two twin
short leaves and 20 sums of three D4 leaves. These differ from the old
59 leaf-difference planes. Their cusp links are locally Coxeter, but an
explicit global ridge is obtuse throughout each interior parameter
interval.

The same work excludes every real interior normal in each of the two
affine C2 and two affine A2 factor spans in the eighteen-cusp inventory.
A [separate SymPy check](../scripts/fresh_h20_single_factor_links_20261013_review.py)
verifies all 156 saved ambient face points against all 210 source walls,
the 79 symbolic interpolation and Gram identities, and the rank-two
auxiliary configurations. Higher-rank arbitrary single-factor directions
remain open; this is not a complete classification of single-factor links.

## Coupled local changes on source edges

The [coupled section report](fresh_h20_coupled_edge_sections_20261013.md)
tests 44 template/source-edge choices and 264 actual sections. The
templates include A1+A2 to B3, A1+A3 to F4, opposite contacts changing
A_n to A1+B_(n−1) for every n from 2 through 20, and a three-contact
A3-to-A2+A1 modification.

All 264 fail exact tests at actual neighboring finite vertices: 261
have an obtuse ridge and three have a squared cosine strictly between
1/4 and 1/2. The exact replay checks every source halfspace at the
witness point. This is a bounded result; it does not exclude all
continuous positions or all coupled local changes.

For the three selected A1+A3-to-F4 families, a separate closure argument
does cover all real positions. A facet bounding the F4 face must be
orthogonal to it. The irrational ratio of the two changed pairings then
forces that facet to be an unchanged source root, fixing the crossing
point to its rational projection onto the source-edge plane. The complete
census has 17, 11, and 14 possible points. Both branches at all 42 points
fail exact actual-ridge checks, excluding those three continuous families.
The F4 face is only required to have finite volume; compactness is not
assumed.

## Compact H9 by changing cusp incidences

The [incidence-surgery report](fresh_h9_incidence_surgery_20261013.md)
starts from three existing noncompact H9 chambers with 15 or 16 facets.
After freeing all original parallel distances, destroying every remaining
affine cusp block still requires at least 8, 7, or 5 changed finite pairs,
respectively. There are 42, 21, and 34 minimum supports. This changes the
cusps themselves rather than adding another bank of truncation walls.

Every replacement labeling from 2 through 6 on those minimum supports
is now exactly rejected: **3,131,392** assignments. Each has a nonzero
11-by-11 principal determinant modulo 241, contradicting the rank-ten
requirement for H9. The principal block avoids all free distance entries,
so the contradiction holds for every real value of those distances.
This finite-field argument includes cases where a chosen numerical
Schur basis was singular.

A different minimum-cover rule splits changed pairs between right and
newly dotted pairs. All **151** resulting incidence patterns have fixed
positive-definite rank-ten principal blocks, independently of their
free distances; those blocks cannot lie in signature (9,1).
Neither result excludes additional mutations, different incidences,
labels outside the stated set, or other compact H9 constructions.

## Next search restrictions

Do not resample the excluded H3/H4 section families, the four complete
rank-two cusp spans, or the 79 new sum planes. Remaining geometric
directions include different coupled edge profiles, higher-rank
single-factor cusp sections, and changed ambient assemblies. For compact
H9, changing more than a minimum affine cover, changing which pairs are
dotted, or changing the source introduces freedoms not settled here.
The existing H22 constructions still have certified infinite-volume
escapes; this checkpoint does not close them.

## Frozen verification

The H3/H4 producer and separate SymPy saved replays pass. The single-factor
producer and separate SymPy saved replays pass. Fresh root replays pass
for the 264 coupled sections, the 84 F4 closure candidates, the complete
modular H9 rank search, and all 151 right/dotted certificates. The
[coupled review](fresh_h20_coupled_edge_sections_20261013_review.md) and
[H9 review](fresh_h9_incidence_surgery_20261013_review.md) distinguish
source/mathematics reviews from independent arithmetic implementations.
The final audit checks all ten recorded source manifests and the local
links in the nine new reports and reviews. No live search jobs remain
from this checkpoint.
