# Independent review of the ideal distance-nine replacement

**PASS: mathematical/source review and independent exact reconstruction.** The
[producer report](fresh_distance8_ideal_corner_replacement_20261017.md), final
SHA256 `83b6c9eb437a0eb39ac3caf1fbd482fba313cae120824b311006100273d1417a`,
proves the stated local quadrilateral replacement without assuming simple
cusps or excluding small affine diagrams. No mathematical gap was identified.
The [review checker](../scripts/fresh_distance8_ideal_corner_replacement_20261017_review.py)
and [saved reconstruction](fresh_distance8_ideal_corner_replacement_20261017_review.json)
import no producer or other project arithmetic code.

The actual ideal corner consists of its varying affine factor and a positive
definite spectator forest. The opposite sides do not contain that ideal
point, so their heights against the positive null combination are strictly
positive even when the full cusp is nonsimple. The four residual side normals
span a Lorentz three-space. Consequently their four-by-four determinant must
vanish. Eliminating the orthogonal spectator forest changes the normalized
determinant by a positive semidefinite quadratic form and cannot invalidate
the strictly positive contradiction.

The one-root extension lists are complete by the stated finite/affine
classification. A long A path admits only the eight displayed oriented
profiles. In the classical tree case, feature distance greater than14 puts
both marks on the simple spine with long B/D outer components. Their retained
forks or doubled bonds prohibit other classical extensions; exceptional
types are too small. Joining both long components gives the corresponding
affine tree through simple long-end bonds. Including the single-component
affine profiles, even where the actual edge condition could exclude them,
only enlarges the necessary domain.

The height bounds use only nonnegative native contacts and the Coxeter gap
between zero and one. In the cycle case, attachments at positions
`j_c,j_d` on a common arc give a white-interior path of length
`P-j_c-j_d+2`; excluding paths of length at most8 gives exactly
`P >= max(9,7+j_c+j_d)`. The determinant expression is nonincreasing in
the reciprocal arc lengths and concave in the two reciprocal heights.
Its rectangle checks therefore cover all real heights and all ranks.
For the feature-free tree, the all-rank endpoint/second-node inverse
diagonals1 and2 give the stated profile table. The separated contacts on
the middle path give the cross term `l*m/P`, and each single-root cost is
at most4, giving the strict gap at `P=9`.

The independent SymPy checker reconstructs all eight cycle profiles and all256
rectangle values, matching the uniform minimum `5/81`. It separately
reconstructs all eight tree costs and all64 ordered pair gaps, matching
the minimum `1/9`. It computes the full symbolic four-root determinant
directly, derives the spectator correction by substitution, and verifies
the four path inverse identities symbolically with an indeterminate length.
These are independent arithmetic checks, not calls to the producer's field
class or determinant routine.

Generation and a fresh saved replay pass. The review retains the published
finite/affine classification and the inherited actual factor-fiber geometry
as explicit inputs; it does not claim an independent classification census.
The conclusion is the local replacement lemma. An ambient dimension bound
requires a separately checked payment and scalar composition. External
refereeing and proof-assistant formalization are not claimed.
