# Reviewed research proof candidate: dimension at most 165

**Candidate theorem.** Every finite-sided, finite-volume real hyperbolic
Coxeter polytope has dimension **n≤165**, without arithmeticity, compactness,
a simple-cusp hypothesis or a facet-count restriction.

The [aggregate proof](fresh_nikulin_height_charge_20261015.md) and its
[independent review](fresh_nikulin_height_charge_20261015_review.md) improve
the frozen [166 candidate](fresh_nikulin_coupled_charge_20261015_dimension166.md).
The [case binder](../scripts/fresh_nikulin_height_charge_20261015_dimension165.py)
and [certificate](fresh_nikulin_height_charge_20261015_dimension165.json)
check the complete branch split, numerical threshold and recursive source
hashes. This is an internally reviewed research proof candidate. External
validation, novelty and current-best literature status remain unestablished.
The primary report's historical pending-review wording remains unchanged.

Nikulin's narrow-facet theorem supplies a facet F whose framing hyperplanes
have pairwise distance strictly below arccosh(7). If F has an ideal vertex,
the reviewed geometric chain supplies affine factors of total rank N=n−1,
actual framing walls x,y, one common height ratio M≥1, and total retained
Schur cost at most four. These are inherited geometric hypotheses; the new
finite calculations do not establish them. See the
[frozen geometric case binder](fresh_nikulin_dimension_181_candidate_20261015.md)
and [Nikulin's official ICM report, printed pp.654–655 and 657–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

For M≤31/30, unit charging fails only in listed E7/E8 support patterns.
Each failure consumes cost at least 3/2 and loses at most 2/9. Since the
total cost is at most four, at most two factors can fail. The aggregate
loss is therefore at most 40/9. Including the unchanged avoid-node loss
and Lorentz identity gives −(x,y)≥M(N−39)/9. For M≥31/30 the inherited
four-thirds estimate gives −(x,y)≥M(N−43)/9, with the same preset omissions.
At N=165 these lower bounds are respectively at least 14 and 14+1/135;
both increase thereafter. They contradict the strict narrow-facet bound
−(x,y)<14 for normals of squared length two. The ideal-vertex branch thus
gives n≤165. The scalar estimate at N=164 is insufficient to exclude that
case; no bound164 follows from this package.

If F has only finite vertices it is compact. The independently reviewed
[compact-facet argument](fresh_nikulin_narrow_20261015_compact_facet.md)
gives n≤143. These alternatives exhaust F, completing the bound165 case
split. The arithmetic ceiling29 and the noncompact arithmetic ceiling21
are separate, smaller published bounds; the new unrestricted result does
not lower them.

The independent checker uses standalone Fraction arithmetic, with no
producer imports or numerical optimizer. It reconstructs all40 exceptional
records,71 single-contact and16 multiple-contact occurrences, all580
complete-support convex-dual certificates,24 classical B/C/D classes and
two symbolic A-family tails. The proof covers arbitrary real contact
magnitudes satisfying the Coxeter gap and every affine rank. Published
classification and geometry inputs remain mathematical dependencies.

Replay: `python3 scripts/fresh_nikulin_height_charge_20261015_dimension165.py --verify`.
