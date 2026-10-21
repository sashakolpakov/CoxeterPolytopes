# Reviewed research proof candidate: dimension at most 164

**Candidate theorem.** Every finite-sided, finite-volume real hyperbolic
Coxeter polytope has dimension **n≤164**, without arithmeticity, compactness,
a simple-cusp hypothesis or a facet-count restriction.

The [residual-cost proof](fresh_nikulin_defect_frontier_20261015.md) and
[independent review](fresh_nikulin_defect_frontier_20261015_review.md) sharpen
the frozen [165 candidate](fresh_nikulin_height_charge_20261015_dimension165.md).
The [case binder](../scripts/fresh_nikulin_defect_frontier_20261015_dimension164.py)
and [certificate](fresh_nikulin_defect_frontier_20261015_dimension164.json)
verify the branch split and recursive source hashes. This is an internally
reviewed research proof candidate; external validation, novelty and current-best
literature status remain unestablished. Historical pending-review wording
in the frozen primary package is superseded by the separate review.

Use the inherited narrow-facet geometry, actual framing walls, common null
height ratio M≥1 and total retained-support Schur budget four. If the narrow
facet has an ideal vertex, its affine factors have total rank N=n−1. For
M≤21/20, the same fifteen possible E7/E8 exceptional support occurrences
are the only factors that can fail unit charging. Three actual failures
would cost at least9/2, so at most two are possible.

Two ordinary E8 failures consume cost at least28/9, leaving strictly less
than one everywhere else. Each retained inverse diagonal is at least1/2,
so this residual product has at most one nonzero retained coordinate. The
new all-rank lemma bounds its energy loss by2/3. The other defect cases are
also checked, giving total loss at most38/9. Including the unique avoid-node
penalty, the Lorentz pairing obeys −(x,y)≥M(N−38)/9. For M≥21/20 retain the
reviewed four-thirds estimate −(x,y)≥M(N−43)/9, using identical preset omissions.

At N=164 the two lower bounds are respectively at least14 and14+7/60,
and both increase with N. They contradict the strict narrow-facet bound
−(x,y)<14 for squared-length-two normals. Thus the ideal-vertex branch
excludes n≥165. If the narrow facet has no ideal vertex it is compact,
and the frozen [compact-facet proof](fresh_nikulin_narrow_20261015_compact_facet.md)
gives n≤143. These alternatives complete the candidate bound164.

The estimate covers every affine rank and every prescribed avoid-node case,
with arbitrary real nonzero Coxeter contacts at least one. It does not use
a finite-label cutoff. The independent review checks the widened height
window, the exceptional singleton domains, symbolic classical families and
aggregate arithmetic. The geometric dependencies are inherited from the
[reviewed 165 case split](fresh_nikulin_height_charge_20261015_dimension165.md);
the local certificate does not itself prove the narrow-facet theorem.

The near-range scalar estimate at N=163 is below14, so this package claims
no bound163. Keeping the two separate endpoint budgets is an ongoing
investigation with its own files. The arithmetic ceiling29 and noncompact
arithmetic ceiling21 are separate published results, unaffected by164.

Replay: `python3 scripts/fresh_nikulin_defect_frontier_20261015_dimension164.py --verify`.
