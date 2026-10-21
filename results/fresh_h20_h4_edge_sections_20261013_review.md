# Review of the bounded A4-to-H4 edge-section construction

**PASS for the nine archived exact sections.** Root read the construction
source, checked its algebra and actual-facet argument, and ran a fresh
saved-result replay:

    python scripts/fresh_h20_h4_edge_sections_20261013.py --verify

All nine pass. The [review record](fresh_h20_h4_edge_sections_20261013_review.json)
binds the source, output, and numerical proposal packet. This is a source
and mathematical review plus a reused-checker replay, not an independent
implementation of the quadratic-extension arithmetic. It does not certify
the other 99 numerical proposals or exclude a continuous family.

At each initial point, the section crosses an actual ambient open edge
transversely. Its twenty incident normals project to a positive-definite
Coxeter Gram. Only one endpoint of an A4 component changes norm, changing
its adjacent angle from pi/3 to pi/5; the remaining local components are
unchanged. The full positive-definite Gram and all 210 old inequalities
justify an actual finite H4-containing section vertex and a nearby strict
interior. Finite volume is inherited from the original Klein polytope.

The field sign implementation is valid in the selected real embedding:
the positive radicand has negative norm from Q(sqrt(5)) to Q, so it cannot
be a square in Q(sqrt(5)). The prescribed positive radical therefore gives
a genuine quadratic extension, with signs decided by exact comparisons
of squares in the base field.

The Sherman–Morrison inverse and simplified edge direction are checked
by exact identities. At the next point X, all 210 inequalities hold and
exactly twenty vanish. The retained nineteen projected normals have a
positive-definite Gram. The new normal's Schur complement is

    (D.r_j)^2 / D^2 * X^2 / P^2 > 0.

Both X and P are timelike, and D is spacelike. Thus all twenty incident
normals are independent and give actual facets at a finite vertex. The
positive projected pairing of the reported pair is consequently an
**actual obtuse ridge**, not merely an incompatible pair of possibly
redundant inequalities. This rejects every one of the nine exact sections.

The construction remains a useful bounded probe of generic edge cuts;
other positions, signs, supports, and source edges remain open.
