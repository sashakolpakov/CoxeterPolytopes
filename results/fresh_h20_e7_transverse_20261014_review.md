# Narrow review of the single-weight E7 transverse families

**2026-10-14 — PASS.** Source and proof review of the [transverse gate](fresh_h20_e7_transverse_20261014.md) found no issue. A fresh primary `--verify` passes. The [small separate Fraction checker](../scripts/fresh_h20_e7_transverse_20261014_review.py) and its [saved certificate](fresh_h20_e7_transverse_20261014_review.json) also pass; they bind the frozen primary source and certificate, independently check the rational norm substitutions and null-contact classification, and do not duplicate the ambient multiquadratic reconstruction.

At each certified finite face vertex, releasing one E7 wall leaves nineteen independent positive-definite incident normals. Their two-dimensional orthogonal complement is spanned by the corresponding E7 fundamental weight and the timelike vertex ray. Thus the displayed two-coefficient family exhausts normals with those nineteen right contacts. Acuteness fixes the weight coefficient to be nonnegative, and norm two makes it nonzero. Positivity of every entry of the E7 inverse implies that the weight has negative pairing with every strict interior point. Together with the future orientation of the vertex, this excludes the negative square-root sheet and the zero-root boundary. The positive sheet is therefore necessary for a valid added facet.

The released root and cap both have norm two, so their nonzero contact coefficient belongs to `{1} union [sqrt2,infinity)`. This relaxation covers every integer finite label, parallel pair, and dotted pair. For the positive sheet, `u=sqrt((d-2/c²)/q)` increases with `c`. The three proof branches therefore have the stated direction: an escape bound at `c=1` persists above it; an old-wall lower bound on `u` prevents escape cutting; and in the final three cases the strict escape bound at `sqrt2` removes all larger coefficients while the old-wall inequality removes `c=1`. No finite label or parameter cutoff is used.

For a null residual, the norm equation gives `c²=2/d` independently of its real null coefficient. The seven exact weight squares give ten profiles in `(0,1)`, two in `(1,sqrt2)`, and two with `c=1,d=2`. In each of the last two profiles the pairing with the retained matching filler is identically `+2`. A distinct cap is therefore incompatible; the duplicate filler is already known to miss the escape. This argument does not require a bound on the remaining real parameter.

The saved escape is future spacelike and strict against the retained 87 walls. Nonpositive pairing with every admissible new wall in these families suffices for the finite-volume conclusion: for any finite coherent selection, perturb toward that selection's strict timelike interior point. All inequalities become strict and positive norm persists. The proof properly restricts this perturbation statement to a finite selection and does not assume an unproved uniform margin over infinitely many parameters.

The conclusion covers exactly 308 single-weight/finite-ray and 14 single-weight/null-ray families. Mixed E7 weights, other residual directions, and deletion of retained walls remain open. The review consists of mathematical/source review, a fresh replay using the primary arithmetic engine, and a separate small rational parameter check—not an independent complete chamber or field-arithmetic census.

```sh
python scripts/fresh_h20_e7_transverse_20261014.py --verify
python scripts/fresh_h20_e7_transverse_20261014_review.py --verify
```
