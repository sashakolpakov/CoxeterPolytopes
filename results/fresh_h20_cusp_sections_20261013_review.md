# Independent review of the cusp-section pencils

The [118 oriented pencil certificates](fresh_h20_cusp_sections_20261013.json), representing 59 literal normal planes at 14 saved cusp representatives, pass independent mathematical/source review. The primary saved replay passes. A [separate standard-library Fraction checker](../scripts/fresh_h20_cusp_sections_20261013_review.py), importing no producer code, also replays all 118 actual-ridge constructions and full cusp links; its [certificate](fresh_h20_cusp_sections_20261013_review.json) binds the source and data hashes.

The all-parameter reduction is sound. In each plane, a spacelike normal can be scaled to `e0+t*delta`; the remaining projective normal is null and cannot define an interior hyperbolic section. The positive wall combination at `t=1/2` excludes that endpoint and every larger parameter. Reversing the leaves covers negative parameters. The interpolated ridge points establish genuine interior sections for `0<t<1/2`; at zero the common strict center lies in the section.

For positive `t`, the projected leaf difference is `-2t*delta`. Its positive pairing with every future timelike point makes the retained leaf inequality stronger, so the omitted leaf is globally redundant. The connected acute projected Gram components have positive kernels. Rescaling by those kernels gives weighted graph Laplacians, hence the asserted corank one and total rank 19. Their positive null multiples remain positive throughout the open interval. This certifies a bounded Euclidean Coxeter cusp link, not just formal allowed pairwise angles.

I checked both auxiliary-wall cases independently. The actual ridge is obtained by positive interpolation of two future timelike ambient face points, whose zero sets intersect in exactly the asserted two walls. Its projected off-diagonal pairing is `-(1/2-t)`, and the auxiliary squared norm is `7/4+s*t-t²`, with `s=1` or `-1`. Therefore subtracting four times the squared-cosine numerator from its denominator gives

```text
3/4 + (4+s)*t - 5*t².
```

This concave polynomial is positive at both interval endpoints. Thus the actual ridge has squared cosine strictly between zero and `1/4`, excluding every integer Coxeter label. At `t=0` the two leaf projections coincide positively, and the same ambient face gives a genuine ridge after removing the duplicate, with squared cosine `1/7`. No singleton parameter or boundary exception is discarded numerically.

The separate replay checks the literal supporting identities, all face-point wall signs, the interpolation equations, both cosine formulas, the zero-parameter duplicate, and every affine component's kernel and connectedness. It also checks both orientations of all 59 leaf pairs. The source review checks the twin census and the real-parameter completeness argument; the replay does not independently reproduce the previous full H21 cusp enumeration.

```sh
python3 scripts/fresh_h20_cusp_sections_20261013.py --verify
python3 scripts/fresh_h20_cusp_sections_20261013_review.py --verify
```

Both pass. The scope is exactly the 59 declared normal planes and their stated symmetry transports. Other cusp-adapted directions and H20 constructions remain open. No new Coxeter polytope is established.
