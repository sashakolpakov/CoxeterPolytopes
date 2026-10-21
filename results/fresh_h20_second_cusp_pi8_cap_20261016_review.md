# Independent review of the pi/8 cap and cusp-to-finite edge

**PASS.** The [review source](../scripts/fresh_h20_second_cusp_pi8_cap_20261016_review.py)
and [certificate](fresh_h20_second_cusp_pi8_cap_20261016_review.json)
independently reconstruct the [alternative 152-wall chamber](fresh_h20_second_cusp_pi8_cap_20261016.md).
Generation and a fresh saved replay both pass. This branch starts at 151
walls and does not retain the optional all-dotted 152nd cap.

The preceding 151-wall physical geometry is independently reconstructed.
A general interval Gaussian inversion of the Lorentz three-frame
95,142,145 defines its exact orthogonal projection. Projecting wall 150
and then omega31 supplies an orthogonal basis of the positive two-dimensional
complement. The prescribed pi/8 product and squared norm 2 give two exact
line-circle branches. Both are checked; exactly one satisfies all old wall
gates, retains the center and cuts E5.

The accepted normal has nineteen right contacts, angle pi/8 to wall 150,
and 131 strictly dotted contacts. Projection identities establish the
prescribed equalities and norm; interval overlaps are consistency checks,
not the justification for those equalities. Strict intervals certify every
dotted relation, all 152 actual facets and preservation of five proper cusps.

A separate general four-by-four inverse, with positive-definite pivots,
constructs the new finite corner after projecting away the sixteen-root
support. Its diagram is A10+A4+A2+2A1+I2(8). All nonincident inequalities
are strict. The twenty finite old neighbors give exactly twenty candidate
supports on the new facet; the other nineteen each contain a certified
dotted pair. Hence there is exactly one new finite vertex, bringing the
known total to at least 172.

The finite corner and completed cusp nu2 share nineteen independent
positive-definite roots. Their open hyperbolic segment satisfies every
other inequality strictly, so it is an actual finite-to-ideal edge.

A different exact strict spacelike escape survives all 152 walls. Thus the
chamber still has infinite volume. The old isolated parameter and chamber
acuteness remain inherited inputs; the ambient finite-vertex census is
incomplete. No finite-volume example or new dimension is claimed. The
separate numerical residual-cone diagnostic is not a premise of this review.

```sh
python3 scripts/fresh_h20_second_cusp_pi8_cap_20261016_review.py --verify
```
