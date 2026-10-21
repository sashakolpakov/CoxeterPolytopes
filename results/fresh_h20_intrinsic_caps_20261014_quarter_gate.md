# Exact obstruction to the uniform quarter-vertex layer

The 252 proposed caps with centers `z_v/4` and common radius `3/2` cannot close the golden H20 chamber: no Coxeter-compatible subset covers two specified points. This exact statement replaces the numerical infeasibility result for this finite bank. It does not apply to face-barycenter placements; a separate search has already found a compatible set of those that covers these points.

The [independent Fraction checker](../scripts/fresh_h20_intrinsic_caps_20261014_quarter_gate.py) reads the coordinates of the [275-facet construction](fresh_h20_intrinsic_caps_20261014_golden_seed.json), without importing its arithmetic or the numerical search. Its [saved certificate](fresh_h20_intrinsic_caps_20261014_quarter_gate.json) verifies the following in `Q(phi)`.

Take the quarter points of vertices 0 and 3, whose omitted source walls are `[4,23,22]` and `[4,23,33]`. A quarter-vertex cap covers either point exactly when the corresponding old vertices have squared distance at most 36. Checking all 252 centers gives the covering sets `{0,7,14}` and `{3,10,17}`. All 505 spherical caps of the old 527-wall chamber miss both points strictly.

Every cross-pair of the two covering sets is forbidden. If the triangle-factor vertex agrees, the old squared separation is 48, giving normalized cap cosine `-1/3`, an obtuse angle. Otherwise the squared separation is `60+12phi`, giving cosine `(phi-1)/6`, strictly between 0 and 1/2. No integer Coxeter label has that cosine. These nine exact conflicts prevent simultaneous coverage of the two points.

Consequently at least one probe remains strictly outside all caps in any compatible subset of this finite bank. Both probes lie strictly inside the vertical seed alcove; the remaining open ideal-boundary region implies infinite volume. This is a finite-family obstruction, not a general obstruction to the seed, quarter-scale caps, or changed chambers.

Generation and saved replay pass:

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_quarter_gate.py --verify
```
