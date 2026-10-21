# Review of the integer-gap inner escape

**PASS.** Source review and a fresh saved replay of
[the primary checker](../scripts/fresh_h20_integer_inner_gap_20261014.py)
found no issue. A [separate Fraction checker](../scripts/fresh_h20_integer_inner_gap_20261014_review.py)
reconstructed all 651 archived small-power recipes using Gaussian inverses
of the three Cartan blocks, independently of the producer's scaled inverse
formula. Its [generation and saved replay](fresh_h20_integer_inner_gap_20261014_review.json)
passed. Completeness of the omitted-recipe filtering is source-reviewed;
this is not a second enumeration of all factor compositions.

The family is every radius-two cap whose nonnegative integral seed gaps
sum to `7,7,8`, and which is Coxeter-compatible with the retained
radius-two central cap. It imposes no endpoint-layer or mutual-cap
compatibility condition. The center norm is rational. Therefore its
normalized negative pairing with the central cap is rational, and the
complete finite-angle possibilities are 0 and `1/2`; parallel/dotted
relations have value at least 1. In the producer's scaled notation this
is exactly `norm28 in {224,336}` or `norm28>=448`. A duplicate of the
central cap adds no inequality and is separately harmless.

At `x=3z0/20`, the squared distance from a candidate center is a sum of
three nonnegative factor distances. Hence any candidate with ball power
at most `1/100` must pass each individual factor cutoff `44912`, in the
common denominator 11200. The factor enumeration contains all
compositions, and the full Cartesian product of the surviving factors
is checked including equality. No candidate capable of lowering the
claimed bound is discarded by this reduction.

There are 651 small-power recipes. Exactly two pass the central gate;
they exchange the two A6 factors and both have center norm 8,
center pairing 28 with `z0`, and power `1/100`. Every recipe omitted from
the final list has greater power. Thus the uniform lower bound is exact.

The lift `E=(3z0/20,881/400,1)` has square `1/200` and strictly positive
seed gaps at least `17/20`. The central power is `41/100`. Every
vertex-centered cap of radius at most 11 also has power at least
`2061/100`: all old vertices have norm 14, so reverse Euclidean triangle
inequality gives distance at least `14-21/10=119/10` from the probe.
All these power bounds exceed the chosen Lorentz square. The resulting
strict spacelike escape therefore survives the whole stated union,
independently of a particular endpoint layer or mutual compatibility.

This excludes finite-volume completion by finite selections from this
dictionary. It does not exclude other radii, nonintegral seed gaps,
nonvertex endpoint centers, or replacing the central cap.
