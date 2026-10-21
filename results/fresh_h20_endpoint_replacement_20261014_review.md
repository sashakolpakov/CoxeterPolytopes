# Review of the continuous endpoint/face replacement gap

**PASS.** The [source](../scripts/fresh_h20_endpoint_replacement_20261014.py)
and [report](fresh_h20_endpoint_replacement_20261014.md) have no identified
proof or scope gap. A fresh primary `--verify` passed. The additional
[small Fraction checker](../scripts/fresh_h20_endpoint_replacement_20261014_review.py)
and its [result](fresh_h20_endpoint_replacement_20261014_review.json)
independently reconstruct the escape and its 22 seed inequalities, the
radial minimum attainer, and the exact triangle-power cancellation. This
is source/formula review and targeted independent arithmetic, not a second
full family enumeration.

For a radial center `gamma*z` with `gamma<1`, the nineteen positive
retained seed gaps are `1-gamma`. A Coxeter nonright cosine is at least
`1/2`, so the cap radius is at most `2(1-gamma)`. The ball-power lower
bound is the quadratic

`|x|²-4 - 2gamma(x·z-4) + 192gamma²`.

Minimizing over all real `gamma` is a safe relaxation. The exceptional
`gamma=1` endpoint centers are covered by the separately bounded endpoint
radius family; larger `gamma` violates retained seed acuteness. The report
does not extend the endpoint radius bound to arbitrary endpoint caps.

For the shifted triangle, write
`w=z1+r*d`, with `d=(z0-2z1+z2)/n²`. Direct independent reconstruction gives
`d²=4/n²`. The radius square is `4r²/n²-r`, so the quadratic coefficient
in `|x-w|²-R²` vanishes exactly. The remaining power is affine in `r`;
checking both interval endpoints proves the whole closed interval
`16<=r<=49/3`, independently for every cap. The midpoint-cap power is
also affine. The zero-radius A7 triangle at the lower endpoint can be
discarded and creates no missing positive-radius case.

The positive metric root is correctly isolated: its quadratic leading
and linear coefficients are negative, its constant is positive, and its
signs at 16 and `49/3` are opposite. The exact local contact identities
are valid; the much larger numerical all-pair metric scout remains
explicitly separate and is not certified by this review.

At the archived point, all six uniform ball-power bounds exceed 1 and
all vertical seed gaps are at least `1/100`. Thus the displayed Lorentz
lift has square 1, future orientation, and strictly negative pairings
with the entire relaxed dictionary. Any finite subcollection also has
a strict timelike point: choose old spatial coordinate zero and a
sufficiently large upper-half-space height. The strict spacelike escape
therefore proves infinite volume for every finite compatible selection,
including deletion of optional walls. It does not exclude mixed-factor
face caps or differently positioned caps outside the displayed formulas.
