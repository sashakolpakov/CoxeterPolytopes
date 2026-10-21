# Independent review of the arbitrary-height H18 lift obstruction

**PASS.** The stated fixed-source-direction family cannot have finite volume in any target dimension `18+k`, `k>=2`. The restriction to this family is essential: new projected wall directions and deletion of the retained long anchors are outside the theorem.

The [separate checker](../scripts/fresh_h18_general_height_lift_20261014_review.py) and its [certificate](fresh_h18_general_height_lift_20261014_review.json) bind the frozen primary source/data. They import no producer or project arithmetic. Sparse squarefree-radical Fraction arithmetic with rational `isqrt` sign bounds independently rechecks all 40 center pairings, all 120 null-pair wall pairings, the six fixed anchors, and the escape-square coefficients. The source Gram inertia and the original chamber's complete geometric provenance remain inherited dependencies; this is not a second original-source census.

## The all-real reduction

For each of `L18,L25,L26`, a retained norm-two anchor has pairing `-1`. A normalized lift `r+v` therefore has squared cosine `1/(2(1+|v|^2))`. It is positive and at most `1/2`; the only Coxeter possibilities are labels three and four, giving `|v|^2=1` and zero respectively. Neither a right angle nor a divergent pair can occur at that anchor.

For two distinct lifts of the same short root, acuteness requires `1+v.w<=0`. Thus a zero height cannot coexist with another height, and distinct unit heights must be antipodal. Together with the assumed nonempty height set and zero in its convex hull, the only sets are `{0}` and `{u,-u}` with `|u|=1`. Duplicate identical walls may simply be removed.

For two nonzero core height sets, the two cross cosines are `(1+d)/2` and `(1-d)/2`, where `d` is the dot product of their unit directions. Both lie in `[0,1]` and sum to one. A positive finite Coxeter cosine is at least `1/2`; consequently either both are `1/2`, or they are zero and one. The lines are orthogonal or equal. This covers every integer finite label, and the endpoint one is the parallel case.

In two added dimensions, three such lines cannot be pairwise distinct and orthogonal. A zero set or a repeated line therefore yields a pair whose height span has dimension at most one. For `k>=3`, every pair has height span at most two. In either case a nonzero perpendicular added direction exists.

## Strict global escape

For the selected pair, the independent data verify that `delta=r_i+r_j` is future null, pairs to zero with its two short roots, is strictly negative on all other short roots, and is nonpositive on every retained long root. The old center `h` is timelike and strict. Set `tau=-1/(4 delta.h)>0` and choose a unit added vector `u` perpendicular to both selected height sets. Then

```text
X = delta + epsilon*u + epsilon^2*tau*h
X^2 = epsilon^2/2 + epsilon^4*tau^2*h^2 > 0
```

for sufficiently small positive `epsilon`. All walls with a negative constant pairing remain strict because the entire height bank is finite. The selected short lifts and the zero-pairing long walls become strict at order `epsilon^2`; their linear terms vanish. Also `X.h<0`. This supplies a strict spacelike point in a chamber already containing the strict timelike center `h`, and hence an open patch on the future ideal boundary and infinite volume.

The polynomial positivity bound in the certificate controls only `X^2`; the required bound for every individual wall can depend on that finite bank. The theorem correctly uses existence of a sufficiently small parameter and does not claim one universal numerical `epsilon` for arbitrarily large finite heights.

Generation and saved `--verify` both pass. No mathematical or scope issue was found.
