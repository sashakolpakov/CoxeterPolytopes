# Review of the fixed-W paired-curve obstruction

**PASS within the stated curve family.** The frozen [source](../scripts/fresh_h20_structured_extension_20261014.py) and [certificate](fresh_h20_structured_extension_20261014.json) support the theorem in the [report](fresh_h20_structured_extension_20261014.md). A fresh producer `--verify` passed. The [separate small checker](../scripts/fresh_h20_structured_extension_20261014_review.py) and [review certificate](fresh_h20_structured_extension_20261014_review.json) independently check the four rational frames, all 120 predecessor center/escape signs and the analytic identities below. This is not a second implementation of the full interval census.

The two-sheet parameterization is complete for unit spacelike vectors in the specified Lorentz plane. Adding and subtracting the fixed orthogonal unit vector W gives two norm-two, mutually orthogonal caps. Requiring **both** caps makes their old-wall acuteness equivalent to the stated quadratic inequality containing `2|W·r|t`. The negative-sheet witness excludes every positive parameter on that sheet. On the positive sheet, concavity and strict positivity at both interval endpoints exclude the whole closed interval; no singleton or boundary value is dropped.

The extra pairing used for supports 0 and 2 is strictly decreasing for `t>1`, because both its E and T coefficients are negative. Its two endpoint bounds therefore place the entire remaining interval in `(-1,0)`. Both root norms are two, so the normalized cosine lies in `(0,1/2)`, which excludes every finite Coxeter label as well as parallel and dotted pairs. This normalization is correct. The archived acute but non-Coxeter cutter also demonstrates why acuteness alone is insufficient.

The strengthened source explicitly checks a common timelike center: `H²=-4296`, `E·H<0`, `T·H<0`, `W·H=0`, and every height reflection fixes H. All necessary lower parameters exceed one. Hence every cap on every retained parameter tail contains H strictly. The common spacelike witness also pairs strictly negatively with all sixty predecessors and all admissible paired caps. Thus any finite selection remains a nonempty infinite-volume chamber; deleting previous optional paired caps cannot remove that witness.

The cross-family argument is valid. Differentiating a pairing with respect to `log(t)` gives the two brackets recorded by the first four derivative gates; the other four control the second parameter. The absolute values cover both choices of cap sign, while the negative scalar coefficients and `cosh>=1`, `sinh>=0` keep both derivatives negative on the entire quadrant. The strictly dotted lower-corner values therefore persist for all larger parameters. Within one family, the positive same-sign pairing at distinct parameters prevents two complete pairs from coexisting acutely.

The height-reflection permutation check is sound with the new explicit assertion that every predecessor radicand is two. The separate checker confirms rank 18 of each rational support, exact frame orthogonality and positive/negative squares, using no producer arithmetic. It also checks the hyperbolic normalization, acuteness polynomial and both cross-pair derivative formulas symbolically.

The negative conclusion requires fixed W, one of the four saved supports or its height reflections, both members of each pair, and retention of all sixty predecessor facets. It does not exclude individual unpaired caps, rotations within the support complement, changed supports, or replacement of predecessor walls. It is not a global obstruction to the H20 construction.

```sh
python scripts/fresh_h20_structured_extension_20261014_review.py --verify
```
