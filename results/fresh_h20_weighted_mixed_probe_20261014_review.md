# Independent review of the mixed weighted-box covering obstruction

The mathematical argument passes review. It excludes finite-volume completion of the retained weighted-box/central core using any coherent finite combination of **all seed-compatible radius-two caps and all seed-compatible vertex-centered caps of arbitrary radius**. It does not exclude other nonvertex radii or a changed core.

This review concerns the [producer](../scripts/fresh_h20_weighted_mixed_probe_20261014.py), [certificate](fresh_h20_weighted_mixed_probe_20261014.json), and [report](fresh_h20_weighted_mixed_probe_20261014.md). The separate [Fraction checker](../scripts/fresh_h20_weighted_mixed_probe_20261014_review.py) and [review record](fresh_h20_weighted_mixed_probe_20261014_review.json) use no producer or project imports. Independent generation and saved replay both **PASS**. The record binds the reviewed files and prior classification source by SHA-256.

## Complete radius-two coordinate alphabet

The smaller opposite-wall gap is at most the inradius `b<2`; therefore it cannot be parallel or dotted. Since the gap is `2cos(pi/m)`, its label is at most 3, 4 or 5 for inradius `1`, `sqrt(2)` or `phi`, respectively. This leaves the old finite coordinate list, with just two complementary gaps needing an all-label exclusion:

`2sqrt(2)-1` and `2phi-sqrt(2)`.

Both lie strictly between `1.82` and `1.84`. The independent checker proves this from rational isolating intervals for the radicals. Also

`2cos(pi/7) < 2*(1-(3/7)²/2+(3/7)^4/24) < 1.82`,

using `pi>3` and the alternating cosine bound, while `2cos(pi/8)=sqrt(2+sqrt(2))>1.84`. Monotonicity excludes every finite label for either gap; they are also below the parallel/dotted threshold two. Thus no finite-label cutoff remains in the radius-two family.

The prior exact classification of the 88 possible cutters of each half-vertex probe, and its adjacent-probe compatibility rule, are explicit dependencies. Their source was inspected; this review does not repeat the 7,744 recipe-pair checks. Every cutter designates one full golden coordinate. For probes differing only in golden coordinate `j`, compatibility requires exactly one endpoint to designate `j`. The additional zero-unit-coordinate conditions can only strengthen this necessary rule.

## Independent covering arithmetic and graph argument

For the weighted box, `N=20+8phi`. A radius-four ball centered at a box vertex covers a half-vertex probe exactly when

`H=h1+2h2+(1+phi)h5 <= 11/2-phi`.

The independent checker exhausts all 198 weighted Hamming cases. With no changed golden sign there are `176+11=187` probes; with one changed golden sign there are `8*11=88`; all larger golden differences fail. Hence each radius-four ball covers exactly **275** probes. Equality is impossible, so this count has no boundary ambiguity.

Fix the eleven unit/heavy signs. In the resulting golden cube `Q8`, let `A` be all probes covered by the selected radius-two caps. Choose one radius-two cutter for each vertex of `A`. Orient each induced edge toward its unique endpoint designating that edge coordinate. Each vertex receives at most one edge, so `e(A)<=|A|`.

Degree counting also gives `e(A)>=8*(|A|-128)`. Therefore `|A|<=floor(1024/7)=146`, or at most `2048*146=299008` radius-two-covered probes in total. The prior compatible-family bound is 512 radius-four balls, whose total coverage is at most `512*275=140800`. The union leaves at least `524288-299008-140800=84480` probes uncovered; overlaps cannot increase the upper bound. Equivalently, full coverage would require at least `2048*110=225280` radius-four-covered probes, exceeding 140,800.

For clarity, the prior 512 bound follows from a weighted Plotkin bound of four caps per golden sign pattern and the incompatibility of adjacent golden patterns, giving at most 128 represented patterns. It concerns pair-compatible caps and remains valid when other cap families are added.

## Arbitrary vertex radii and the finite-volume implication

At any box vertex an opposite unit wall has distance two. Seed compatibility therefore forces either `R=4` at label 3, or `R<=sqrt(8)` at every other finite label or parallel/dotted contact. A right contact is impossible for this positive gap. Every half-probe has squared distance at least `N/4=5+2phi>8` from every box vertex. Thus all non-radius-four vertex caps miss every probe strictly, regardless of their other contacts.

The half-probes are strictly inside the vertical box and strictly outside the central ball. A radius-two covering recipe cuts its probe strictly, with power `8-5phi<0`; the complete classification has no boundary-only exception. Together with the strict radius-four Hamming comparisons, the counting contradiction gives a probe strictly outside every selected cap. Its neighborhood is an uncovered open ideal-boundary patch. This proves infinite volume for every coherent finite selection from the stated union.

The independent arithmetic checker verifies the Hamming counts, graph integer bounds, radical enclosures, strict margins and radius reduction. The proof review supplies the orientation/counting argument and its geometric implication. This is not a second exhaustive classification of all radius-two cutters or an exclusion of arbitrary H20 caps.

```sh
python scripts/fresh_h20_weighted_mixed_probe_20261014_review.py --verify
```
