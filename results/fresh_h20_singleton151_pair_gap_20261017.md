# One singleton151 pair has no Coxeter angle at any label

This exact calculation handles one pair left outside the fixed-profile endpoint theorem. The first cap has only a singleton151 contact of label m≥3 at the old edge159/160. The second is fixed profile37 at edge160/161, with contacts151 of label4 and161 of label3. Both are right to ambient wall158, and both use the coherent future norm-completion branch.

For every integer m≥3, their mutual product is either obtuse or strictly between consecutive finite Coxeter angles. This excludes this one pair, not the entire singleton151 family or arbitrary H3 closures. The [source](../scripts/fresh_h20_singleton151_pair_gap_20261017.py) and [exact certificate](fresh_h20_singleton151_pair_gap_20261017.json) pass generation and fresh saved replay. This is a primary certificate using the frozen exact162 interval engine, with the rational cosine routine from the [fixed-profile gate](fresh_h20_e8_endpoint_belt_20261017.md).

Write d=2cos(π/m). Direct projection onto the two old finite rank20 corners, using the literal21-root Gram, gives the first root as

\[
 r(d)=d p+\sqrt{8d^2-2}\,V,
 \qquad p^2=8,\quad p\cdot V=0,\quad V^2=-1.
\]

If s is the second cap, its exact norm is2 and

\[
 f(d)=-\frac{r(d)\cdot s}{2}
     =-4\sqrt2\,d+\beta\sqrt{8d^2-2},
 \qquad \beta=-\frac{V\cdot s}{2}>0.
\]

The certificate saves the algebraic expression for β, reconstructed from the Gram. It checks f′(2)>0. Since

\[
 f''(d)=-\frac{16\beta}{(8d^2-2)^{3/2}}<0,
\]

f is strictly increasing throughout1≤d≤2. At m=3 its value is negative, so the two walls are obtuse. At m=4,5 it lies strictly between the label2 and3 values; at m=6,7,8,9 it lies strictly between labels3 and4. Finally the exact interval checks give

\[
 \cos(\pi/4)<f(2\cos(\pi/10))
 \le f(d)\le f(2)<\cos(\pi/5)
\]

for every m≥10. The limiting value is approximately0.79948145046, with a gap of approximately0.00953554 below cos(π/5). The proof uses the saved rational bounds, not these decimal summaries. It therefore covers the unbounded label tail without sampling or extrapolation.
