# Independent replacement-cap escape review

Generation and saved-result replay **PASS**. The [independent checker](../scripts/fresh_h20_e7_replacement_escape_20261014_review.py) imports no producer or project arithmetic. Its [certificate](fresh_h20_e7_replacement_escape_20261014_review.json) binds the two reviewed sources, their data, and its own source by SHA-256.

The checker reconstructs the physical E7 and A11+A2 roots, both dual bases, and the timelike vector from the saved ambient form. All thirteen physical A roots have norm two; the scaled intrinsic coordinate basis is not substituted for them. It decodes all 42 source-index/value-index recipes, reconstructs their norm radicals, checks norm two and their prescribed contacts, and independently checks strict negative pairings with the common timelike center and the saved escape. General squarefree-monomial Fraction arithmetic handles the extra square root of three; outward rational enclosures handle the positive norm square roots without an independence assumption. Every sign resolves at 128 bits.

The selected cap 34 is exactly

\[
r=-W_2-W_4-W_6-\sqrt2L_0-2L_{12}+5T.
\]

It cuts the previous escape and is strictly negative on the new one. The new escape has square `4374521734010421633/1000000000000000000` and negative center pairing. It is strict against all 87 original walls, and in particular against the retained 85. The 42 recipes therefore cannot remove this escape, even in their simultaneous halfspace relaxation. This does not certify that those recipes are mutually Coxeter compatible or actual facets.

The separate [continuous gate](../scripts/fresh_h20_e7_replacement_escape_20261014_continuous.py) also passes independent reconstruction. With

\[
V=-W_2-W_4-W_6-\sqrt2L_0,\qquad
r(x)=V-xL_{12}+\sqrt{(71+x^2)/3}\,T,
\]

the three orthogonal vectors have squares `148/3`, `2/3`, and `-2`, so every real parameter gives norm two. Its escape pairing is

\[
C+ax-B\sqrt{(71+x^2)/3}.
\]

The checker independently reconstructs and proves `C,a,B>0`, `D=B²−3a²>0`, and `71D/3−C²>0`. Lorentzian Cauchy gives

\[
B\sqrt{(71+x^2)/3}-ax\ge\sqrt{71D/3}>C,
\]

which proves a strict negative pairing for **every real x**, including parameters not producing Coxeter caps. A finite coherent selection from this one family, retaining any subset of the old walls and defining a nonempty chamber in the chosen future component, still has an open ideal escape. Other E7 profiles and arbitrary replacements are outside this statement.

Replay:

```sh
python scripts/fresh_h20_e7_replacement_escape_20261014_review.py --verify
```

The next constructive incidence attempt uses a five-cube with four adjacent positive-facet ridge cuts, followed by wedges at four mutually meeting negative facets. Its missing faces should have orders two and three only. This changes both the base and the missing-face structure from the previously excluded fifteen-facet Gale and product-section batches; it is a proposal, not a surviving Coxeter metric.
