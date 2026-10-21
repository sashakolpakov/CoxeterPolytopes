# Independent geometric and literature audit of the Prokhorov refinement

**PASS, with a literature-status qualification.** The proposed upper bound 990 follows from the cited published geometric lemmas and the new additive diagram estimate. It requires no arithmetic, compactness, ideal-vertex simplicity, finite-label cutoff, facet count, or source-chamber restriction. It is an explicit refinement of the published 995 argument. This review does **not** establish publication priority or a new strongest known bound: an older primary source announces a stronger bound below 300.

The [small standalone checker](../scripts/dimension_bound_prokhorov_20261015_proof_review.py) and [certificate](dimension_bound_prokhorov_20261015_proof_review.json) bind the [primary construction](dimension_bound_prokhorov_20261015.json) and [independent arithmetic audit](dimension_bound_prokhorov_20261015_review.json). They check the elementary identities below. The geometric audit is a source/proof review, not another computerized proof of the imported topological classification.

## Why the edge estimate applies, including at cusps

[Prokhorov's original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=1484&what=fullteng), printed pp. 408 and 410, identifies the relevant edge diagram as elliptic of rank/order (n-1). Its unordered marked-node pairs determine the three-face angle incidences. Lemma 2 on p. 406 permits any positive edge constant (C); Proposition 5 supplies the good-face lower bound. The proposed argument changes only that constant, retaining every angle weight and all geometric lower bounds.

There is also a direct local explanation. An edge has points in hyperbolic space even when one or both endpoints are ideal. At a relative interior point, its incident reflections act on a positive-definite transverse space of dimension (n-1). Its diagram is spherical. The affine diagrams at ideal endpoints are different objects; they must not replace this edge diagram. Edge simplicity gives the local Boolean interval of containing faces, so deleting an unordered pair of incident facet constraints selects the corresponding three-face. In fact, the upper estimate only needs an injection into node pairs.

Consequently there is no omitted affine edge family, cusp multiplicity, or ordered-pair factor. Distances are ordinary numbers of graph edges, ignoring labels. Disconnected components contribute no pairs. Unbounded finite labels occur in (I_2(m)), whose unlabelled graph is always the two-node path. The complete finite real Coxeter classification has no exceptional connected family of unbounded rank beyond paths and (D_r).

## An all-rank argument without a small-diagram census

Put (W=(2f_7+f_{15})/3), where (f_k) counts unordered node pairs at graph distance at most (k). For a path of rank (r\ge17), (f_k=kr-k(k+1)/2) for (k=7,15). For (D_r), each of these counts is larger by one. Thus the deficits (29r-3W) are respectively 176 and 173.

For any component of rank (1\le r\le16), even the complete-graph bound suffices:

\[
29r-3W\ge29r-3\binom r2=\frac{r(61-3r)}2\ge\frac{13r}2>0.
\]

For a spherical forest of total rank (N\ge27), a large component already supplies deficit at least 173, while all remaining deficits are positive. If every component is small, their combined deficit is at least (13N/2>173). Hence (W\le(29N-173)/3) for every such forest. This weaker threshold than the primary's (N\ge17) is sufficient for the dimension theorem, and avoids depending on any finite exceptional-graph calculation.

For a hypothetical dimension (n\ge991), choose the single positive constant

\[
C=C_n=\frac{29}{3}-\frac{173}{3(n-1)}.
\]

It is uniform over all edges of this polytope. Dependence on the fixed dimension is permitted. The published lemma yields

\[
n<996-\frac{5536}{n-1},\qquad
n^2-997n+6532<0.
\]

Writing (n=991+t) turns the last polynomial into (t^2+985t+586>0), a contradiction for every real (t\ge0). At 990 it equals (-398), so this argument stops at 990. The primary's separate parity calculation agrees; the simpler substitution already proves the conclusion.

## The older 300 announcement must remain visible

The official [ICM 1986 proceedings](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf) contain Nikulin's *Discrete Reflection Groups in Lobachevsky Spaces and Algebraic Surfaces*, printed pp. 654–671. His definition on p. 654 makes CHRG mean finite-covolume reflection groups, including nonarithmetic groups. Printed p. 659 explicitly asserts dimension below 300. Pages 657–659 outline a narrow-face method and work out an affine-(A) cusp example, but do not supply the full numerical treatment of all diagram types needed to verify the final 300 claim from this text alone. The local audit read the original pages, not merely a search snippet.

[Dolgachev's 2008 survey](https://sites.lsa.umich.edu/idolga/wp-content/uploads/sites/1334/2024/08/reflections08.pdf), printed p. 11, explicitly calls that result an announcement and cites Nikulin. More recent papers still quote 995, including the [2026 paper on strongly reversible geodesics](https://www.cambridge.org/core/journals/journal-of-the-institute-of-mathematics-of-jussieu/article/counting-and-equidistribution-of-strongly-reversible-closed-geodesics-in-negative-curvature/176CFA32FE84C00D6F098FF8D3FCA29F). This disagreement in citation practice does not establish that the announcement is false, nor does the announcement alone constitute an independently checked complete proof here.

The accurate status is therefore: **990 is a fully explicit, independently audited consequence of the stated published Prokhorov inputs; a stronger historical primary announcement exists, and priority/current-best status is unresolved.** Recovering the full narrow-face constants is a substantive next task. No external announcement should erase that distinction.

For reproducible provenance, the source above was downloaded as the complete official 974-page first volume. Its PDF pages 756–773 contain printed pages 654–671; the local review excerpt was extracted from precisely those pages, without OCR substitutions. The source URL and page mapping, rather than a transient local path, are the durable reference.

The missing verification target is concrete: derive uniform explicit separation inequalities for every irreducible affine cusp type and every product of affine components, with the possible outgoing edge completions fully covered. Then compare those bounds with the narrow-face distance constraint and treat the alternative where the narrow facet has only finite vertices. The affine-(A) illustration alone does not establish those omitted cases. This is a proposed proof-reconstruction task, not an assertion that the announcement cannot be completed.

## Separate conditional compact audit

The geometric argument in the [compact-component report](dimension_bound_lane_20261015.md) also passes source/proof review. At a finite vertex, inverse positive-definite Coxeter Gram blocks count nonright intrinsic two-face angles, giving the stated average component-pair inequality and the conditional rank bound (d\le4r). This counts all pairs within each component, not only graph edges.

Its stronger quadrilateral charge requires the full **flag/3-free** premise: every empty facet family contains a disjoint pair. Merely having no triangular two-faces is insufficient. Under that premise, opposite quadrilateral side facets are globally disjoint. If all adjacent pairs were right, the four normals would give two negative directions. Thus one direct spherical-diagram edge can be charged at a corner, and the corner plus released pair determines the two-face, making the charge injective. This justifies (q_4\le\sum_v(d-c(v))) and the conditional bounds 11 or 10 stated there. The proof uses [Alexandrov's definition and Lemma 4.3](https://arxiv.org/html/2203.07248v2#S4). These are precisely scoped consequences, with no unrestricted compact-dimension or novelty claim.

## Scope of this review

The proof imports the published good/bad three-face and cusp analysis without attempting to re-prove every quasi-Lannér case. It checks that the unchanged construction meets the exact hypotheses needed for the improved edge constant. The separate arithmetic reviewer recomputes graph distances and scalar algebra independently. No floating-point evidence or finite search cutoff enters the dimension conclusion. The checker binds source/data only, allowing the primary report to add review links without invalidating the audit.

```sh
python scripts/dimension_bound_prokhorov_20261015_proof_review.py --verify
```
