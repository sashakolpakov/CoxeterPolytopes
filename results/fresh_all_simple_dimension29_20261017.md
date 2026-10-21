# Dimension 29 for finite-volume Coxeter polytopes simple at every vertex

**Complete internally checked proof composition; generation and fresh saved replay PASS.** The conclusion is:

> A finite-sided, finite-volume hyperbolic Coxeter polytope which is simple at every finite and ideal vertex has dimension at most29.

There is no arithmetic hypothesis. The simplicity requirement at cusps is essential to this proof. This is not an unrestricted finite-volume bound29, a claim of a new example, or a novelty/current-best claim. The unrestricted bound163 and the independent compact-facet ambient bound127 are unchanged.

The [exact composer](../scripts/fresh_all_simple_dimension29_20261017.py) and [saved certificate](fresh_all_simple_dimension29_20261017.json) bind the separate geometric proofs and verify their common hypotheses, the exhaustive case partition, and the rational final inequality. They do not infer a geometric theorem from a numerical sample or from a certificate status alone. The [limited primary-source scope check](fresh_simple_dimension29_literature_scope_20261017.md) found the cited dimension29 statements under compact hypotheses; it is not an exhaustive priority search and makes no novelty determination.

## One rule at finite and ideal vertices

It suffices to exclude n>=30, so every following ingredient's minimum dimension, at most18, is satisfied. At finite vertices retain Vinberg's original nonnegative corner weights. At ideal vertices use the [affine rule](fresh_simple_cusp_vinberg_weight_frontier_20261017.md): weight1 for marked roots at graph distance at most3 or in a permitted end segment of order at most6, and weight0 otherwise. Fork-end segments respect interchange of the two fork leaves. A cycle has no end segment.

At a finite vertex the total weight is at most3n. The finite connected-component classification and counts are the published argument in Vinberg §6.2, also checked in the [finite-component package](fresh_compact_dimension_bound_20261015.md). All nonzero payments stay inside one component. The only unbounded-rank components are classical; at rank r>=10, their r−1 adjacent pairs and at most20 further end-segment pairs give at most r+19<=3r. Rank2 includes every finite dihedral label, without a label cutoff.

An ideal vertex in this class has a connected affine diagram on n nodes. For n>=18 it has one of the classical cycle, path, single-fork or double-fork shapes. The ideal rule has total exactly3n: the short-distance contributions are respectively3n,3n−6,3n−5,3n−4, and the additional end contributions are0,6,5,4. These are proved all-rank formulas; the saved graph checks confirm them rather than supply a rank cutoff. Thus the total corner payment over all vertices is at most3n f0.

## Every two-face receives its required payment

The [global affine-exclusion lemma](fresh_simple_cusp_affine_exclusion_20261017.md) says that every connected affine subdiagram has exactly n nodes. In particular the minimal indefinite cores relevant to the triangle proof are Lannér diagrams of order at most5. The weight interface transfers the literal cases of Vinberg §6.4 and proves total payment at least2 on every triangle, including those with ideal vertices.

Quadrilaterals have three exhaustive branches:

| Corner types | Bound input | Minimum dimension |
|---|---|---:|
| All four vertices finite | [Vinberg finite-corner transfer](fresh_simple_finite_quad_vinberg_transfer_20261017.md) |18|
| At least one affine-cycle vertex | [Complete affine A proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md), [independent review](review_simple_cusp_affine_a_quadrilateral_20261017.md) |14|
| At least one classical affine-tree vertex | [Complete affine-tree proof](fresh_simple_cusp_affine_tree_completion_20261017.md), [terminal/E8 audit](fresh_simple_cusp_affine_tree_terminal_audit_20261017.md) |18|

Each branch gives total payment at least1. If a quadrilateral has any ideal vertex, that vertex is a connected classical affine n-node diagram, so the second or third branch applies. Otherwise its four vertices are finite and its face is compact, so the first applies. Affine exceptional diagrams have order at most9 and do not leave another branch in the dimensions being excluded.

The finite-corner transfer and its [separate source review](fresh_simple_finite_quad_vinberg_transfer_20261017_review.md) check the hypotheses in Vinberg §6.5–11. They do not apply the compact ambient theorem to a noncompact polytope. In particular an affine diagram of order n omitting a white supporting wall is not dismissed: indefinite carriers of order at most n are affine-free for a valid size reason, pieces orthogonal to an actual indefinite diagram are positive definite by Lorentz geometry, and the full quadrilateral dependence proves the determinant polynomial used in §6.10 is nonzero. Published finite classification and local determinant calculations remain explicit proof inputs.

For a k-gon with k>=5, nonnegative weights already give payment at least0>=5−k. Therefore, writing

\[
D_5=\sum_{F\text{ a two-face}}(5-k_F),
\]

the complete face payment statement is

\[
3n f_0\ \ge\ \text{total corner payment}\ \ge\ D_5.
\]

The earlier frozen weight interface correctly says its quadrilateral obligation was open at that stage. The individual affine reports similarly leave the other branches open. The present composition binds the separate completed proofs; it does not change those historical reports.

## The strict combinatorial inequality

The projective closure of a finite-volume hyperbolic polytope is the convex hull of its finite and ideal vertices, inside the closed Klein ball. Here that closure is a simple convex polytope, because simplicity is assumed at all its vertices. Its dual simplicial h-vector satisfies h_i>=0, h_i=h_(n−i), and h_0=h_n=1. The standard face identities give

\[
f_0=\sum_i h_i,\qquad f_2=\sum_i\binom{i}{2}h_i,
\qquad \sum_F k_F=\binom n2 f_0.
\]

In particular D5=5f2−binom(n,2)f0. Set

\[
\beta_n=\begin{cases}
n(n-6)/8,&n\text{ even},\\
(n-1)(n-5)/8,&n\text{ odd}.
\end{cases}
\]

After pairing h_i with h_(n−i), the exact slack is

\[
D_5-\beta_n f_0=
\begin{cases}
\frac52\sum_i(i-n/2)^2h_i,&n\text{ even},\\[2mm]
\frac52\sum_i((i-n/2)^2-\frac14)h_i,&n\text{ odd}.
\end{cases}
\]

All coefficients are nonnegative. The coefficients at i=0,n are strictly positive and both endpoint h entries equal1. Hence D5>β_n f0, including the first threshold dimension where the scalar budgets are equal. This is the strict simple-polytope inequality used by Vinberg's weighted-face argument, written as an exact h-vector identity rather than a rounded average.

For every integer t>=0,

\[
8(\beta_{30+2t}-3(30+2t))=4t^2+60t\ge0,
\]
\[
8(\beta_{31+2t}-3(31+2t))=4t^2+64t+36>0.
\]

Thus every n>=30 would satisfy the contradiction

\[
3n f_0\ge D_5>\beta_n f_0\ge3n f_0.
\]

The dimension is at most29. Dimension29 itself is not excluded by this calculation: β_29=84<87=3·29.

## Verification and scope

The composer checks the finite-component costs, all348 saved affine cost rows, the h-vector coefficient identities, and the parity polynomial identities coefficientwise. Finite-rank examples in the output are verification aids; the displayed formulas and nonnegative polynomial coefficients supply the all-rank conclusion. It also checks all recursively named source manifests, without claiming to rerun every inherited census.

Generation and a fresh saved replay both terminate with PASS. The recursive integrity check covers13 manifests,85 hash references and42 distinct inherited files; the composer also binds29 direct source/report inputs, including both transfer reviews' recorded scope and the limited literature check. Syntax and local report links pass. Separate root and frontier source/mathematical reads find no remaining gap in the finite-corner transfer; root independently rederived the h-vector identity, strictness and both threshold polynomials. The affine cases have their separately linked exact reviews. This package is final and frozen.

The proof uses published finite Coxeter classifications, Vinberg's finite local lemmas, the standard simple-polytope h-vector facts, and the separately reviewed geometric extensions linked above. Internal review is distinct from external refereeing or formal proof-assistant verification. No claim concerning nonsimple ideal vertices follows from this composition.

Primary source for the finite rule and weighted argument: [Vinberg, Trudy MMO47 (1984), §§1 and6](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). The primary theorem there concerns compact polytopes; the mixed-simple conclusion here depends on the new affine-corner and finite-corner transfer proofs, not on broadening the wording of that theorem.
