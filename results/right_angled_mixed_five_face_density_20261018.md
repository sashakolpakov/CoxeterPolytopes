# Necessary five-face distributions in right-angled dimensions nine to twelve

18 October 2026. **Derived universal constraints; two full independent internal reviews completed.** Every entirely ideal right-angled hyperbolic five-polytope has at least 56 facets. Consequently, in a finite-volume right-angled polytope of dimension 11 or 12, fewer than one tenth of its five-faces can be entirely ideal. This is a forbidden boundary distribution, not a dimension exclusion, and no optimality or literature priority is asserted.

## An ideal five-polytope needs at least 56 facets

The frozen [five-face carrier theorem](right_angled_five_face_carrier_20261018.md) supplies a four-facet Q with positive $\kappa(Q)=3f_0(Q)-4f_3(Q)+24$. The [second-value equality theorem](ideal_ra4_second_volume_equality_20261018.md) implies $\Delta(Q)=f_0(Q)-f_3(Q)\ge5$. Kolpakov's published facet-minimum theorem gives $f_3(Q)\ge24$, with equality only for the regular 24-cell; positive kappa excludes that equality. Hence

\[
 f_3(Q)\ge25,\qquad f_0(Q)\ge30,\qquad f_3(Q)+f_0(Q)\ge55.
\]

In any ideal right-angled five-polytope P containing Q as a facet, its $f_3(Q)$ ridge-neighbor facets are distinct. At every cusp of Q the rectangular four-box link has a unique facet opposite Q. The hyperplane of such an opposite is asymptotic to Q's hyperplane at exactly that cusp. It cannot be opposite at another cusp, and cannot equal a ridge neighbor, which is orthogonal to Q's hyperplane. Thus Q, all ridge neighbors and all cusp-opposites are distinct. This proves

\[
 \boxed{f_4(P)\ge1+f_3(Q)+f_0(Q)\ge56.} \tag{1}
\]

The positive-kappa sixteen-cell construction shows that such carrier geometry is not excluded in dimension four. Its one-layer attachment is a separate problem; (1) does not construct a five-polytope.

## Explicit published premises and access limits

We use the Nikulin–Khovanskii strict bound for a polytope simple at edges:

\[
 \frac{1}{f_k(P)}\sum_{F\in\Omega_k(P)}f_l(F)
 <\binom{n-l}{n-k}
 \frac{\binom{\lfloor n/2\rfloor}{l}+\binom{\lceil n/2\rceil}{l}}
 {\binom{\lfloor n/2\rfloor}{k}+\binom{\lceil n/2\rceil}{k}}.
 \tag{2}
\]

Only $k=5$, $l=1,4$, and $n=9,10,11,12$ are used, within the range $0\le l<k\le\lceil n/2\rceil$. Finite-volume right-angled polytopes are simple at edges, and every positive-dimensional face is itself finite-volume and right-angled. Thus the hypotheses apply to the ambient polytope and to the five-faces used below. The formula and these right-angled facts are available in [Alexandrov's 2023 primary preprint, Sections 2 and 6, Proposition 6.1](https://arxiv.org/html/2303.09533v1#S6).

We also use the published Dufour lower bound of 16 facets for a finite-volume right-angled five-polytope: for a noncompact one it follows from $f_4\ge1+\nu_4$ and $\nu_4=15$, explicitly reported in [Egorov's 2026 primary preprint, Lemma 2.3](https://arxiv.org/html/2603.29925v1#S2). There are no compact right-angled five-polytopes, by the Potyagailo–Vinberg theorem stated as [Alexandrov's Theorem 1.1](https://arxiv.org/html/2303.09533v1#S1.SS1). These are external published inputs, not new theorems of the exact script.

The [original Dufour article](https://link.springer.com/article/10.1007/s10711-009-9454-2) was accessible only as a subscription preview on 18 October 2026. Its full proof and its equality mechanisms were therefore **not** independently read here. Alexandrov's accessible Section 5 gives an alternative thirteen-dimensional exclusion using weighted two-face incidences; that argument does not itself exclude dimension 12 or classify its potential extremal boundaries. The present calculation does not claim to complete either argument. The [Kolpakov 2012 facet theorem](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf) and the internally proved carrier/equality inputs remain separate.

## Few five-faces can be entirely ideal

Let $F=f_5(P)$ and let I count the entirely ideal five-faces. Every five-face has at least 16 facets; the ideal ones have at least 56 by (1). If U is the strict upper bound in (2) for $k=5,l=4$, then

\[
 16F+40I\le\sum_{E\in\Omega_5(P)}f_4(E)<UF.
\]

Therefore $I/F<(U-16)/40$. Exact substitution gives:

|Ambient dimension|Strict average facet bound U|Strict ideal five-face fraction upper bound|Strict nonideal five-face fraction lower bound|
|---:|---:|---:|---:|
|9|30|7/20|13/20|
|10|30|7/20|13/20|
|11|20|1/10|9/10|
|12|20|1/10|9/10|

In particular a proposed dimension-11 or dimension-12 boundary with at least one tenth of its five-faces entirely ideal is impossible. This is a coupled incidence restriction and remains valid although kappa is unbounded on individual four-polytopes.

## A consequent relation between finite and ideal vertices

Let $v_f,v_i$ count the finite and ideal vertices. Every nonideal five-face contains a finite vertex, and each finite vertex belongs to $\binom n5$ five-faces. Write c for the strict nonideal fraction lower bound in the table. Then

\[
 \binom n5v_f\ge F-I>cF. \tag{3}
\]

For an optional further consequence, use (2) with $k=5,l=1$. An edge belongs to $\binom{n-1}{4}$ five-faces. The global edge identity is

\[
 2f_1(P)=nv_f+2^{n-1}v_i;
\]

the ideal-vertex degree is the **number of vertices** of its $(n-1)$-box link, not its $2(n-1)$ facets. If the edge average bound is $U_e$, this gives

\[
 F>\frac{\binom{n-1}{4}}{U_e}
 \left(\frac n2v_f+2^{n-2}v_i\right)=a v_f+b v_i.
\]

Combining with (3) gives the strict necessary relation

\[
 v_f>\frac{cb}{\binom n5-ca}\,v_i. \tag{4}
\]

The saved exact fractions are:

|n|Strict lower bound for $F$|Strict lower bound for $v_f/v_i$|
|---:|---:|---:|
|9|$v_f/2+128v_i/9$|$3328/45243$|
|10|$v_f+256v_i/5$|$3328/25135$|
|11|$7v_f/2+3584v_i/11$|$3072/4807$|
|12|$6v_f+1024v_i$|$512/437$|

The algebraic form (4) is used even if $v_i=0$; the ratio notation is only shorthand when $v_i>0$. These are derived coupled inequalities, not claims to improve every previously published absolute vertex bound. For dimension 12, for example, $437v_f>512v_i$ is a necessary relation, not a contradiction.

## Exact scope

The [producer](../scripts/fresh_right_angled_mixed_five_face_density_20261018.py) binds the two internal theorem reports and checks all rational substitutions, signs, four density rows and resulting vertex ratios. It does not reprove Nikulin–Khovanskii, Dufour or the compact-dimension theorem, infer universal geometry from a bounded census, or perform an optimizer run. The geometric counting arguments are written above. Generation and a separate-process read-only replay both completed with exit 0/PASS in terminal chunk `32cd3c`.

The order16 agent completed a full mathematical read of this producer and report (`fe68e2`/`21215c`), including a reread of the frozen carrier/equality consequence. The actual-facet injection, NK substitutions, simple-edge and finite-vertex incidences, ideal box degree and strict rational elimination were checked by proof/source review. That reviewer did not execute the arithmetic, independently reprove the frozen equality theorem, or independently open/reprove Dufour/NK.

The root agent independently completed a full mathematical source/report read, conditional on the explicitly listed published inputs. That review also checked the NK formula and range in Alexandrov Proposition6.1/Egorov Theorem2.1, Dufour's recurrence as stated in Alexandrov Proposition6.5, and the constant nu4=15 in Egorov Lemma2.3. It did not execute the arithmetic or read the inaccessible original Dufour proof. Both reviews found no correction. The [four-page manuscript](../manuscripts/right_angled_five_dimensional_constraints.tex) and [reproduction/lemma map](../docs/right_angled_five_dimensional_constraints_reproduce.md) keep the actual infinite-volume attachment separate from these universal incidence constraints. No full formal-kernel or external human-review claim is made.
