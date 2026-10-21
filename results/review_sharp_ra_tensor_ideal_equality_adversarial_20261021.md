# Adversarial mathematical review of the sharp right-angled paper

21 October 2026. **PASS at the explicitly stated geometric and published-input boundary. No mathematical defect or missing internal proof found.** This is a full manual mathematical source review, not a new checker execution, PDF build, external referee report or end-to-end Lean verification. The reviewer previously inspected some published inputs and verification interfaces; all calculations below were rederived rather than inferred from earlier PASS labels. No manuscript, guide or proof source was changed.

Reviewed source: `manuscripts/right_angled_sharp_dimension_bounds_publication.tex`, SHA-256 `0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722` (382 lines). References below are to that source.

## Projection, actual faces and strictness

Lines 93–105 and 118–164 supply the necessary actual-face bridge. With (q^2=-1), (n^2=1), (q\cdot n=-t), the vector (q+tn) has square (-1-t^2), product zero with (n), and negative product with future (q). For every other actual wall (m), its product is (q\cdot m+t(n\cdot m)<0). Thus the normalized point is future, lies on the actual facet, and is strict against every other wall. This uses the stated **all-wall** nonpositivity, not merely orthogonality of ridge neighbors. Intrinsic facet normals are the actual ridge neighbors of (n); their orthogonality to (n) preserves unit scaling and gives exactly (t'_i=t_i/\sqrt{1+t^2}).

The monotone update (f(t)=t/\sqrt{1+t^2}) gives (A_j\le A_0/\sqrt{1+jA_0^2}<1/\sqrt j). Its strict final inequality follows from finite (A_0), including ties in the maximizing choice. There is no assumed limiting position. In the radial split (lines 197–211), all-wall acuteness gives (u_i\cdot u_j\le f(t_i)f(t_j)). Coincident directions would instead have normal product (\sqrt{(1+t_i^2)(1+t_j^2)}-t_it_j>0), so facet indexing is injective.

The face-descent assertions are valid in the dimensions actually used. In the Klein closure a face has vertices among the original finite/ideal vertices; simplex triangulation gives finite volume. At a finite relative-interior point of an intrinsic ridge, the right-angled local chamber is an orthant and restriction preserves right angles. Every face used has dimension at least four; no ideal polygon-angle inheritance is being assumed. Higher-dimensional exclusions use actual face flags to dimensions nine or five, not inheritance of the theorem by arbitrary Gram subdiagrams.

## Both tensor contractions and the strict 26 bound

For unit (u,v\in\mathbb R^6), put (x=u\cdot v), (C(u)=u^{\otimes3}), and use the paper's (A,S,T) (lines 228–247). The complete contractions are

\[
\begin{aligned}
\langle u\otimes u,v\otimes v\rangle&=x^2,&
\langle u\otimes u,I\rangle&=1,& \langle I,I\rangle&=6,\\
\langle C(u),C(v)\rangle&=x^3,&
\langle C(u),S(v)\rangle&=3x,&
\langle S(u),C(v)\rangle&=3x.
\end{aligned}
\]

Of the nine terms of (\langle S(u),S(v)\rangle), the three matching-slot terms equal (6x) each and the other six equal (x) each. Hence that contraction is (24x). Consequently

\[
\langle A(u),A(v)\rangle=x^2-\tfrac16,
\qquad
\langle T(u),T(v)\rangle=x^3-\tfrac{6x}{8}+\tfrac{24x}{64}
=x^3-\tfrac38x.
\]

Both normalizations (1/6) and (1/8) are correct. Substitution gives exactly the displayed polynomial identity and, after finite double summation, the three nonnegative Euclidean squared norms in (SOS), lines 248–260. No positive-kernel theorem or supplied cardinality conclusion is needed.

Here (F(1)=27/16) and (F(x)=(x+1/2)^2(x-1/4)\le0) when (x<1/4). Thus (N^2\le27N). For (N=27), equality of the total bounds forces every nonpositive off-diagonal summand to be zero. Strict (x<1/4) removes the other root, leaving every distinct-pair product equal to (-1/2). Any four selected vectors would then have squared sum (4+12(-1/2)=-2). This proves (N\le26), including the strict equality exclusion. The manuscript does **not** extend this conclusion to the closed condition (x\le1/4).

Three projections from dimension nine give a six-face with (t_i<1/\sqrt3); therefore (f(t_i)<1/2) and every direction product is strictly below (1/4). Lines 277–291 correctly contradict the universal six-dimensional 27-facet input. Compactness, arithmeticity and ideality are not used in this branch.

## Ideal 24-cell model and the metric equality case

The 24 normals ((1,r)), (r=\pm e_i\pm e_j), in lines 295–318 have Lorentz square one. This is the correct ideal 24-cell model: in Klein coordinates its inequalities are (\lvert y_i\rvert+\lvert y_j\rvert\le1), and its vertices are the eight (\pm e_i) and sixteen ((\pm1/2,\pm1/2,\pm1/2,\pm1/2)), all on the ideal sphere. The facets are regular ideal octahedra. Distinct spatial root products are (-2,-1,0,1), giving normal products (-3,-2,-1,0); finite ridges are the orthogonal pairs. Unit-normal scaling is therefore consistent with the projection lemma.

For each spatial opposite pair ((1,r),(1,-r)), a future unit point (q=(a,w)) has heights (a-r\cdot w) and (a+r\cdot w). If both are at most one then (a\le1), while (a^2-\lVert w\rVert^2=1), (a>0), gives (a\ge1). Thus (a=1,w=0), and **all** heights equal one. The pair is spatially opposite, not a normal and its Lorentz negative. The calculation works even without an interior-point premise and so certainly covers the projected interior point.

A projection from an entirely ideal five-polytope gives an entirely ideal four-facet with all heights (<1), hence direction products (<1/2). Musin bounds its facet count above by 24; Kolpakov bounds it below by 24 and identifies the metric equality case. Isometry then transports its actual outward **unit** normals to the displayed model. The opposite-pair contradiction is valid. Neither combinatorial 24-cell recognition alone nor uniqueness of a kissing configuration would suffice; the paper explicitly uses neither replacement (lines 320–336).

## Precise external inputs and sharpness

| Input | Primary-source check and exact scope |
|---|---|
| 27 facets in dimension six | [Alexandrov, §4 and §6](https://arxiv.org/html/2303.09533v1) explicitly records the six-dimensional finite-volume bound and attributes it to Dufour Lemma 1/Proposition 4. In §6 the noncompact statement follows from the five-dimensional facets-plus-cusps lower bound 26 and (f(P^6)\ge1+\nu(P^5)). Compact right-angled dimensions (>4) are excluded in [Potyagailo–Vinberg, printed p. 65](https://ems.press/content/serial-article-files/43032). Thus the all-vertex-types application is supported. Dufour's paywalled original was not directly reread; the manuscript accurately discloses the Alexandrov route. |
| Spherical (\pi/3) code in (S^3) | [Musin, printed pp. 2, 4](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf) identifies kissing configurations with such spherical codes and proves (k(4)=24). The closed pair-product bound is valid and covers the paper's strict input. No uniqueness theorem is imported. |
| Ideal four-dimensional facets and equality | [Kolpakov, Theorem 4, printed p. 1945](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf) states the **facet minimum** and isometric equality characterization used here. This is not inferred from the separate minimum-volume theorem. |
| Attainment in dimension eight | [Potyagailo–Vinberg, §3, Proposition 3.1, printed pp. 66–67](https://ems.press/content/serial-article-files/43032) constructs finite-volume right-angled (P^n), (3\le n\le8), and lists the eight-dimensional example. The ideal right-angled 24-cell attains dimension four. |

The standard acute-polytope all-wall theorem and face foundations remain explicit published geometric inputs (lines 93–105), not formal consequences of the tensor calculation. The two pointwise tensor contractions are proved in the manuscript, but remain hypotheses in the currently selected Lean `CoordinateWitness`; this is a formalization boundary, not a missing written proof. No revision of the manuscript is required by this review.
