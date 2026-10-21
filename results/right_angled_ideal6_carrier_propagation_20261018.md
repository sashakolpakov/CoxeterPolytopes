# Ideal six-dimensional size and mixed six-face density

18 October 2026. **Necessary incidence consequences; full independent internal review completed.** Every entirely ideal finite-volume right-angled hyperbolic six-polytope would have at least 227 facets and 3,859 vertices. In a finite-volume right-angled polytope of dimension 11 or 12, fewer than 3/40 of the six-faces can be entirely ideal. Neither statement proves existence, excludes a dimension, or asserts an optimal published bound.

## The ideal five-dimensional input, with its proof

The frozen [carrier theorem](right_angled_five_face_carrier_20261018.md) gives a four-facet Q with positive kappa in every ideal RA5 polytope. The frozen [second-volume equality classification](ideal_ra4_second_volume_equality_20261018.md) then gives integer delta at least5. Kolpakov's published facet minimum24, with equality only for the regular24-cell, implies that this Q has at least25 facets and30 vertices.

For an actual ideal facet Q of an ideal RA polytope P, its ridge-neighbor facets and its cusp-opposite facets are distinct, including across the two lists. A cusp opposite has supporting hyperplane asymptotic to Q at exactly that ideal point, whereas a ridge neighbor is orthogonal to Q. The rectangular cusp link supplies one opposite for each vertex of Q. Thus

\[
 f_{d-1}(P)\ge1+f_{d-2}(Q)+f_0(Q),\qquad \dim P=d.
 \tag{1}
\]

Applied to the positive carrier, this proves $f_4(P^5)\ge56$, as also recorded in the separate [five-face density packet](right_angled_mixed_five_face_density_20261018.md). The logical inputs here are the two frozen internal theorems and the explicit injection (1), rather than an arithmetic inference from the density report.

Write $v=f_0(P^5)$ and $f=f_4(P^5)$. The frozen carrier proof derives

\[
 f_3(P^5)=3v+3f-6.
\]

Every four-facet is an ideal RA4 polytope, so has at least24 three-facets. Every three-face has exactly two four-cofaces. Hence $2f_3(P^5)\ge24f$, and therefore

\[
 v\ge3f+2\ge170. \tag{2}
\]

No assumption that all four-facets are 24-cells has been made; the inequality allows all larger types, including unbounded positive kappa.

## Consequences for an ideal six-polytope

Every five-facet of an ideal RA6 polytope P is itself ideal RA5 and satisfies the two lower bounds above. Applying (1) to any such facet gives

\[
 f_5(P)\ge1+56+170=227. \tag{3}
\]

An ideal vertex of P has a rectangular five-box link, with exactly ten facets. Counting vertex/facet incidences and using (2) gives

\[
 10f_0(P)=\sum_{Q\in\Omega_5(P)}f_0(Q)\ge170f_5(P).
\]

Thus $f_0(P)\ge17f_5(P)\ge3859$. These are universal necessary bounds conditional only on the explicitly listed published and internally proved statements. They are not derived from a bounded search for ideal six-polytopes.

## Mixed dimensions eleven and twelve

Use the published Nikulin–Khovanskii inequality for polytopes simple at edges, with $k=6,l=5$ and ambient dimension $n=11,12$. Both are in the permitted range $0\le l<k\le\lceil n/2\rceil$. The strict average number of facets of a six-face is

\[
 \frac1{f_6(P)}\sum_{E\in\Omega_6(P)}f_5(E)<
 \binom{n-5}{n-6}
 \frac{\binom{\lfloor n/2\rfloor}{5}+\binom{\lceil n/2\rceil}{5}}
 {\binom{\lfloor n/2\rfloor}{6}+\binom{\lceil n/2\rceil}{6}}=42.
\]

Every finite-volume RA6 polytope has at least27 facets: Dufour's published inequalities give $f_5\ge1+\nu_5\ge27$ for the noncompact case, and the published compact-dimension theorem excludes compact RA6. Every entirely ideal six-face has at least227 facets by (3). If $F=f_6(P)$ and I counts its entirely ideal six-faces, then

\[
 27F+200I\le\sum_E f_5(E)<42F,
 \qquad \boxed{I/F<3/40.}
\]

This forbids a specific boundary distribution, without forbidding all possible boundaries. It is separate from the five-face restriction $I_5/f_5<1/10$ in these dimensions; neither inequality assumes that the two sets of ideal faces are independent.

## Published premises, verification and scope

The external inputs are the [Kolpakov2012 facet theorem](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf), the absence of compact RA polytopes in dimensions at least5, and the Nikulin–Khovanskii/Dufour statements in [Alexandrov2023, Theorem1.1 and Propositions6.1,6.4,6.5](https://arxiv.org/html/2303.09533v1). Finite-volume RA faces inherit right angles and finite volume, and the ambient polytope is simple at edges, so the published average theorem applies. The original subscription Dufour proof was not independently accessible; no equality case from it is used. These source-access limits are the same as in the preceding density packet.

The [exact producer](../scripts/fresh_right_angled_ideal6_carrier_propagation_20261018.py) directly binds the frozen carrier/equality reports. It checks the rational consequences (56,170,227,3859), both NK substitutions, the strict fraction3/40 and all saved no-claim flags. It does not computationally reprove actual hyperplane distinctness, the carrier/equality theorem, or the external literature. Generation and a separate-process read-only replay both completed with exit 0/PASS in terminal chunk `f22478`.

The root agent completed a full mathematical read of the producer and report, checking the actual-facet injection in dimension6, the ideal5 Euler and facet-incidence derivation, the 227/3859 bounds and the NK range/42/general27 calculation giving3/40. This was a proof/source review, with no execution or full reproof of external theorems. No optimality, full formal-kernel or external human-review claim follows. A further peer review has not been completed and is not claimed.
