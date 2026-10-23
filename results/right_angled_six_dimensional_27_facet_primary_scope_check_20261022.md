# Primary-source scope check: six-dimensional minimum of27 facets

22 October 2026. **Scope confirmed for the selected sharpRA external input.**
This is supporting literature documentation, not a proof execution, new
dimension bound or independent reproof of the cited literature. The selected
statement is in `manuscripts/right_angled_sharp_dimension_bounds_publication.tex`,
lines63–72. No source was changed and no scientific, Lean, TeX or replay job ran.

## Directly inspected primary text

Stepan Alexandrov, *On ideal vertices of right-angled hyperbolic polyhedra*,
**arXiv:2303.09533v1,16 March2023**, seven-page original PDF:
[versioned PDF](https://arxiv.org/pdf/2303.09533v1),
[matching versioned HTML](https://arxiv.org/html/2303.09533v1).
Page numbers below are printed PDF pages, also the PDF page numbers.

On **p4, Section4**, the exact statement is:

> Every 6-dimensional finite volume right–angled hyperbolic polyhedron has at least 27 facets.

This is the only verbatim quotation here (12 words), immediately attributed
to **Dufour, Lemma1 and Proposition4**.

On **p6, Section6**, define `nu(P)=facets(P)+ideal_vertices(P)`.
Proposition6.4 attributes `nu_5>=26` to Dufour Proposition4.
Proposition6.5 attributes `a_(n-1)(P^n)>=1+nu_(n-1)` to Dufour Lemma1.
At n6 these give `a_5(P^6)>=27`, printed again on **p7**. HTML confirms
the weak inequality signs that PDF text extraction can misread.

The family on **p2, Section1.2**, comprises finite-volume noncompact
right-angled polyhedra, allowing both mixed finite/ideal and entirely ideal
vertices. No arithmetic or special cusp assumption is imposed. The compact
six-dimensional case is separately empty by **p1, Theorem1.1** (Potyagailo–
Vinberg: no compact right-angled polyhedra in dimensions at least5).
Thus the numerical input covers every possible finite-volume vertex type.

## Dufour access boundary

Guillaume Dufour, *Notes on right-angled Coxeter polyhedra in hyperbolic
spaces*, *Geometriae Dedicata*147 (2010),277–282,
[publisher record](https://doi.org/10.1007/s10711-009-9454-2).
The publisher record was accessible, but its original PDF remained
subscription-only/inaccessible; no author copy was located. Dufour's own
internal lemma/proposition page locations and proof were therefore **not
directly rechecked**. The confirmation above is direct inspection of
Alexandrov's primary text, explicit attribution and reproduced inequalities.
It does not claim direct access to Dufour's original argument.
