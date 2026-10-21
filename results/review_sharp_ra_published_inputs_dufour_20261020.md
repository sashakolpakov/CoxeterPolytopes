# Primary-source check of the sharpRA facet and acuteness inputs

20 October 2026. Audited manuscript:
`manuscripts/right_angled_sharp_dimension_bounds.tex`, SHA-256
`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`,
especially lines 41–45, 63–72 and 93–102. These statements were read
before the external search. No manuscript edit, mathematical worker,
Lean job or TeX job was run; no PDF assets were saved locally.

**Finding:** the stated mixed-vertex facet input is supported explicitly
by Alexandrov's original paper, and the cited Vinberg passage directly
supports the all-wall implication. No mathematical discrepancy was found.
**Access limit:** Dufour's own six-page article was not obtained in full.
This is not a direct verification of its original proof or exact internal
page locations.

## Facet minimum and the vertex scope

Read Alexandrov, *On ideal vertices of right-angled hyperbolic polyhedra*,
[original PDF](https://arxiv.org/pdf/2303.09533v1) and
[arXiv HTML](https://arxiv.org/html/2303.09533v1), 16 March 2023.

Section 1.2, PDF p.2, defines `P^n` as finite-volume **noncompact**
right-angled hyperbolic polyhedra; it does not require every vertex to be
ideal. Thus mixed finite/ideal vertices are included. Section 4, PDF p.4,
explicitly states the 27-facet minimum for every finite-volume
right-angled six-polyhedron and attributes it to Dufour Lemma 1 and
Proposition 4.

Section 6, PDF p.6, defines
`nu(P)=a_(n-1)(P)+v_infinity(P)` and minimizes this over `P^n`.
Propositions 6.4 and 6.5 restate respectively Dufour Proposition 4 and
Lemma 1 as

`nu_5 >= 26`, and `a_(n-1)(P^n) >= 1+nu_(n-1)` for `n>=3`.

Their combination gives `a_5(P^6)>=27`; PDF p.7 prints this consequence
again. Crucially, 26 is a bound on **facets plus cusps**, not just facets.
For the compact case, Alexandrov Theorem 1.1, PDF p.1, recalls the
Potyagailo–Vinberg exclusion above dimension four. Hence the noncompact
definition creates no uncovered compact-six case. No ideal-only or
arithmetic restriction appears in the needed conclusion.

## What was and was not read from Dufour

The [publisher record](https://link.springer.com/article/10.1007/s10711-009-9454-2)
confirms Guillaume Dufour, *Notes on right-angled Coxeter polyhedra in
hyperbolic spaces*, *Geometriae Dedicata* **147** (2010), 277–282,
DOI `10.1007/s10711-009-9454-2`. Its abstract concerns finite volume.
The accessible page offers subscription access; the direct publisher PDF
was unavailable. No readable author copy was located. Therefore the
numbering and formulas above are verified in Alexandrov's original text,
not independently in Dufour's original text. Dufour's exact definition and
the individual page numbers of Lemma 1/Proposition 4 remain unchecked.
The manuscript already identifies Alexandrov as the accessible statement
of this consequence, so that narrower bibliographic verification does
not contradict its current wording.

## All-wall acuteness: original Vinberg text

Read Vinberg, *Hyperbolic reflection groups*, *Russian Mathematical
Surveys* **40**:1 (1985), 31–75, in the
[original English PDF hosted by Math-Net](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&paperid=2140&what=fullteng&option_lang=eng)
([bibliographic record](https://www.mathnet.ru/eng/rm2140)).

Page 42, Section 2.1 (PDF page 12), first formulates acuteness for all
intersecting bounding hyperplanes and identifies it with nonpositive
off-diagonal Gram entries. The following remark cites Andreev [4] for
the implication from nonobtuse angles at actual codimension-two faces:
supporting hyperplanes of nonadjacent facets do not intersect, so the
all-wall condition follows. This is exactly the conversion needed by
manuscript lines 93–98; it is not merely a definition silently imposed
on right-angled polytopes. No compactness, ideal-only or arithmetic
qualification is attached there. Page 39, Section 1.1, states the finite
volume/finite ordinary-or-ideal convex-hull criterion and assumes finitely
many bounding half-spaces, consistent with the manuscript's framework.

The conversion is explicitly present in Vinberg; this audit does not
claim to have separately checked Andreev's original proof. OCR renders
some weak inequality symbols incorrectly, so the exact Dufour-recurrence
signs were also checked against the arXiv HTML math rather than copied
from PDF extraction alone.
