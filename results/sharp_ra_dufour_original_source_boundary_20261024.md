# Sharp RA: Dufour source boundary

24 October 2026. Author: `audit_sharp_ra`. This is the bounded follow-up
to Section 7 of `results/audit_sharp_ra_fresh_20261024.md`. No proof job,
test or compilation was run, and no canonical manuscript was edited.

**Status:** Dufour's original full article was not obtained. Its exact
original hypotheses, internal numbering and printed page locations are
therefore still **not independently authenticated**. The stated
27-facet premise is explicitly supported by Alexandrov's accessible
paper; that attribution is secondary evidence about Dufour, not a
substitute for reading Dufour himself.

## Precisely authenticated accessible statement

In [Alexandrov, arXiv:2303.09533v1](https://arxiv.org/html/2303.09533v1),
Section 4 explicitly gives at least 27 facets for every finite-volume
right-angled hyperbolic six-polyhedron, citing Dufour's Lemma 1 and
Proposition 4. Section 6 defines the noncompact finite-volume family,
which includes mixed finite/ideal vertices, and
`nu(P)=number_of_facets(P)+number_of_ideal_vertices(P)`.
Its Propositions 6.4–6.5, attributed to those Dufour results, give
`nu_5>=26` and `facets(P^n)>=1+nu_(n-1)`. Thus `facets(P^6)>=27`, also
printed there. The number 26 is not asserted to count facets alone.
The PDF locations are Section 4, p. 4, and Section 6, pp. 6–7.
Alexandrov also recalls the compact right-angled exclusion above
dimension four, so there is no separate compact-six case to cover.

These scopes match Proposition `prop:27` at lines 63–72 of
`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`.
The manuscript discloses Alexandrov as the accessible statement of the
Dufour consequence. Its application is to a finite-volume RA6 polytope
of unrestricted vertex type. Neither entirely ideal vertices nor an
equality classification for 27 facets is required.

## Original-source retrieval assessment

The [publisher's original record](https://link.springer.com/article/10.1007/s10711-009-9454-2)
confirms G. Dufour, *Notes on right-angled Coxeter polyhedra in hyperbolic
spaces*, Geometriae Dedicata 147 (2010), 277–282, DOI
`10.1007/s10711-009-9454-2`. It exposes the abstract and subscription
access, not the needed propositions. The direct publisher PDF route
was unavailable. No subscription action, payment or contact with
another person was attempted.

The repository's available reference files contained no Dufour original.
An exact DOI query to the HAL public metadata API returned zero records;
the author-name query located the distinct 2012 cubulation thesis but
no deposit of this article. OpenAlex's record W1986003067 listed no
open-access copy or repository full text. Those metadata observations
describe retrieval only; they do not prove nonexistence of an author copy.
The attempted historical Orsay author hosts could not be resolved.
The inaccessible full article is not treated as read from search snippets.

This bounded search ends with the original-source question open. Release
documentation should retain one exact external premise — **RA6,
finite volume, all vertex types, at least 27 facets** — and distinguish:

- mathematical use: exact scope explicitly stated in Alexandrov;
- bibliographic attribution: Dufour Lemma 1 plus Proposition 4, as reported
  by Alexandrov;
- original-source authentication: pending access to the complete article;
- proof/formalization: this remains an external geometric input to the
  conditional internal proof.

This is an access/provenance limitation, not a newly found mathematical
counterexample or a silently narrowed ideal-only hypothesis. It does
not change the current bound or supply a new runtime verification.

The existing Lean interface already isolates this precise statement as
`PublishedDufour.minimumSixFacets`; no additional arithmetic, entirely
ideal or equality hypothesis belongs in that interface. The release
proof map and assumption register should link this access assessment
without upgrading its authentication status.

Bound local sources: publication TeX SHA256
`0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`;
fresh sharp-RA audit SHA256
`ed69a9ff3dcbbbfe10595e3fb4d0d496783241ba83896e89a26775b48f505584`;
earlier published-input access audit SHA256
`a76fb89b17ad71721916aba75b71a6c641fe319c68b5450ce2422d65107058e7`.
