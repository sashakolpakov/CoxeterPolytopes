# Fresh adversarial mathematical review of the sharp right-angled paper

Reviewed source: `manuscripts/right_angled_sharp_dimension_bounds.tex`.
SHA256: `3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`.

Verdict: no mathematical countercase or gap found in this fresh full read.
This conclusion is conditional on the three published inputs stated in the
paper and the standard acute-polytope foundations; it is not an assertion
that those external proofs have been independently reproved. No manuscript,
frozen evidence, catalogue, or assumption file was edited, and no proof,
scientific, or typesetting job was run. I have previously worked on the
scalar Lean bridge, but did not use its execution or earlier PASS reports
as evidence for this review.

## Actual geometry and strictness

I recalculated the projection with the paper's outward-normal convention.
For **every** other wall, the numerator is
`q·m + t(n·m) < 0`: the first summand is strictly negative, the second
nonpositive. Thus the normalized point is in the relative interior of the
actual facet, including when the maximizing wall is not unique. Its future
orientation follows from its negative product with the original future
point. There is no replacement of an actual facet by a formal hyperplane.

At a facet's actual ridge, the other defining wall is perpendicular to the
facet wall, so it retains unit norm after restriction. All intrinsic facet
heights therefore obey exactly `A -> A/sqrt(1+A^2)`. Iteration gives
`A_j <= A_0/sqrt(1+j A_0^2) < 1/sqrt(j)` for finite positive `A_0`.
Neither an attained infimum over points nor a limiting argument occurs.

Nonsimple ideal vertices do not invalidate this restriction. An intrinsic
ridge in every face used here has positive dimension and hence a finite
relative-interior point. Its transverse all-right spherical link supplies
the mutually perpendicular incident walls. The nonsimple horospherical
box link at an ideal vertex is irrelevant to that local calculation.
Finite volume and entire ideality pass to these actual faces: in the Klein
model a face is the convex hull of its finite/ideal vertices, a subset of
the ambient vertex set, and has a finite triangulation by finite/ideal
simplices. The restriction to face dimension at least three correctly
avoids ideal polygon angles. Higher-dimensional reduction to an actual
nine-face or five-face uses the same facts.

## Spherical code, tensor identity, and endpoint

Coincident spatial directions would give normal product
`cosh(asinh(t_i)-asinh(t_j)) >= 1`, contrary to all-wall acuteness. Thus
facet count equals code cardinality, even for unequal heights. The strict
height bounds yield the claimed strict products, not merely non-strict
thresholds.

I independently contracted the tensors: the two cubic--trace cross terms
are each `3x`, the trace--trace product is `24x`, and hence
`<T(u),T(v)> = x^3-3x/8`. Together with
`<A(u),A(v)> = x^2-1/6`, these give the displayed positive sum with constant
`N^2/16`. The diagonal bound is exactly `27N/16`. If `N=27`, every
off-diagonal polynomial value must vanish separately; strict `x<1/4`
leaves only `x=-1/2`. The squared norm of the sum of any four distinct
vectors is then `-2`. This correctly excludes the equality cardinality;
it does not assume uniqueness of a non-strict optimal code.

For the ideal argument I checked the unit scaling of all 24 normals
`(1, +/-e_i +/-e_j)`. One opposite spatial pair already gives `a<=1` when
both heights are at most one; future unit timelikeness forces `a>=1` and
then `w=0`. Therefore metric 24-cell rigidity, rather than kissing-code
rigidity, supplies the necessary contradiction.

## External scopes checked and remaining input boundary

- [Vinberg, p.42, Section 2.1](https://www.mathnet.ru/links/8684a695fdf254b4efd6aa83de004ad9/rm2140_eng.pdf):
  the retrieved indexed primary page explicitly records Andreev's passage
  from actual nonobtuse ridge angles to nonintersection of nonadjacent
  supporting planes and the nonpositive Gram condition. The direct PDF
  request redirected, so this check used the indexed original page text.
- [Alexandrov, Section 4](https://arxiv.org/html/2303.09533v1): the text
  explicitly gives the six-dimensional finite-volume right-angled minimum
  of 27 facets and attributes it to Dufour's Lemma 1 and Proposition 4.
  It does not restrict this sentence to entirely ideal polytopes. Dufour's
  original full text was not obtained in this bounded check; that published
  theorem remains an external input, exactly as declared in the paper.
- [Musin, pp.3--4](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf):
  the code formulation and Main theorem bound arbitrary finite subsets of
  `S^3` with distinct-pair products at most `1/2` by 24. No equality
  classification is needed or imported.
- [Kolpakov, p.1945, Theorem 4](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf):
  the original theorem gives both the minimum of 24 facets for an ideal
  right-angled four-polytope and isometry to the regular ideal 24-cell at
  equality. Its hypotheses match the descended four-facet.

I found no hidden arithmetic, compactness, vertex-simplicity, equal-height,
or code-uniqueness assumption. The universal eight-dimensional argument
uses neither Musin nor Kolpakov; the ideal four-dimensional argument uses
both, with their distinct roles preserved. The existence examples are
inherited published inputs, not constructions certified by this audit.
