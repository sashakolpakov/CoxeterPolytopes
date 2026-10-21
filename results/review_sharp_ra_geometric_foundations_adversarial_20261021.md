# Sharp right-angled dimensions: adversarial geometric audit

21 October 2026. Foundations reviewer. Audited all382 source lines of
`manuscripts/right_angled_sharp_dimension_bounds_publication.tex`, SHA256
`0c29a7869fa46e4e1a542d8ae6f4dfc7e616136c7a1aad3984d941a95cf89722`.
**No geometric gap or counterexample found at the explicit published-input
boundaries.** I reconstructed the implications below without using prior
PASS labels as evidence. I have contributed earlier reviews/editorial
clarifications of this paper; this is a fresh adversarial reconstruction,
not a claim of new outside authorship or external refereeing.

No scientific, Lean or TeX job ran. No manuscript, guide, receipt or PDF
was changed. Primary sources were read through the web, without keeping
downloaded assets.

## Actual input, not an implicit definition

Lines42–45 assume a convex, finite-sided finite-volume hyperbolic
polytope, with right angles at actual ridges and arbitrary finite/ideal
vertex mixture. The global wall sign is an additional geometric theorem,
not a formal consequence of the definition of a unit normal.

I rechecked [Vinberg1985, printed p42, Section2.1](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng).
The passage explicitly supplies the implication from nonobtuse actual
ridge angles to disjoint supporting planes for nonadjacent facets and
the nonpositive off-diagonal Gram condition, attributing the implication
to Andreev. Thus manuscript93–98 uses the stronger global conclusion
legitimately. This audit does not reprove Andreev's theorem. Vinberg's
finite-volume convex-hull criterion at printed p39 also supports the
finite/ideal vertex description used below. No ideal-only or arithmetic
condition was added to these inputs.

## Projection is onto an actual facet, with every intrinsic wall retained

For119–149 I independently compute, with t=-<q,n>,

    (q+tn)^2=-1-t²,   <q+tn,n>=0,   <q+tn,q>=-1-t²<0.

Thus normalization gives a finite future unit point, not the past sheet
or a spacelike point. For EVERY other original wall m,

    <q+tn,m>=<q,m>+t<n,m><0.

The first term is strict because q is interior; the second is
nonpositive by the cited theorem. Hence q' lies in P and is strict at
every wall except n. It is a relative-interior point of the ACTUAL facet
Q. The maximum choice is unnecessary for feasibility; its purpose is
the simultaneous bound on all new heights.

Every facet of Q is a codimension-two ridge of P and comes from the
other bounding facet at that ridge. At a finite relative-interior point
of the ridge, the two original facets meet at their actual right angle.
Consequently that neighbor's unit normal m satisfies <m,n>=0. It already
lies in the intrinsic Lorentz space n-perp, so no projection/rescaling
or missing normal is concealed in formula129–131. The finite half-space
description supplies the complete intrinsic facet list. Redundant
restricted inequalities do not add facets; a duplicate intrinsic normal
would give the same original normal because both are orthogonal to n.

The explicit finite-sided convention makes the maximum finite and
attained. After each projection the new relative-interior point is again
strict at every intrinsic wall. Thus the construction can be repeated
on the actual face, not on a formal supporting plane.

## Face inheritance survives ideal nonsimplicity

I tested99–105 against the potentially problematic mixed/ideal vertices.
A ridge of an actual face of dimension at least three has a FINITE
relative-interior point. All original facet planes containing that
point meet inside hyperbolic space. The cited nonadjacent-plane property
therefore makes each pair adjacent in P, hence orthogonal. The incident
normals are independent, and restriction to the orthogonal supporting
face preserves the right angles between its intrinsic walls. This also
justifies induction to deeper faces. No simplicity assertion at an ideal
vertex is needed: none of these local points is ideal.

In the closed Klein model, the finite-volume polytope is the convex
hull of finitely many finite/ideal vertices. The closure of an actual
face is an ordinary exposed face of that Euclidean polytope; its vertices
are a subset of the original vertices. Thus an entirely ideal polytope
cannot acquire a finite vertex in such a face. A finite triangulation of
the face uses only finite/ideal vertices and gives finite hyperbolic
volume. Mixed simplices are included. The dimensions used are at least
four (or six in the general contradiction), so the infinite length of
an ideal edge is irrelevant. No claim that an entirely ideal polygon
inherits right angles is used.

## Both descent chains and the actual code

The actual maxima satisfy A_(j+1)<=A_j/sqrt(1+A_j²).
Monotonicity or inverse squares gives152–164 exactly:

    A_j<=A0/sqrt(1+j*A0²)<1/sqrt(j).

Strictness uses finite A0 and survives tied maximizing facets. At j3
the intrinsic six-face has every t<1/sqrt3. In its Lorentz space the
spatial length of a unit normal is sqrt(1+t²), so the radial directions
in197–211 lie in R6 and their products are strictly below1/4. If two
directions coincided, their wall product would be

    sqrt((1+t_i²)(1+t_j²))-t_i*t_j>0,

contrary to the all-wall sign. The code is finite and has EXACTLY one
point per actual facet. Thus no multiplicity, ambient-dimension or
infinite-wall exception prevents application of the displayed tensor
bound. I also checked its contractions and strict endpoint: N=27 forces
all distinct products=-1/2, and four such vectors have squared sum-2.

Accordingly277–290 really produces an actual finite-volume right-angled
six-face with at most26 facets. The universal27-facet minimum is the
explicit external numerical premise that contradicts it. Any higher
dimension has a flag down to a nine-face with the same inherited
hypotheses. No monotonicity of existence or arithmetic face theorem is
being assumed.

For321–336 one projection instead gives an ACTUAL entirely ideal
four-dimensional facet with every t<1, hence a code in R4 with products
strictly below1/2. Musin's stated cardinality premise applies even with
the weaker nonstrict bound. Higher-dimensional entire ideality is reduced
by a flag to dimension five using the same inheritance just justified.

## Metric24-cell equality is the needed endpoint

I rechecked [Kolpakov2012, printed p1945, Theorem4](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf#page=5).
It gives the24-facet minimum for an ideal right-angled four-polytope
and hyperbolic isometry to C at equality. Section2 uses ideal to mean all
vertices ideal. This is the metric conclusion required in327–331, not
combinatorial recognition or uniqueness of a kissing configuration.

The model normals(1,epsilon*e_i+delta*e_j) have square-1+2=1 and
occur in opposite SPATIAL pairs. Their two heights at(a,w) are
a-r.w and a+r.w. Bounding both by1 forces a<=1, whereas future unit
timelikeness gives a>=1, with equality only when w=0. All heights then
equal1. A hyperbolic isometry preserves outward unit normals and their
heights, so Kolpakov's equality transports the projected point into
precisely this contradiction. There is no free normal rescaling that
could remove the strict t<1 obstruction.

## Limits and verdict

The manuscript proves the projection, simultaneous descent, finite-code
translation, cubic bound and metric height obstruction. It assumes the
explicit all-wall geometric theorem and standard Lorentz/Klein convex
geometry, the27-facet statement attributed to Dufour via Alexandrov,
Musin's cardinality theorem, Kolpakov's metric theorem and the published
sharp examples. I did not newly retrieve Dufour's original article or
re-referee any of these published proofs. There is no unpublished
finite-data census needed to fill one of the geometric steps above.

These inputs suffice for the advertised finite-volume8 and entirely
ideal4 conclusions. No source repair is required by this audit. This
verdict is a written mathematical assessment, not new execution evidence
or an end-to-end geometric Lean proof.
