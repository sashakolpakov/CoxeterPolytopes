# Full mathematical audit of the sharp right-angled dimension paper

Date: 20 October 2026. Reviewer: resource_cleanup.

**Verdict: PASS.** I found no mathematical gap in either exclusion, their
stated scopes, or sharpness. No source correction is requested.

Audited source: `manuscripts/right_angled_sharp_dimension_bounds.tex`, all
378 lines, SHA256
`3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`.
The source and its frozen verification packets were not changed. I
re-derived the geometric and tensor arguments before consulting the primary
literature below; I did not read earlier project reviews or use the exact
workers. This is a fresh full-paper audit, not independent discovery: I
participated in the earlier development of this argument and know its history.
No scientific, Lean, or TeX job was run.

## Foundations and external inputs

The hypotheses throughout are a full-dimensional convex polytope, finitely
many actual facets, finite hyperbolic volume, and right angles at all actual
ridges. The finite/ideal vertex mixture is unrestricted in the first theorem.
The following are explicit external inputs, not conclusions of the audit.

* Vinberg, *Hyperbolic reflection groups*, printed p.42, Section 2.1,
  explicitly records Andreev's implication from nonobtuse actual ridge angles
  to disjoint nonadjacent supporting hyperplanes and nonpositive
  off-diagonal outward-normal Gram entries. I read the relevant passage in
  the [primary English text](https://www.mathnet.ru/php/getFT.phtml?jrnid=rm&option_lang=eng&paperid=2140&what=fullteng).
  Thus lines 93–110 apply to every pair of actual walls; they do not assume
  that all pairs share a ridge.
* The universal six-dimensional minimum is verified in
  [Alexandrov's primary paper](https://arxiv.org/html/2303.09533v1), Section 4,
  and again by Propositions 6.4–6.5: the quoted Dufour inequalities are
  `nu_5 >= 26` and `f_(n-1) >= 1 + nu_(n-1)`, giving at least 27 facets in
  dimension six. The statement includes all vertex types. Dufour's original
  article is [subscription content](https://link.springer.com/article/10.1007/s10711-009-9454-2);
  I did not reread its original proof in this audit. This is the same declared
  external numerical input as Proposition 3, not an additional assumption.
* [Musin's published Main theorem](https://annals.math.princeton.edu/wp-content/uploads/annals-v168-n1-p01.pdf),
  printed p.4, gives kissing number 24 in dimension four. The paper's
  reformulation as a finite spherical code in S^3 with distinct-pair
  products at most 1/2 is exact. No uniqueness result is imported.
* [Kolpakov's published Theorem 4](https://msp.org/agt/2012/12-4/agt-v12-n4-p03-s.pdf),
  printed p.1945, gives at least 24 facets for an ideal right-angled
  four-polytope and says equality is isometry to the ideal hyperbolic
  24-cell. This is metric rigidity, not only a combinatorial identification.
  The theorem number and page in the manuscript are correct.
* Sharpness in dimension eight is furnished by the actual finite-volume
  right-angled P_8 in
  [Potyagailo–Vinberg, Section 3 and Proposition 3.1](https://ems.press/content/serial-article-files/43032),
  printed pp.66–67. That construction is distinct from the later universal
  27-facet minimum; the manuscript attributes the latter correctly.

Here is also an explicit check of the compressed face-inheritance statement
at lines 99–102. A finite-volume finite-sided hyperbolic polytope is the
convex hull of its finite and ideal vertices; each actual face is the convex
hull of the vertices on its supporting plane. Hence faces of the dimensions
used here (at least three) have finite intrinsic volume, and an entirely ideal polytope has
entirely ideal faces. Each ridge is contained in exactly two actual facets.
For a facet Q with normal n, its intrinsic facet normal comes from that
unique other facet m and satisfies n.m=0. If two such intrinsic walls meet
along a ridge of Q, their supporting hyperplanes meet at a finite point of
that ridge. The all-wall theorem makes the corresponding ambient facets
adjacent, hence their normals orthogonal. Their restrictions therefore meet
at a right angle in Q. Iterating proves precisely the face inheritance used
here. No general-angle assertion that every metric face is Coxeter is used.
Nonsimplicity at ideal vertices does not affect this argument.

## Projection, strict descent, and the facet code

I checked lines 106–207 with signature (-,+,...,+) and outward *unit*
spacelike normals. If t=-q.n is the maximum positive height, the numerator
q+tn has square -1-t^2, is perpendicular to n, and has negative product
-1-t^2 with the future point q. Thus its normalization is future unit
timelike. For every other actual wall m,

    (q+tn).m = q.m + t(n.m) < 0.

The first term is strictly negative and the second nonpositive. This proves
actual relative-interiority of the projected point, including all walls
that are not neighbors of n. For a ridge neighbor m_i, restriction to
n-perp preserves its unit norm and gives height
`t(m_i)/sqrt(1+t^2) <= t/sqrt(1+t^2)`.
There is no hidden rescaling or need to choose a finite vertex.

The increasing function `a/sqrt(1+a^2)` gives by induction the bound
`A_0/sqrt(1+j*A_0^2)`. A_0 is a positive finite maximum because there are
finitely many walls and the starting point is finite interior. Consequently
the final inequality is strict against `1/sqrt(j)`, even when a maximum is
tied. Three descents from dimension nine end in an actual six-face; one
descent from dimension five ends in an actual four-facet. Neither chain
uses a two-face or a projection from a cusp.

At the resulting point, a unit normal is
`m_i=(t_i,sqrt(1+t_i^2)*u_i)`. The all-wall inequality gives exactly the
displayed bound `u_i.u_j <= t_i*t_j/sqrt((1+t_i^2)(1+t_j^2))`.
Equal directions would give a positive Lorentz product, since
`(1+t_i^2)(1+t_j^2)-t_i^2*t_j^2 = 1+t_i^2+t_j^2 > 0`.
Thus facets inject into the spherical code. The strict bounds on each
`t_i/sqrt(1+t_i^2)` give distinct-pair products below 1/4 after three
descents, and below 1/2 after one descent. Cardinalities are actual facet
counts, not counts of a selected subset of walls.

## The six-variable cubic and general dimension bound

I independently contracted every tensor in lines 224–255. For x=u.v,

    <A(u),A(v)> = x^2 - 1/6,
    <u^3,S(v)> = <S(u),v^3> = 3x,
    <S(u),S(v)> = (18+6)x = 24x,
    <T(u),T(v)> = x^3 - 3x/8.

Substituting into the displayed kernel identity gives exactly
`F(x)=x^3+3x^2/4-1/16`, and summing gives the four nonnegative square
terms with constant `N^2/16` shown in the tensor-square identity. In the strict pair domain,
every off-diagonal value of F is nonpositive; the diagonal is 27/16.
Hence N<=27. At N=27, equality of the two sums forces every off-diagonal
value individually to vanish. The only allowed root is -1/2, and already
four such unit vectors have sum of squared norm -2. Thus N<=26.
The strictness and equality argument are essential and valid; no spherical
linear-programming theorem, numerical optimization, or code uniqueness is
needed.

This contradicts the external 27-facet minimum for the inherited
six-face. Every higher-dimensional polytope has an actual nine-face with
the same hypotheses, so the reduction proves the complete stated
finite-volume bound, including compact and mixed polytopes.

## Entirely ideal bound and metric equality

The four-facet obtained from an entirely ideal five-polytope inherits all
three required properties. Musin gives at most 24 facets and Kolpakov
gives at least 24, with the required isometry at equality. In the standard
24-cell model the 24 normals `(1,+/-e_i+/-e_j)` are unit spacelike;
their spatial signs occur in opposite pairs. For a future unit point
`(a,w)`, the paired heights are `a-r.w` and `a+r.w`. Both being at most
one forces a<=1, while unit timelikeness forces a>=1. Hence a=1 and w=0,
and all heights equal one. The strict heights produced by projection are
therefore impossible. Isometry transports outward unit normals without
changing these heights. No assumption about uniqueness of a 24-point
kissing configuration enters this step.

Entirely ideal higher-dimensional polytopes have an actual ideal five-face,
so the quantifier in Theorem 2 is justified. The regular ideal right-angled
24-cell supplies equality in dimension four. The diagram and caption track
the exact descent and do not introduce a further geometric premise.

## Audit boundary

This is a complete mathematical read of the displayed proofs, relying on
the named published inputs and standard face geometry identified above.
It is not a new audit of the full proofs of Musin, Kolpakov, Dufour, or
Vinberg, a novelty assessment, a typesetting review, or a fresh formal or
computational verification. Within that boundary both theorems and their
sharpness pass without a requested change.
