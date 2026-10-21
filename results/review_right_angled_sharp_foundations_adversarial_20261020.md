# Independent adversarial source audit: sharp right-angled bounds

Date: 20 October 2026. Reviewer: separate source-only mathematical agent.

Audited source: `manuscripts/right_angled_sharp_dimension_bounds.tex`,
SHA-256 `3c3ddd8540593004bcd5fb57832afc77ffc639268e66cb0ab29302b286654310`,
15,449 bytes, 378 lines. The source was read and the arguments below
reconstructed before opening any earlier review; no earlier review was
opened. No scientific worker, Lean job, TeX job, or external source lookup
was run. This is mathematical source evidence, not an executable or
end-to-end formal verification claim.

**Finding:** no internal counterexample or proof gap was found at the
source's declared geometric and published-theorem boundaries. The
dimension exclusions follow from those inputs. In particular, the
projection is into an actual facet, its complete intrinsic facet list
has the stated heights, the strict spherical endpoint is justified, and
the ideal equality argument uses metric rigidity rather than kissing-code
uniqueness.

## 1. Foundations and the exact imported scope

Lines 93–102 explicitly import the acute-polytope all-wall theorem:
actual nonobtuse dihedral angles imply nonpositive Lorentz products of
all distinct outward unit facet normals, including nonneighbors. This
is essential. Local right angles alone must not be substituted directly
for all-wall orthogonality; disjoint walls generally have negative
products. The proof consistently uses the nonpositive version until an
actual ridge supplies product zero. This audit accepts the cited
Vinberg foundation, rather than claiming to have reproduced or checked
its published proof.

The numerical input in lines 63–72 has the exact required scope:
**every finite-volume right-angled six-polytope**, with arbitrary mixtures
of finite and ideal vertices, has at least 27 facets. The general proof
does not assert that its selected six-face is entirely ideal and would
not follow from an ideal-only minimum. Conditional on the stated
Dufour/Alexandrov proposition, there is no scope mismatch. The other
external inputs are Musin's cardinality bound for an arbitrary finite
subset of S^3 with inner products at most 1/2, and Kolpakov's metric
rigidity for the equality case of the entirely ideal four-dimensional
24-facet minimum. Neither a uniqueness theorem for spherical kissing
codes nor an arithmeticity assumption is needed. Sharpness also imports
the cited eight-dimensional example and the ideal 24-cell.

Here is an independent justification of the inherited-face assertions
used after these foundations. Near an ideal vertex, a sufficiently
small horosphere cuts a finite-volume polytope in a bounded Euclidean
cusp section; a face cuts it in a face of that bounded section. Thus
the face has finite cusp sections. The parts away from the finitely
many ideal vertices are compact, so faces of intrinsic dimension at
least two have finite intrinsic volume. Vertices of a face are vertices
of the original compactified convex polytope: hence an entirely ideal
polytope has entirely ideal faces. No new finite vertex is created by
passing to a supporting face.

For right-angle inheritance in an intrinsic face of dimension at least
three, an intrinsic ridge has positive dimension and therefore contains
a finite point. Locally one can take its finite spherical normal link.
For one facet restriction, the link relevant to an intrinsic ridge is
a spherical polygon with every actual interior angle pi/2. Its area is
`(2-k/2) pi > 0`, so it has three sides. The three incident original
facets are pairwise perpendicular, and restriction to one leaves a
right intrinsic dihedral angle. Repeating facet restriction establishes
the required inheritance. Ideal vertices need not be simple and are
not being assumed simple in this argument. A two-dimensional face
would require a different discussion at ideal vertices; no chain in
the manuscript reaches dimension two.

## 2. Actual-facet projection and completeness of intrinsic normals

For the data in lines 115–147, put `v=q+t n`, with `q.n=-t`,
`q^2=-1`, and `n^2=1`. Direct calculation gives

`v.n=0`, `v^2=-1-t^2`, and `v.q=-1-t^2<0`.

Consequently `q'=v/sqrt(1+t^2)` is finite and future unit timelike in
the supporting hyperplane. For **each** other actual wall `m`,

`q'.m=(q.m+t(n.m))/sqrt(1+t^2)<0`,

because `q.m<0`, `t>0`, and the imported all-wall theorem gives
`n.m<=0`. This checks every original inequality, including nonadjacent
walls. Hence the foot is in the relative interior of the actual facet;
the proof has not silently replaced that facet by its supporting plane.

Every intrinsic facet of that actual facet is a ridge of the original
polytope. A relative-interior finite point of a convex-polytope ridge
has exactly two incident facets: the chosen facet and its unique ridge
neighbor `m_i`. Their angle is pi/2, so `n.m_i=0`. Thus `m_i` is already
a unit vector in the intrinsic Lorentz subspace `n^perp`, with the same
outward orientation, and its inequality is the intrinsic facet
inequality. Restrictions of other original walls can be redundant or
meet lower-dimensional faces; they supply no missing intrinsic facets.
This proves the completeness needed for the height maximum to control
**all** facets of the next face, even when the original polytope has
nonsimple ideal vertices.

It follows exactly, without an estimate on nonneighbor heights, that

`t'_i=t(m_i)/sqrt(1+t^2)<=t/sqrt(1+t^2)`.

If the old maximum is at most `A`, monotonicity of
`f(x)=x/sqrt(1+x^2)` gives new maximum at most `f(A)`. Repeated use of
`f(x)^(-2)=x^(-2)+1` gives

`A_j<=A_0/sqrt(1+j A_0^2)<1/sqrt(j)`.

There are finitely many facets and the starting point is finite, so
`A_0` is finite. This is the strictness source; unique maximizing
facets, generic position, and a limiting point are unnecessary.

## 3. Spherical translation and the strict 27-point endpoint

In a frame whose time axis is the selected finite point, an intrinsic
unit outward normal is `(t_i,sqrt(1+t_i^2) u_i)`. Nonpositive Lorentz
product gives

`u_i.u_j <= [t_i/sqrt(1+t_i^2)] [t_j/sqrt(1+t_j^2)]`.

Equal directions cannot merge two actual facets: if `u_i=u_j`, writing
`t_i=sinh(a_i)` gives their Lorentz product
`cosh(a_i-a_j)>=1`, contradicting nonpositivity. The code cardinality
therefore equals the facet number. After three descents in dimension
nine, every height is strictly below `1/sqrt(3)` and each bracket is
strictly below `1/2`. Hence all distinct code products are strictly
below `1/4` in six Euclidean variables.

I independently expanded the tensor calculation in lines 216–270.
For dimension six, `A(u)=u tensor u-I/6` has kernel `x^2-1/6`.
For `S(u)` as in the source, the contractions are
`u^3.S(v)=3x` and `S(u).S(v)=24x`; therefore
`T(u)=u^3-S(u)/8` has kernel `x^3-3x/8`.
The displayed decomposition of
`F(x)=(x+1/2)^2(x-1/4)=x^3+3x^2/4-1/16`
then has constant term `1/16` and positive tensor-square coefficients.
Summing yields `N^2/16 <= sum F(u_i.u_j) <=27N/16`.

At the potentially problematic endpoint `N=27`, both inequalities are
equalities, so every nonpositive off-diagonal summand must vanish.
Strictly below `1/4`, the only zero of `F` is `-1/2`. Four of the
27 vectors would then have sum of squared norm
`4+12(-1/2)=-2`, which is impossible. Thus at most 26 facets remain,
contradicting precisely the mixed-vertex six-dimensional input above.
Replacing strict `<1/4` by `<=1/4` would invalidate this endpoint
argument; the manuscript does not make that replacement.

## 4. Ideal endpoint and reduction of higher dimensions

One descent from a hypothetical entirely ideal five-polytope produces
an actual entirely ideal right-angled finite-volume four-facet with
all heights `<1`. Its directions have products `<1/2`; Musin still
applies because this is a subset of the permitted `<=1/2` domain.
The facet has at most 24 facets, whereas Kolpakov gives at least 24
and identifies equality up to hyperbolic isometry with the regular
ideal right-angled 24-cell.

The displayed model normals `(1, +/-e_i +/-e_j)` have Lorentz norm
one and come in opposite spatial pairs. For a future unit point
`(a,w)`, the two heights for a pair are `a-r.w` and `a+r.w`.
If every height were at most one, their sum would force `a<=1`;
the unit equation forces `a>=1`, so `a=1,w=0`. All heights are then
exactly one, contradicting the actual projected heights `<1`.
Hyperbolic isometry preserves the unit normalization and these scalar
heights, so the metric rigidity input applies to this calculation.
The argument does not infer spherical-code uniqueness from Musin.

Finally, repeatedly taking actual facets produces a nine-face in any
higher-dimensional general counterexample, or a five-face in any
higher-dimensional entirely ideal counterexample. Every intermediate
face has dimension at least the target dimension, so the preceding
inheritance arguments apply throughout. No ideal two-face, nonactual
supporting-plane section, or finite-volume inference from a few vertices
is involved.

## Status and limits

The reconstructed arguments support both stated upper bounds at the
manuscript's declared boundary. No canonical source edit is requested.
The published input statements, the Vinberg all-wall foundation, and
the cited examples were treated as the explicit external premises
already isolated by the manuscript; this audit is not a fresh external
bibliographic verification or external refereeing. It makes no new
runtime, PDF, or Lean claim.
