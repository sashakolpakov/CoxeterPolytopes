# Fresh high-support H20 sections of the Borcherds chamber

This search finds no new Coxeter polytope. It establishes three new,
continuous restrictions on sections of the frozen finite-volume,
210-facet Borcherds polytope in H21:

1. Every interior section through either of its two displayed all-right
   finite vertices is non-Coxeter. The section normal may have **arbitrary
   real coefficients in a 21-dimensional space**.
2. Every interior section whose normal lies in the plane spanned by the
   sums of the 21 point walls and 21 line walls is non-Coxeter.
3. For **each of the 210 original walls**, every interior section with
   normal in the plane spanned by that wall and the chamber center is
   non-Coxeter. All real parameters and all possible angle labels are covered.

Each rejection is an actual forbidden dihedral ridge, not a raw projected
Gram entry between possibly redundant inequalities. Other section normals,
other ambient polytopes and general H20 constructions remain open.
Rational sections were not rejected for arithmetic reasons: rational
reflection data alone do not imply integrality of the cyclic products.

The starting [wall packet](../data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json)
and [complete corner walk](fresh_h21_source_corner_walk_20261013.md) are
frozen dependencies. The two new checkers use only Python's standard
library, integer arithmetic and `Fraction`. A numerical linear program
was used only to propose two endpoint/dual certificates; their final replay
does not import or run an optimizer.

## Why these sections were new search directions

The earlier [two-wall theorem](fresh_h20_global_20261013.md) covers normals
in a span of two original walls; the [three-wall results](fresh_h20_triple_sections_20261013.md)
and [elliptic-span theorem](fresh_h20_bivariate_sections_20261013.md) cover
other specified sparse normal families. The present constructions use a
whole orthogonal vertex span, the two 21-wall sums, or the full chamber
center. They are not another coefficient enlargement of a prior sparse
pencil. Their claims concern the declared spaces, without asserting a
minimum original-wall support for every specialization.

An ordinary symmetry-fixed codimension-one section offers no untested
alternative: the [complete 28-class cyclic census](borcherds_char2_all_cyclic_element_completion.md)
has no nonidentity fixed rank 21. That does not exclude hyperplanes merely
invariant under a subgroup, or generic sections.

Every genuine interior hyperplane section used here has finite volume.
In the Klein model the original chamber is a polytope with finitely many
finite or ideal vertices. A hyperplane intersection is again a polytope;
strict convexity of the ideal sphere forces its ideal vertices to be
inherited ones. Its finite/ideal simplex triangulation has finite volume.
Thus the obstruction is the Coxeter-angle condition, not a failed volume
test or an arithmetic-dimension shortcut.

## All sections through an all-right vertex

The 42 short roots split into the point and line sets of `PG(2,4)`. Within
either set the 21 norm-two roots are mutually orthogonal. Write their sums
as `a,b`, and use the strict center `h=(a+b)/3`. Exact source identities are

```text
a²=b²=42,  (a,b)=-105,  h²=-14.
```

The two points `v=h+a/2` and `h+b/2` have square `-49/2`. Each is
orthogonal to exactly its 21 short walls, and strict on every other wall.
These 21 positive orthogonal directions and the timelike vertex also
certify the ambient signature `(21,1)`.

Fix one vertex and write the section normal in its perpendicular space as
`e=sum c_i r_i`. A section meeting the chamber interior requires both
positive and negative nonzero coefficients. Let their counts be `p,q`;
zero coefficients are unrestricted.

If `p,q>=2`, choose one wall from each sign class. Their projected pairing
is `-4 c_i c_j/e²>0`. It is an **actual** ridge: after setting these two
local wall coordinates to zero, both signs still occur in the remaining
equation. Choose every remaining local wall coordinate negative while
satisfying that equation, then perturb the old finite vertex slightly in
this direction. The point stays future timelike and strict on every old
wall except the chosen two. The angle is obtuse.

If one sign class has one member and the other at least three, choose three
walls from the larger class. Each pair is an actual ridge by the same
local construction. Their projected Gram is positive definite, since the
section normal has a component outside their positive span. Every mutual
pairing is strictly negative. If all three angles were Coxeter, every
normalized off-diagonal entry would be at most `-1/2`. Its value on the
all-ones vector would then be nonpositive, contradicting positive
definiteness. This argument covers every label, with no degree bound.

It remains to check supports two and three. The exact incidence census
checks all 7,980 ordered triples in each bipartition: for distinct `p,q,r`
there are three or four opposite-side short walls `l` incident to `q`
and orthogonal to `p,r`.

* For `e=A p-B q`, with `A,B>0`, choose such an `l`. The ambient
  `{p,q,l}` positive-definite face gives an actual section ridge: the
  projected `p,q` are positive multiples of one another. Its squared
  cosine is `A²/(4A²+3B²)`, strictly between zero and `1/4`.
* For `e=A p-B q-C r`, with `A,B,C>0`, choose `l` incident only to `q`.
  There are strict timelike ambient face points `z1` on `{r,p,l}` and
  `z2` on `{r,q,l}`. Their e-pairings have opposite signs. The positive
  combination `(e,z1) z2-(e,z2) z1` lies in the section and vanishes only
  on `r,l`. Their projected pairing is `2BC/e²>0`, an actual obtuse ridge.

The required ambient face points follow from strict-center projection
onto an acute positive-definite Gram with nonnegative inverse. The checker
also saves explicit representative points and all 210 wall signs. Reversing
the overall normal handles the opposite choice of singleton sign. This
exhausts every real normal through either vertex, including vanishing
coefficients and all higher-support cases.

## The point-sum/line-sum plane

Every interior normal in `span(a,b)` is, up to nonzero scale,
`e(t)=a-t*b`, `t>0`. Same-sign coefficients cannot vanish at a strict
interior point. A strict section center is

```text
c(t)=(5+2t)a+(2+5t)b.
```

For two point walls `p1,p2`, put

```text
A=2a-2p1-2p2,  B=2b+5p1+5p2,
x(t)=(95+17t)A+(38+95t)B.
```

For `0<t<=1`, `x(t)` is a future timelike section point, zero on precisely
`p1,p2`, and strict on all other 208 walls. Its square is
`-33516(38+190t+17t²)`. The exact cosine of this actual ridge is

```text
(2+5t)² / (80+400t+59t²).
```

It is strictly between zero and `1/2`; subtracting twice its numerator
from the denominator gives `72+360t+9t²>0`. Thus its acute angle is strictly
between `pi/3` and `pi/2`. The symmetric construction using two line walls,
with parameter `1/t`, covers `t>=1`. All polynomial identities and all
210 parameter-dependent wall signs are checked exactly.

## Every center/facet pencil

Let `r_i²=n` and `(h,r_i)=-d`. The exact symmetry generators have two
wall orbits, `(n,d)=(2,1)` and `(4,8)`. Every possible spacelike normal
in `span(h,r_i)` can be written, up to scale, as `e=r_i-t*h`.

The complete interval in which its hyperplane meets the strict chamber
interior is

| Wall type | Number | Interior interval | Actual ridge walls, zero based |
| --- | ---: | --- | --- |
| short, representative 0 | 42 | `0<t<1/3` | 2,4 |
| long, representative 19 | 168 | `0<t<1` | 3,5 |

These are certified intervals, not LP estimates. For each representative
the archive supplies `T*h-r_i=sum lambda_j r_j`, all `lambda_j>0`, and a
feasible future point `u` with `(h,u)=-1` and `(r_i,u)=-T`. The dual identity
excludes every `t>=T`; positivity of the old inequalities excludes `t<=0`.
Mixing `u` with the relative-interior point of facet i proves existence
for every `0<t<T`.

The two listed ridge walls are short, mutually orthogonal and orthogonal
to `r_i`; they also vanish at `u`. Let `z` be the strict timelike ambient
triple-face point, normalized by `(h,z)=-1`. Then

```text
x(t)=(T-t)z+t*u
```

is an actual section-ridge point throughout `0<t<T`: precisely the two
ridge inequalities vanish, while the former wall i and all other walls
are strict. The projected Gram gives

```text
cos(theta)=t²/(2n+4dt-29t²).
```

The concave polynomial `2n+4dt-31t²` is positive at both interval endpoints
(`4,17/9` in the short case; `8,9` in the long case), and hence positive
throughout. Again `0<cos(theta)<1/2`, so no Coxeter label is possible.
All six recorded ambient symmetry generators are checked as exact form
isometries preserving the entire wall set and h; orbit traversal covers
all 210 walls. No unproved symmetry transfer is used.

## Artifacts, replay and next restriction

* [Orthogonal-vertex and bipartition checker](../scripts/fresh_h20_borcherds_sections_20261013.py),
  [exact certificate](fresh_h20_borcherds_sections_20261013.json).
* [Center/facet checker](../scripts/fresh_h20_borcherds_sections_20261013_facet_verify.py),
  [exact certificate](fresh_h20_borcherds_sections_20261013_facet_verify.json).
* [Small LP proposal script](../scripts/fresh_h20_borcherds_sections_20261013_facet_scout.py),
  [frozen rational recipes](fresh_h20_borcherds_sections_20261013_facet_recipes.json).

```sh
python3 scripts/fresh_h20_borcherds_sections_20261013.py --verify
python3 scripts/fresh_h20_borcherds_sections_20261013_facet_verify.py --verify
```

Both saved replays pass. A separate agent completed a
[source and theorem review](fresh_h20_borcherds_sections_20261013_review.md)
and reran both checkers; that review reuses the exact checkers and is not
a second arithmetic implementation. These are exact restrictions on the stated
continuous families, rather than a finite failed angle search. A section
continuation must avoid the two all-right vertices and use a normal outside
the bipartition-sum plane and every center/facet plane. No claim is made
that such remaining sections are impossible or that any new dimension
has been constructed.
