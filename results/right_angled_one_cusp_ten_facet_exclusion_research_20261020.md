# A one-cusp ten-facet exclusion and the exceptional-incidence budget

20 October 2026. **New paper-only local lemma, pending independent
review.** No script, census, optimizer, Lean or build ran. This does not
give a further dimension bound. In particular the improved budget below
leaves the existing cusp-eliminated dimension-eleven intercept unchanged.

## 1. The excluded three-polyhedron

**Lemma.** A finite-volume right-angled hyperbolic three-polyhedron
with exactly one ideal vertex has at least eleven polygonal facets.

Euler and polygon area give, with `F` facets, `s` cusps and `p` finite
vertices,

    p+2s=2F-4,
    sum_polygons (u+2v-5)=F+2s-12>=0.

Thus `s=1` implies `F>=10`. Suppose equality holds. Every polygon
has zero excess. A polygon without the cusp is a finite pentagon;
one with the cusp is a quadrilateral with three finite vertices.
The cusp has four incident facets. Consequently there are exactly
four such quadrilaterals and six finite pentagons.

Label the cusp quadrilaterals `Q_i`, cyclically for indices modulo
four. Their common cusp is `u`, and write their boundary cycles as

    Q_i = (u,A_i,B_i,A_(i+1)).                  (1)

Here the four edges from `u` have distinct finite endpoints `A_i`.
Consecutive `Q_i` meet along these edges, and opposite `Q_i` have
asymptotic supporting planes meeting only at `u`: this is also immediate
from the rectangular Euclidean cusp link. Each `B_i` belongs to just
`Q_i`, so the eight displayed finite vertices are distinct. The union
of the four quadrilaterals is a disc with boundary the alternating
eight-cycle in (1).

At `A_i`, let `P_i` be the third incident facet. It is a pentagon,
and its boundary contains the two edges

    B_(i-1)--A_i--B_i.                          (2)

At `B_i` the three facets are `Q_i,P_i,P_(i+1)`. In particular
`P_i != P_(i+1)`, and `P_i,P_(i+1)` share an edge through `B_i`.
We also need to rule out a repeated opposite pentagon label; it cannot
be dismissed merely by assuming the complementary dual boundary is
already simple.

Suppose `P_i=P_(i+2)`. Its intersection with `P_(i+1)` is one
convex face and contains `B_i` and `B_(i+1)`. Since these are distinct
polytope vertices, that intersection is the single edge joining them.
Together with (2) for `P_(i+1)` it gives the three-edge boundary cycle

    B_i--A_(i+1)--B_(i+1)--B_i.

A convex polygon boundary is one simple cycle. This would make
`P_(i+1)` a triangle, contradicting that it is a pentagon. Hence
the four `P_i` are distinct. Their four cyclic adjacencies are four
distinct pentagon--pentagon edges.

The four quadrilaterals have exactly eight edges adjoining pentagons,
namely the finite boundary edges in (1); all of these adjoining
pentagons are among `P_1,...,P_4`. Let `R,S` be the two remaining
pentagons. Neither adjoins a quadrilateral. Every one of its five edges
therefore adjoins another pentagon. Two convex polytope facets share
at most one edge, so `R` and `S` each adjoin all the other five
pentagons. These adjacencies give nine distinct pentagon--pentagon
edges. The four edges in the `P_i` cycle are additional, giving at
least thirteen.

But counting the six pentagon degrees gives exactly

    number of pentagon--pentagon edges=(6*5-8)/2=11.

This contradiction proves the lemma. The proof uses actual face
intersections, the finite/cusp links, planarity of the cusp star and
polygon area; it does not assume a classification or metric rigidity.

## 2. Improved universal exceptional budget

Use the notation from the reviewed polygon/three-face coupling:
`nu=F+s`, and call a three-face exceptional unless it is the canonical
`P3` or the zero-excess type `Q8`. Then

    p <=18*(nu-10)+4*(8-F)                     (3)

for every exceptional three-face. Indeed, Euler gives
`p=4F-2nu-4`, so the right side minus `p` equals

    20nu-144-8F.

At `nu=10` the exceptional case is `F=7,s=3,p=4`, giving equality.
At `nu=11`, the area bound gives `F<=10`; the lemma excludes
`F=10,s=1`, so `F<=9` and the difference is at least four.
For `nu>=12`, `F<=nu` makes the difference at least `12nu-144>=0`.
The only `nu=9` face is `P3`, excluded. The budget is zero on `Q8`.
The compact right-angled dodecahedron has `nu=F=12,p=20` and attains
(3), so this particular coefficient of `nu-10` cannot be reduced below
18 while retaining the coefficient four of `8-F`.

Thus the total exceptional finite-vertex incidence obeys

    X <=18 E+4 R,

where `E=sum_non-P3(nu-10)>=0` and `R=sum_non-P3(8-F)` remains
a signed defect, not a nonnegative count.

## 3. Exact dimension-eleven consequence and limitation

Retain the already reviewed actual-corner inequality

    X >=(11/23)*(165f-2B),

and the identities

    E=B-10a3+9a2+11520c,
    R=8a3-9a2-2B,
    a3=(9/2)a2-(165/4)f-5760c.

Substituting (3) gives

    7B-345a2+3850f+677120c>=0.                  (4)

Using the strict weighted row `a2>(25/2)f+1930c` yields

    7B>(925/2)f-11270c.                         (5)

Together with `B>(75/2)f+360c`, this still gives exactly

    2B>(117825/1379)f.

Therefore the local improvement alone does **not** raise the previously
established ordinary degree intercept. Keeping the excess instead of
discarding it gives the more informative strict row

    2758B>117825f+2254E,                        (6)

where the old coefficient from `22E+4R` was 2070. To check (6)
without a new elimination, the exact identities give

    R=-(5/4)B-(165/8)f+5760c-(3/4)E,
    360c<B-(75/2)f-E.

Substitute the second into `18E+4R >=(11/23)*(165f-2B)`;
clearing the positive denominator 46 yields (6).

This can be useful together with a genuine quantitative excess density
or an excess-dependent actual-graph opposition bound. The qualitative
fact `E>0` does not supply a uniform positive lower bound on `E/f`.
Neither a new graph row cap nor a new dimension exclusion is claimed.
