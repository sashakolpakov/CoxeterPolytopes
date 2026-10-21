# Root review: extra triangle credit at the equality carrier

19 October2026. **PASS for the new written geometric/counting argument.**
The reviewed source is
`right_angled_ideal5_carrier_equality_strictness_20261019.md`, SHA256
`df51cb742378972cce496360c3224ee9194cc16a26fd278b970c8aafafdfa09b`.
I read the complete proof and independently checked its normalization,
exceptional-triangle accounting and integer consequences. No scientific
process, enumeration, build or Lean job was run. The predecessor
outside-capacity, spinor, actual-wall and eleven-vertex carrier theorems
remain the explicit previously reviewed inputs.

## Equality and the distinguished triangles

For a36-vertex four-facet Q with an11-vertex three-facet R, the exact
defect is3. The existing bounds force q_R=0,N0_R=3,G_R=6. I checked
each alternative q_R=1,2,3 and q_R>=4: its required defect is at least
4,4,5,5 respectively. Thus every R edge has internal null product-2,
and each triangular face has all internal coefficients1.

The five nongood R vertices have at most two small polygon corners
each. At least eight triangles supply at least24 triangle corners,
so at least14 of them are at good R vertices. A triangle has only
three corners. At least five distinct R triangles therefore contain
an R-good vertex. This does not assume a classification of
eleven-vertex ideal three-polyhedra.

## The two scales and the essential loss cancellation

The internal R scale is tau=r+g; the ambient Q-in-P scale is t=n+h.
Both use actual unit walls. At a vertex of R their relation is
tau=alpha0*t. The two internal R coefficients are alpha1/alpha0
and alpha2/alpha0. Because the supporting walls are orthogonal,
restriction preserves their unit normalization.

At any R triangle corner, internal coefficient1 means that the
corresponding ambient coefficient equals alpha0. Hence if this
triangle corner is nonsmall, at least TWO ambient Q coefficients
are at least2. This justifies the B variable; an arbitrarily small
lambda=1/alpha0 is not discarded. Such a vertex lies on at most
four R triangles, so at most4B distinguished triangles fail to have
all corners small.

If a distinguished triangle is unexceptional relative to P, its
ambient coefficients are1. At its R-good corner alpha0=1, and the
other two ambient coefficients are less than2. Thus that corner is
Q-good. At most4G distinguished triangles can be unexceptional.
The remaining at least K-4B-4G are all-small exceptional triangles.
They are disjoint from the usual unexceptional-triangle credit.

Counting all polygon corners gives

    S>=6m+6f-3q+K-4B-4G,
    S<=8m+4G-4B.

The latter inequality accounts for at most4 small corners at a B
vertex rather than the usual8 at a nongood vertex. The B terms
cancel, leaving8G+3q>=6f-2m+K. The lower estimates remain valid
if T-3q or K-4B-4G is negative; they bound actual nonnegative
counts from below and are not asserted to be exact counts.

## Universal consequences

The inherited9G<=2N0 implies8G<=16N0/9<=3N0. With
V(P)=5m+N0+q and K>=5 this proves

    3V(P)>=13m+6f+5

for the equality carrier. Every positive carrier has m>=36,f>=25.
If m>=37 the older propagation gives V(P)>=211; if m=36, the
internal defect bound forces the displayed R to have exactly11
vertices and the new inequality gives3V(P)>=623. Therefore every
ideal right-angled five-polytope has at least208 vertices. Its
facet minimum62 is unchanged.

Facet/cusp injection then gives at least271 facets in dimension6.
The already reviewed universal transport68V>=1485F yields at least
5919 vertices, since68*5918=402424<402435=1485*271 and
68*5919=402492. These are stronger necessary size bounds. They
do not change the ideal dimension ceiling6 or the mixed ceiling12.

The new extra-credit geometry has written independent verification
only. The existing conditional Lean modules still check the older
206/269/5875 consequences; they have not been rerun or extended
for208/271/5919. Manuscript integration is likewise pending.
