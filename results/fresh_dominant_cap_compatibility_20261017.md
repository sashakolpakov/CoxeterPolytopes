# A geometric corner constraint on the dominant raw23 prisms

At a finite vertex of type

\[
Y=A_1+A_3+H_3+H_4+I_2(5),
\]

at most one of the selected corners of two specified raw23 prism classes
can occur. The resulting necessary whole-face incidence inequality is

\[
N_{26448}+2N_{26474}\le N_Y.
\tag{1}
\]

It is violated by the current exact partial-model witness. Its left side
is310256599/10000000000, whereas the saved allowance for vertex3137 is
200076613/10000000000. The gap is **55089993/5000000000**. A monotone
repair of that single allowance reduces the fixed witness's ratio from
12.11501285634343 to **11.986190572394824**. This crosses12 for this witness
only: no optimized feasible rule, complete geometric coverage, or new
dimension theorem is asserted.

The [source](../scripts/fresh_dominant_cap_compatibility_20261017.py) and
[certificate](fresh_dominant_cap_compatibility_20261017.json) reconstruct
the common finite vertex, all six literal corner orientations, every
provenance map, and all15 pair exclusions. They bind the
[preceding witness](fresh_sqrt5_scaled_core_cap_edge_repair_20261017.md).
No optimizer or historical large scan is run.

## Complete wall contacts fix the missing normal

Let v be the future unit timelike vector of a finite vertex of an acute
hyperbolic Coxeter polytope. Normalize its13 incident outward normals
to norm2, with positive definite Gram matrix V. If an external actual
wall has complete contact covector p with those normals, its projection
s(p) onto their span satisfies

\[
s(p)^2=p^{\mathsf T}V^{-1}p.
\]

The wall is not incident to v, so its outward normal r has r·v<0.
Consequently

\[
r=s(p)+\sqrt{p^{\mathsf T}V^{-1}p-2}\,v,
\qquad p^{\mathsf T}V^{-1}p-2>0.
\tag{2}
\]

Both the sign and magnitude are fixed. Two external actual normals thus
have the forced product

\[
r\cdot r'=p^{\mathsf T}V^{-1}p'
-\sqrt{(p^{\mathsf T}V^{-1}p-2)
        (p'^{\mathsf T}V^{-1}p'-2)}.
\tag{3}
\]

All covectors must be expressed in the same complete labelled vertex
frame. Equal covectors describe the same outward normal and are not
treated as an incompatible pair. Formula(3) adds geometric information
that is absent from abstract disjoint port packing.

## The two actual prism classes and their orientations

The frozen actual-provenance map identifies Q26448 with provenance490
and Q26474 with provenance497. Both have canonical prism core19 and
literal raw source23, a five-dimensional prism with ten quadrilateral
subfaces and ten finite vertices. Their full stars have15 normals after
the specified orthogonal spectator forests are included.

The producer reconstructs all ten corners of each prism. Exactly one
corner of Q26448 has forest Y: the omitted raw roots are(0,4). Exactly
two corners of Q26474 have forest Y: omitted roots(0,3) and(2,3).
These are the coefficients1 and2 in(1); they count actual vertex
occurrences, not automorphism presentations.

In the canonical Y frame the components occupy

\[
A_1:\{0\},\quad A_3:\{1,2,3\},\quad
H_3:\{4,5,6\},\quad H_4:\{7,8,9,10\},\quad
I_2(5):\{11,12\}.
\]

The component types are distinct. Their complete weighted automorphism
group consists of A3 reversal and I2(5) reversal, independently: four
permutations. The source explicitly constructs these permutations,
checks the weighted graph, and verifies that the complete graph
isomorphism enumeration agrees. Applied to the three corner records,
they produce four distinct orientations of the Q26448 corner and one
orientation of each Q26474 corner: six literal full stars in total.

The full external contact pairs are retained, not just the occupied
port pairs. The source also verifies their positive timelike squared
heights and the prescribed product between the two omitted normals
using the common future sign in(2).

## Seven pairs share an actual quadrilateral

All four Q26448 orientations occupy the port pair{1,3}, the endpoints
of the A3 component. Both Q26474 orientations occupy the four pairs
{0,i}, i∈{7,8,9,10}. Hence six pairs within the first class and the one
pair within the second class share an actual quadrilateral port pair.

At a simple finite vertex a port pair determines one actual two-face.
A quadrilateral has a unique complete set of its supporting and side
walls. The [canonical full-star aggregation theorem](fresh_quad_product_aggregation_20261017.md)
therefore assigns it to just one quad aggregate. The two literal
contact pairs in each comparison differ; by(2) their full wall stars
are distinct. They cannot both realize that same quadrilateral.

This excludes seven of the15 pairs. It does not incorrectly exclude
two normal vectors merely because they coincide.

## The other eight pairs force a forbidden angle

Every Q26448 orientation contains an external actual wall touching
only an A3 endpoint, by the dotted squared normalized product

\[
\tau=(7+\sqrt5)/2.
\]

Every Q26474 orientation contains an external wall touching only the
orthogonal A1 root, with the same dotted square. The relevant inverse
Gram diagonal entries are3/4 and1/2. Their squared timelike heights are

\[
h_{A_3}=\frac{3\tau}{4}-2=\frac{5+3\sqrt5}{8},
\qquad
h_{A_1}=\frac{\tau}{2}-2=\frac{\sqrt5-1}{4}.
\]

The two positive projections lie in orthogonal finite components.
Equation(3) gives their actual wall product

\[
-\sqrt{h_{A_3}h_{A_1}}
=-\sqrt{\frac{5+\sqrt5}{16}}\in(-1,0).
\tag{4}
\]

Distinct normalized norm2 Coxeter wall normals cannot have product
strictly between−1 and0: orthogonal walls give0, a nonright finite
Coxeter angle gives at most−1, and disjoint or parallel walls give
at most−2. This excludes all eight cross-class pairs. The exclusion
uses actual hyperbolic geometry and the specified exact contacts;
it needs no additional arithmetic integrality hypothesis.

Together the two arguments prove at most one selected corner at any
Y vertex. Counting all such corners yields(1).

## Fixed-witness repair and scope

The corresponding new positive charge column has vertex cost1 at3137
and face coefficients[[26448,1],[26474,2]]. The exact diagnostic changes
only

\[
\mu_{3137}:\frac{200076613}{10^{10}}
\longmapsto\frac{310256599}{10^{10}}.
\]

All face masses are fixed. Every inherited signed or cap load is
unchanged, every complete local packing cover remains valid, and every
old local/joint/forced-incidence inequality has unchanged load and a
nondecreasing allowance. Thus the old checks persist without repeating
the94722 signed checks or the millions of joint subset checks. The new
normalizer, exact repaired ratio, and its comparison with12 are saved
as fractions. This is a valid repair of the existing witness, not a
search for the best witness after adding(1).

**FINAL/FROZEN:** final generation and fresh saved replay passed. Root's
independent source/proof read accepted the common-frame formulas,
automorphism and provenance accounting, the15 pair exclusions, and the
monotone repair; it did not duplicate the arithmetic replay. A separate
agent's proof/scope read accepted the fixed-sign argument, exact height
product, shared-quad uniqueness, and occurrence-versus-orientation
distinction; it did not reconstruct the15 pairs or provenance maps.
Direct input hashes, source syntax, and local report links were checked.
