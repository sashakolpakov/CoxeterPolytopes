# A D25 short-eight family financed at an off-face A16/E8 cusp

21 October 2026. New manual mathematical supporting result, pending
full root and outside-author review. No scientific, symbolic, enumeration,
Lean or TeX job ran. No manuscript, current rule or dimension bound changed.
The author supplied the earlier A15-cusp and context662 financing proofs;
those are inherited inputs, not independent reviews of this result.

**Result.** The explicit full-star family below has two D25 short-leaf
distance-eight units and no other current corner payment. Deleting those
units cannot be justified by a local replacement corner. Nevertheless each
actual face in this family determines an ideal vertex OUTSIDE the face,
with whole affine A16 and E8 factors. An injective link-rectangle allocation
charges one replacement unit per face to these cusps. The actual A16 cycle
saving finances every such allocation, with arbitrary other cusp factors.
This discharges this particular direct finite-eight producer family. It
does not classify all D25 short-eight quadrilaterals or authorize deleting
the entire entry, a new finite maximum, or a ceiling56.

## 1. Pinned inputs and precise family

Paths in the table are relative to `results/`, except the manuscript.

| Input | SHA256 |
| --- | --- |
| `finite57_context662_opposite_vertex_incidence_financing_20261021.md` | `21aab45e9020be06613301d2a5ab679e97407b43d01be5d0da84a5340d872cad` |
| `review_finite57_context662_opposite_vertex_incidence_financing_root_20261021.md` | `fd03ce8034a1aead74f0db636efdf0ad829a9efaaecafbb8dbe36e262883203d` |
| `review_finite57_context662_opposite_vertex_incidence_financing_outside_author_20261021.md` | `03adaf48658aa90ec223227e7074ac7c3f76028af4e8481c4c4687d99fbfe6a6` |
| `finite57_remaining_a17_d25_producers_and_rank16_financing_obstruction_20261021.md` | `80632655e62e39fa7fc60ef28506a6364908cdcb450a4cee2ba2bedb7d30b7d6` |
| `finite57_full_cusp_deletion_rigidity_and_curvature_financed_priorities_20261021.md` | `ed578a84178a982c988bebcad8d2c6212b650c3d2e48928b83622f9ad8b54c4c` |
| `fresh_cusp_factor_payment_framework_20261017.md` | `c65cfc5a4409a8cc83fcf0a083a6c1a7a29bb7bc7d649c11435773877e70ccab` |
| `order11_fractional_cusp_closure_research_20261018.md` | `3cdd8e76fb472d45a36ff4ad4ff75af17ffaf4cd0e7eaf965759566a3f46fb0c` |
| `../manuscripts/finite_volume57.tex` | `c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0` |

Work in a hypothetical finite-volume Coxeter57-polytope. Use norm-two
roots and the manuscript's D25 numbering: spine0..22, fork leaves23,24
at22. An actual quadrilateral Q has opposite side-root pairs {a,c},{b,d},
where a=15,b=23. Its complete supporting-white forest is

    W = (D25 minus {15,23}) + Z,       rank(Z)=32.

At its a,b corner D25 is WHOLE. Require the following EXACT contacts:

    c.0=c.15=c.16=-1,       d.23=-1;

all other c/d products with D25, with each other, and with Z are zero.
The Z roots are orthogonal to the entire27-root displayed core. These
are full actual-star hypotheses, not a claim that every D25 producer
has this form. Faces, rather than embeddings of the profile, are counted.

## 2. Exact local geometry and the remaining local payment

Let J be the inverse norm-two D25 Gram. For spine indices i,j, written
as one-based i+1,j+1, its entries are min(i+1,j+1); the spine/fork entries
are (i+1)/2, and J[23,23]=25/4. Put x=e0+e15+e16 and y=e23. Then

    x.J.x=70,       y.J.y=25/4,       x.J.y=17.

Eliminating D25 from its extension by c,d leaves

    [ -68    -17  ]
    [ -17   -17/4 ],

negative semidefinite of rank1. Thus the27-root core has signature
(25 positive,1 negative,1 zero). All four corner diagrams are finite:

| Corner marks | Complete finite components, apart from Z | Current payment |
| --- | --- | ---: |
| a,b | D25; marked distance8, b a short leaf | 1 |
| a,d | A24 + A1(d) | 0 |
| c,b | D25; marked distance8, b a short leaf | 1 |
| c,d | A24 + A1(d) | 0 |

In the second A24 path the order is14,13,...,0,c,16,...,22,24;
c occupies position15. The two zero payments are unlisted A1/A24
cross roles, and no presence guard changes them. Extra Z components
cannot change these whole marked factors under the stated orthogonality.
Neither marked D25 pair satisfies the deep-pair deletion or the
distance-nine short-leaf half condition. There is no white-interior
path of length<=7 between the two opposite black pairs: the shortest
ones are a to b and c to b, both length8; d has no white contact.

For completeness the Gram profile has genuinely feasible finite corner
points in a Lorentz realization; positive corner diagrams alone are not
being called an ambient completion. Let omega_i be the D25 duals and
t a future unit timelike vector orthogonal to them. Set

    c=-(omega0+omega15+omega16)+2sqrt(17)t,
    d=-omega23+(sqrt(17)/2)t.

Taking timelike coordinate1, the four corner rays have spatial parts
0; -2sqrt(17)omega23/25; -2sqrt(17)omega15/33;
and -(16omega15+2omega23)sqrt(17)/281. Their Lorentz squares are
-1,-8/25,-1/1089,-200/78961, respectively. The two unselected side
products are strictly negative at each ray: at a,d the c product is
-16sqrt(17)/25; at c,b the d product is -sqrt(17)/66; the remaining
ones follow directly from the displayed negative dual coefficients.
Thus this is an exact local obstruction to deleting the two units
without replacement, not merely an indefinite matrix with infeasible
corners. No finite-volume ambient realization is asserted.

## 3. An actual off-face ideal vertex forced by the relation

The roots0..15,c form an affine A16 cycle. The roots17..24,d form
affine E8: its branch22 has arms of lengths5,1,2. Their positive null
vectors satisfy in the Lorentz realization

    nu = sum_(i=0..15) alpha_i+c
       = 2(alpha17+2alpha18+3alpha19+4alpha20+5alpha21+6alpha22
           +4alpha23+3alpha24+2d).                    (1)

Both sides are future null, since they have negative product with the
finite a,b corner. Moreover nu.alpha16=-2, whereas every other core
root is orthogonal to nu. Hence nu is OUTSIDE Q: alpha16 supports Q.

For an ACTUAL Coxeter polytope all distinct wall products are
nonpositive. The positive expression on either side of(1) therefore
has nonpositive product with every actual facet normal. Consequently
the ideal ray nu lies in the closure of the actual polytope. The
finite-volume vertex theorem makes it an actual ideal vertex. Every
ADDITIONAL facet through nu, not among the two displayed affine
components, has zero product with each summand in both positive null
expressions. The displayed A16 and E8 components are mutually
orthogonal. Thus they are WHOLE factors of its cusp diagram, not
proper subdiagrams being charged as complete factors.

This is independent of the A15 source and context662. For example a
full-rank native local cusp completion with no A15 is obtained by
letting Z=4D8 and adding g_j=-theta_j+nu/2 for their highest roots.
Then g_j^2=2, g_j.alpha16=-1, and every other core product is zero;
the source cusp is affine A16+affine E8+4 affine D8, of rank56.
The new g_j are NOT supporting whites of Q. Their finite-corner
products are negative, because these corners are orthogonal to Z and
future timelike. This remains local cusp/face data, not a polytope.

## 4. Injection into cusp link rectangles

Drop the supporting facet alpha16 from W. The remaining54 supporting
facets determine the actual THREE-FACE R containing Q and nu. At nu
its link is the product of an edge from the A16 simplex factor and
an edge from the E8 simplex factor; all other simplex factors are at
vertices. Thus this is one specified rectangular two-face of the
product link. Its four link vertices specify four edges from nu.

Each corner of Q lies on the corresponding one of those four actual
edges. Indeed its segment to nu lies in the polytope by convexity,
and keeps the56 supporting facets other than alpha16. It is the
unique finite endpoint of that edge: no intervening truncation can
separate two points already in a convex face. Therefore two members
of the family mapped to the SAME nu and link rectangle have the same
four corner vertices, hence are the same two-face. This proves an
injection of actual quadrilaterals, not of core embeddings.

Let a be the number of whole affine A16 factors and e the number of
whole affine E8 factors at nu. Write V=product(r_i+1) for its link
vertex count. For each chosen A16/E8 factor pair, there are17 allowed
adjacent pairs in the A16 cycle. In the E8 graph the marked pair b,d
is the UNIQUE terminal edge on the length-two arm (the other arm
lengths are1 and5). Thus there is ONE allowed E8 pair, not all its
eight edges or all36 pairs. The remaining factors contribute their
ordinary vertex ports. The resulting count is

    number of charged Q at nu <= (a*e/9)V.             (2)

If a face has more than one description, choose one; the same bound
holds, because a fixed link rectangle still determines at most one Q.
Since16a+8e<=56 and a>=1 whenever a debt occurs, e<=5.

## 5. Curvature financing and exact scope

Use the current manuscript envelope c(r) and its inherited exact
ambient57 bound, as in ed578:

    kappa/V - sum c(r_i) >= -261/2662 > -1/10.

An actual affine A16 cycle has exactly85 retained distance<=5 pairs,
so normalized charge5. The manuscript c(16)=227/34 therefore overcounts
each such factor by57/34. (The895/136 in806 is an actual D16 charge,
not this all-family envelope.) Charge ONE unit for each Q to nu,
even though nu is not incident to Q. By(2), total ideal charge plus
these debts is at most

    V[sum c(r_i) - (57/34)a + a*e/9].

Hence its normalized curvature surplus is strictly greater than

    a(57/34-5/9)-1/10 >= 781/765 >0.                 (3)

Here V>=17*9=153, so the absolute surplus is >781/5, in particular
larger than the inherited uniform43/128. Cusps receiving no debt keep
their existing bound. The old A15-conditioned rank<=11 restoration can
also be retained: at a cusp having A15 its separate saving91/64 exceeds
its entire extra cost41/33, and these savings concern different factors.
No A15 saving has been reused to prove(3). The context662 debt goes to
its opposite finite vertex, not to this source, and is unchanged.

Delete the two D25 short-eight entries provisionally. Each Q in this
family now has zero corner payment, so credit it one unit and allocate
that unit to its uniquely selected source nu. Equations(2)--(3) bound
the TOTAL added charge, without pretending nu is a corner of Q or
counting its omitted roots as supporting whites. All finite vertex
charges only decrease. Thus this family has a complete replacement
in the global face-payment argument.

Other D25 short-eight profiles have NOT been exhausted. In particular
contacts of c,d with the original spectator forest, other native
attachments to A15/A8, and the other independent direct producers
require their own proof or classification. The present result identifies
and finances a genuinely different off-face source; it is not a full
new payment rule, a new forest calculation, or an ambient56 theorem.
