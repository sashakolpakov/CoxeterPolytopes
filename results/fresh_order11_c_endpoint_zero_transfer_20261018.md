# The order11 affine C endpoint half can be deleted

In the whole varying affine C10 factor of order11, assign weight0 to its
two doubled terminal leaves instead of the previous1/2. Every affected
quadrilateral still receives at least1 at another corner. The C total
decreases from109/2 to54, with no new finite coupon and no change to any
other weight. Combined with the completed order11 rules, the family totals
are A99/2, B54, C54, D53, so the normalized rank10 envelope is54/11.

The [source](../scripts/fresh_order11_c_endpoint_zero_transfer_20261018.py)
and [certificate](fresh_order11_c_endpoint_zero_transfer_20261018.json)
check the full native profile calculation. This is a universal local
theorem, with arbitrary supporting spectators and all native labels.
It is compatible with the [completed dimension63 rule](fresh_order11_order17_dimension63_20261018.md),
including cycle11 halves, B/D11 halves, the order17 refinement, finite
depth11 zeros and the D16/D17 exchange. No new ambient bound is claimed.

## Complete physical end profiles

Use an actual quadrilateral star W+a+b+c+d with opposite pairs{a,c},{b,d}.
At the affected ideal corner let F be the whole varying C10 factor, with
path0,...,10 and doubled bonds0--1 and9--10. Its marks are a=0,b=10.
The remaining supporting diagram Z is positive definite and orthogonal
to F, by the [actual factor-fiber theorem](fresh_cusp_factor_payment_framework_20261017.md).
All physical roots have squared norm2. The positive kernel nu has
coefficients1/sqrt2 at0,10 and1 at1,...,9. The opposite heights
h=−c·nu and k=−d·nu are strictly positive.

If either c.b, c.d or d.a is nonzero, the corresponding actual adjacent
corner has a marked path of length1. If c,d have a common nonzero contact
coordinate in Z, their actual corner has a marked path of length2.
Every such pair retains a full unit in the completed rule and the new
rule. These are corners different from the chosen C endpoint corner.
Thus suppose these alternatives are absent. The three adjacent-black
products then vanish, and the nonnegative native contact vectors x,y
to Z have disjoint supports.

The complete one-root extension classification from the
[coupled order11 theorem](fresh_order11_bcd_half_pair_transfer_20261018.md)
applies without discarding any affine alternatives. Deleting a or b
leaves a whole displayed B10 principal diagram. Its norm-two inverse
diagonals are1 at the simple endpoint near the omitted mark,2 at its next
node, and greater than2 everywhere else. A one-contact PSD extension
therefore allows a label3 or4 at that endpoint, or a label3 at the next
node. Two old contacts create a cycle, impossible in a finite diagram or
in a classical affine diagram containing the old doubled bond. Exceptional
finite and affine diagrams are too small to contain B10. The no-contact
profile is also allowed. These are all possibilities, not a label-limited
census. Spectators cannot repair a non-PSD principal extension.

The magnitude t of the contact with the omitted opposite old mark is
unrestricted except for its native Coxeter values: t=0, or t>=1, including
parallel and dotted possibilities. The four profiles have the following
exact physical quantities, with end cost g=2l/h+c0/h².

|Profile on the old B10|h|l|c0|Native range|Maximum g|
|---|---|---|---:|---|---:|
|isolated|t/sqrt2|0|2|t>=1|4|
|endpoint3|t/sqrt2+1|1|1|t>=0|3|
|endpoint4|t/sqrt2+sqrt2|sqrt2|0|t>=0|2|
|next-node3|t/sqrt2+1|2|0|t>=0|4|

For the isolated profile t>0 follows from positive height, hence t>=1.
In every other profile positive height is automatic at t=0. All four
cost functions decrease with t. The exact source also substitutes
t=u+t_min and checks that maximum−g has a polynomial numerator with
nonnegative coefficients and a strictly positive polynomial denominator
on u>=0. Thus the stated maxima hold on the entire continuous rays,
which contain every native label and every dotted opposite product.

The source independently constructs the norm-two C10 Gram, its positive
kernel, the two B10 inverses and their only admissible contact positions.
It rederives every height, moment and norm constant, and all16 cross
identities. No project arithmetic helper is imported.

## A strict contradiction without spending an affine half

Eliminating the nine white roots in F and scaling the old marks by their
kernel coefficients leaves

    (1/10) [[1,−1],[−1,1]].

The exact quadrilateral Schur determinant is the inherited identity

    0 = 1−(g_c+g_d)/10
          +(x/h−y/k)^T G_Z^(-1)(x/h−y/k)/10.           (1)

The source rechecks its symbolic determinant directly. As in the actual
factor-fiber framework, elimination of all supporting whites leaves the
four side roots of a three-dimensional face quotient, so the determinant
vanishes. The last term is nonnegative because Z is positive definite.

Every profile in the complete table has g<=4. Therefore the right side
of(1) is at least

    1−8/10 = 1/5 >0,

a contradiction. The source checks the16 ordered profile-pair bounds;
their uniform minimum is exactly1/5. The spectator form was only bounded
below by0, so this proof neither bounds spectator rank nor assumes the
spectators empty or sealed.

In particular, endpoint4 and next-node3 are kept inside the strict
calculation. The previous half theorem could use their other whole affine
corners as second halves. That step is not used here: their costs are at
most2 and4 respectively, so they cannot evade the contradiction even if
their new affine endpoint payment has also been deleted.

Consequently at least one of the earlier direct alternatives must occur:
an adjacent black contact or a shared spectator coordinate. A distinct
actual corner then contributes a retained full unit at distance1 or2.
The original endpoint weight is unnecessary and may be0. This proves
quadrilateral coverage for every affected actual C corner, without a
selected-core census or a completion assumption about other finite factors.

## Simultaneous B/D, triangle and later interfaces

The new pair set consists of the unique terminal pair of the whole native
C factor and is invariant under its reflection automorphism. It has
distance10. All54 other pairs remain units, giving total54.

The simultaneous proof uses actual corners, rather than attempting to
reuse the old two-half argument with a missing half:

1. If a quadrilateral has any newly changed C endpoint corner, the proof
   above supplies a full distance1 or2 unit at a different actual corner.
   This unit is unchanged by every completed finite or ideal refinement.
2. Otherwise every actual corner has exactly its completed baseline
   weight. The baseline face-payment theorem supplies its required total.

This directly handles the apparent coupling with the B/D order11 theorem.
That proof sometimes outputs an actual C endpoint half. If such a corner
occurs after its deletion, it lies in the first case and the new independent
unit alternative pays the face. If no such actual corner exists, the old
coupled proof's actual total is unchanged. There is no recursion between
a missing C half and a B/D half, and no inference that the old half theorem
alone justified weight0.

The [finite7/ideal5 triangle theorem](fresh_finite7_ideal5_triangle_payments_20261017.md)
uses finite paths at most7, ideal paths at most5 or protected ideal paths6/7,
and whole affine factors of order at most9. The changed distance10 pair
in an order11 factor is none of those obligations. Triangle coverage
therefore remains at least2. The cycle11 triangle exclusion is unchanged.

The only new terminal payments used here have distance1 or2. They survive
the finite depth11 zero rule, all finite D16/D17 reductions, the order17
ideal half refinement, and every cycle11 or B/D11 half. No whole finite
component or separate-component coupon is consumed. Every finite cost,
including both enlarged coupons in the completed dimension63 budget,
stays exactly as before. The order11 all-family envelope improves to54/11;
the other still-open ambient63 ideal and finite obligations remain separate.

## Verification scope

Generation08ec1b and separate fresh saved replay555738 both terminated
exit0/PASS. Direct audit7e7778 passed all9 source bindings, all6 report
links and source syntax. The parent completed a full source/report
mathematical read, passing the complete B10 principal extension domain,
the exact polynomial native-ray bounds for all four profiles, reflected
cross identities, nonnegative spectator correction and changed-C-first
actual-corner priority. It also passed the triangle and unchanged
distance1/2 interfaces and absence of a new finite coupon. That review
did not rerun the matrix/scalar computations or independently reprove
the inherited finite/affine classification. No concrete gap was found.
The checker covers both deleted B10 inverses, all four exact physical
profiles, every native-ray upper bound, all16 cross/profile-pair bounds,
the symbolic Schur identity, the changed orbit and exact weight total.
The all-label completeness and arbitrary-support interpretation are the
written finite/affine classification and actual-face argument above.
No finite forest optimization, global scalar comparison, external
refereeing, formal verification or new polytope is asserted here.
Source, certificate and report are final/frozen. A requested additional
peer read is not counted as completed in this report.

    python3 scripts/fresh_order11_c_endpoint_zero_transfer_20261018.py
    python3 scripts/fresh_order11_c_endpoint_zero_transfer_20261018.py --verify
