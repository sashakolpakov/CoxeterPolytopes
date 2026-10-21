# Universal order12 costs B62 and C62

The affine order12 family maximum falls from63 to **62**: B63→62,
C63→62, with A60 and D62 unchanged. The B change halves two distance9
pairs; the C change deletes one distance9 pair. The only new finite
payment is a half between a whole A1 and the penultimate node adjacent
to the doubled terminal of a separate whole B11. There is exactly one
eligible marked pair, costing **1/2 per A1/B11 component pair**. This fits
the already counted A1/D11 representative allowance5/2. Every triangle
still receives at least2 and every quadrilateral at least1. No new
ambient dimension or finite forest maximum is asserted.

The [producer](../scripts/fresh_order12_bc_near_terminal_half_transfer_20261018.py)
and [certificate](fresh_order12_bc_near_terminal_half_transfer_20261018.json)
prove universal profile identities and whole-component alternatives;
they do not use a bounded overlap census as a surrogate for arbitrary
supporting roots. The starting rules are the
[B12 center-zero theorem](fresh_order12_b_center_terminal_transfer_20261018.md)
and [B/C12 unit-coupon theorem](fresh_order12_bc_unit_coupon_transfer_20261017.md).
The [BC11 zero theorem](fresh_order11_bc_zero_transfer_20261018.md),
[coupled order11 theorem](fresh_order11_bcd_half_pair_transfer_20261018.md)
and [D13 center/leaf proof](fresh_order13_d_center_leaf_half_transfer_20261018.md)
supply the established profile and independent-priority interfaces.

## Marked pairs and actual quadrilateral

For B11/order12 use spine0–…–9, doubled bond8–9 and fork leaves10,11 at0.
Halve(8,10),(8,11), joining the node just before the doubled terminal to
either opposite fork leaf. For C11/order12 use path0–…–11 with doubled
end bonds0–1 and10–11. Delete(1,10), the pair just inside the doubled
terminals. All three pairs have distance9 and the changed sets are
invariant under every native diagram automorphism. All short8 weights
are unchanged, including previously zero weights.

At an affected ideal corner let F be the whole varying affine factor
containing consecutive face marks a,b. Write the quadrilateral star as
W+a+b+c+d with opposite pairs a,c and b,d. The rest of that corner is
a positive-definite spectator diagram Z orthogonal to F. All physical
root squares are2. Supporting-white and white-plus-one-black principals
are strictly positive definite; corner principals are positive
semidefinite. These are the inherited actual factor-fiber hypotheses,
including nonsimple ideal vertices.

A nonzero adjacent black product gives a separate distance1 unit. If
c,d share a nonzero contact coordinate in Z, their corner gives a
distance2 unit. Otherwise adjacent products vanish and the nonnegative
Z-contact vectors x,y have disjoint supports. Their opposite products
with a,b remain arbitrary native magnitudes, including parallel or
dotted products. A nonzero magnitude is at least1, and above1 the next
possibility is sqrt2. Heights h,k against the positive cusp kernel are
strictly positive.

## Complete native profiles

At a node just before a doubled terminal, deleting the mark leaves an
isolated A1 and a large finite component. For C this large component is
B10; for B it is D10. At the B fork-leaf mark, deletion leaves B11.
For each large component the inverse diagonal is1 at its relevant long
endpoint,2 at the next node and strictly greater than2 elsewhere. These
complete lists are checked exactly.

A positive-semidefinite one-root extension of that connected component
has no contact, a simple or doubled endpoint contact, or a simple
next-node contact. Two contacts would create a cycle in a diagram
already containing a fork or a doubled bond. Finite diagrams are trees;
the only cyclic affine type is A, which cannot contain this old fork or
doubled bond. Exceptional finite/affine types are too small. The single
contact Schur condition q²J_ii≤2 gives the stated all-label list, including
equality; the next native magnitude after sqrt2 already has square>2.
Additional spectators cannot repair a non-PSD principal extension.

For a near-terminal mark m, outer terminal o, first inward node e and
second inward node f, the complete profiles are:

|Profile|Contacts|Height h|Moment l|c0|Cost upper bound|
|---|---|---|---|---|---:|
|isolated|m:t|t|0|2|2|
|tail only|m:t,o:q|t+q/sqrt2|0|2−q²/2|3|
|end3|m:t,e:1|t+1+q/sqrt2|1|1−q²/2|3|
|end4|m:t,e:sqrt2|t+sqrt2|sqrt2|0|2|
|next3|m:t,f:1|t+1|2|0|4|

In the end3 row the optional additional outer contact is o:q. The
isolated profile has t≥1. Tail only has1≤q<2 and t≥0, by the strict
edge principal on o and the new root. End3 has q=0,1 or sqrt2: the
extended large block has remaining Schur norm1, so q²/2≤1. End4 and
next3 contain an affine principal with positive kernel, forcing q=0
and sealing it against every other supporting white in that corner.
Their opposite t-contact remains allowed because the omitted mark is
not in that corner.

For the B fork end use b=10, mate11 and center0. The complete four
profiles are:

|Profile|Contacts|Height k|Moment m|d0|Cost upper bound|
|---|---|---|---|---|---:|
|isolated|b:t|t/2|0|2|8|
|mate3|b:t,11:1|(t+1)/2|1|1|8|
|mate4|b:t,11:sqrt2|(t+sqrt2)/2|sqrt2|0|4|
|center3|b:t,0:1|1+t/2|2|0|4|

The first two attain8 only at t=1 and t=0 respectively, with k=1/2.
All other native instances have cost at most4. Polynomial certificates
and the native gap establish these bounds for unbounded t; no numerical
label cutoff is imposed.

Give the B spine cusp coefficient1, its doubled terminal1/sqrt2 and
fork leaves1/2. Give the C inner spine coefficient1 and both doubled
terminals1/sqrt2. Exact elimination of F minus the two marks yields
resistance R=10 for B and R=9 for C. With end costs

    g=2l/h+c0/h²,   f=2m/k+d0/k²,

eliminating the finite support gives the checked identity

    det(S)/(h²k²) = 1−(g+f)/R + E/R,
    E=(x/h−y/k)^T G_Z^(-1)(x/h−y/k) ≥0.

The actual four-black Schur determinant is zero. The producer verifies
all19 end identities and45 cross identities, both principal inverse
lists per family and the full formal determinant including spectator
terms and unequal heights.

## C: a uniform strict gap gives an independent unit

For C use a=1,b=10. Both sides have the five near-terminal profiles.
Every cost is at most4, so the normalized determinant is at least
1−8/9=1/9. It cannot be zero. Therefore the no-adjacent-contact,
disjoint-spectator branch is impossible: the actual face already has
a distinct distance1 or2 unit. No original C weight or finite coupon
is needed, so its selected unit can be deleted completely.

This is an independent actual-corner theorem. Apply it before the new
B half argument if both changed types occur on one face.

## B: immediate affine alternatives

For B set a=8,b=10. If c has near-end profile end4, the adjacent corner
contains a whole affine B10 of order11 with its doubled-terminal and
opposite short-leaf marks. Its current literal weight is zero. Apply
the independent BC11 zero theorem first whenever such an actual corner
occurs; it supplies a full face payment through its own retained paths
and distinct finite alternatives. It does not consume a new B12 half.

For near-end next3, the adjacent whole factor is affine D10 of order11,
marked at opposite short leaves. Its current weight is1/2. That corner
is distinct from the selected B12 corner; their two retained halves
sum to1 directly. A positive kernel seals these affine components
against all additional supporting roots in their actual corners.
The source reconstructs both kernels, every canonical marked image
and the exact current B/D weights. It does not assume an old all-pairs
unit cutoff for these order11 factors.

## B: arbitrary spectators force a whole finite alternative

After those immediate profiles, the near-end cost satisfies g≤3.
If the fork side is nonextreme, f≤4, giving a strict determinant gap
at least3/10. Hence the fork side is extreme: a single simple contact
to b or its mate, with k=1/2 and f=8. The determinant gives

    E=g−2,   0≤E≤1.

The complete near-end possibilities with g≥2 are:

|Profile|Parameters|Height|Cost|
|---|---|---|---:|
|isolated|t=1|1|2|
|tail only|t=0,q=1|1/sqrt2|3|
|end3|t=0,q=0|1|3|

Indeed isolation has cost2/t². Tail only at t=0 has cost4/q²−1,
so the next native magnitude gives1; if t≥1 and q≥1, monotonicity
bounds its cost by the value at t=q=1, below2. End3 decreases in both
t and q on0≤q≤sqrt2; either t≥1 or q≥1 gives a value below2. Thus
**h≤1** in every remaining case, without bounding the spectator rank.

The entire fork-side vector y vanishes. If y≠0 and x=0, coordinate
Cauchy gives E≥(1/k)²/2=2. If both supports are nonempty, choose i in
the support of y and j in the disjoint support of x. Then

    E ≥ (y_i/k+x_j/h)²/(4−2(G_Z)_ij) > 9/8 > 1.

The numerator is at least(2+1)²; strict positive definiteness of the
two-node spectator principal gives (G_Z)_ij>−2 and the strict denominator
bound<8. Both cases contradict E≤1. No restriction on x is needed.

Choose the actual corner containing a,d and omitting b,c. Since y=0,
all Z roots are orthogonal to its displayed components; c and all its
arbitrary contacts are absent. If d attached to omitted b, d is a
whole A1 and the far mark a belongs to whole finite B11, exactly at
the node next to its unique doubled terminal. Their new half plus
the selected B12 half gives1. If d attached to the mate11, both marks
lie in one whole finite B12 path, at marked distance10; its unchanged
unit already pays. The source checks positive pivots, actual components,
marked degree and doubled-neighbor role in both finite alternatives.
These whole-component assertions hold for arbitrary additional spectators.

## Cost and simultaneous-rule compatibility

The A1 mark and the penultimate B11 mark are intrinsically unique.
There is one new marked half, hence total finite increment1/2. The
current actual A1/B11 charge is0. Its representative pair A1/D11 already
has allowance5/2; all33,489 simultaneous finite-type comparisons through
rank60 continue to dominate the actual rule with this one addition.
No representative weight or residual item changes. This is a pointwise
cost check, not a fresh finite-forest optimization or a global bound.
Other new coupons, such as A4/B9 or D5/D11, concern distinct type pairs
and require their own explicit simultaneous binder; this theorem does
not silently omit their charges.

The proof is applied to actual changed corners, with the independent
BC11 and new C12 unit alternatives before the B12 half branch. The
existing B12 center-zero theorem has different marks and outputs only
its retained D11 center/leaf unit, whole finite D12 unit or typed
A2/B9 and A2/D11 halves. The D13 refinements likewise use smaller
B/D11 roles and D10/D12 coupons. None uses the present B12 or C12
changed distance9 pair. The B/C13 and order14/15/16/17 refinements use
independent short-path alternatives or their stated finite components;
they do not receive a new unpaid output here.

The new B proof itself uses only a retained D11 half, whole finite B12
unit or the new whole A1/B11 half. No changed same-component rank15–17
payment, long affine order13–17 pair, or generic short8 unit is consumed.
All new short1/2 alternatives are literal units under the current rule.
Thus current deep8→short7 and deep10 finite7→short6 replacements remain
independent, and no older long-pair exclusion is called recursively.
If no new changed corner occurs, the previous whole-face proof applies
pointwise; otherwise the universal actual-corner proof above pays it.
This avoids the need to reinterpret any preferred selected-map witness.

Triangles retain their finite7/ideal5, protected ideal6/7 and whole
small-affine-order-at-most9 outputs. The three changed pairs have
order12 and distance9. Therefore triangle payment2 and quadrilateral
payment1 persist. The factor-fiber divisor12 gives new rank11 envelope
62/12=31/6. A new ambient comparison remains a separate obligation.

## Verification scope

Generation chunk `deda55` and separate fresh replay `d170eb` both
terminated exit0/PASS. Root and construction_audit independently read
the complete producer and report and the stated inherited physical
profile, determinant and Cauchy interfaces. Both mathematical reviews
passed: they checked complete native contact reduction, whole-component
alternatives, the distinct retained half, exact new coupon cost and
noncircular current priorities. Neither reviewer repeated the formula,
map, typed comparison, arithmetic or historical census computations. Exact
symbolic end/cross and determinant identities, native bounds, immediate
whole-affine maps, finite whole roles, complete66-pair weight tables,
automorphisms and finite pointwise allowances are checked. Historical
factor-fiber and finite/affine classification inputs are inherited;
no original geometric census or forest optimization is rerun.

    python3 scripts/fresh_order12_bc_near_terminal_half_transfer_20261018.py --verify
