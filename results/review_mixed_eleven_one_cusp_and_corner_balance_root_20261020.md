# Root review: one-cusp exclusion and balanced dimension-eleven counts

20 October 2026. **PASS as written local and necessary-count theorems.**
This review does not exclude dimension eleven. It records a complete
independent read and algebraic derivation, with the explicit inherited
geometric inputs below. No Lean or geometric enumeration was performed.

Reviewed new sources:

* `right_angled_one_cusp_ten_facet_exclusion_research_20261020.md`,
  SHA256 `19590da4197ee04940df141a492ff89d69aea44b49bf87fe3c63699c252b56d4`;
* `right_angled_mixed_eleven_three_face_corner_balance_20261020.md`,
  SHA256 `73b8073843171d4bf11ad0fc11829b0855df36b2a78307fc024cd6d26dd450df`.

## The one-cusp argument

The area and Euler identities force the hypothetical ten-facet equality
case to have four cusp quadrilaterals and six finite pentagons. The cusp
rectangle makes opposite quadrilateral planes disjoint at finite points;
adjacent quadrilaterals meet precisely in their cusp edge. Thus the
eight finite vertices in the displayed boundary cycle really are distinct.
Finite trivalence supplies the four pentagon labels on that boundary.

The proof that opposite labels cannot coincide is essential and correct:
two convex facets have only one common edge. Repeating a label would put
both successive boundary vertices on that edge and force the intervening
pentagon boundary to be a triangle. Hence there are four distinct boundary
pentagons and four distinct cyclic pentagon adjacencies.

The remaining two pentagons have no quadrilateral neighbors. Each has five
edges and can meet each other pentagon along at most one edge, so each
meets all five. Their nine different adjacencies plus the four boundary
adjacencies exceed the eleven available pentagon/pentagon edges. There is
no implicit assumption that an arbitrary dual boundary is simple.

This proves the stated weaker bound F>=11 from first principles. It is
not claimed as a new best three-dimensional classification: the already
identified Nonaka theorem gives the stronger F>=12 for at most one cusp.
That published strengthening is not used in this proof or its budget.

## The exceptional budget and small types

Using p=4F-2nu-4, I get the exact budget difference
20nu-144-8F. The exceptional nu=10 type has F=7 and gives equality.
At nu=11 the newly excluded one-cusp F=10 case leaves F<=9. At
nu>=12 the elementary F<=nu bound suffices. On Q8 the budget is
zero, and P3 is excluded. The compact dodecahedron verifies optimality
of coefficient18 when coefficient4 of the signed defect is held fixed.

I reread the complete small-type and degree proofs in
`right_angled_mixed_nu10_edge_transition_research_20261018.md`, including
the attachment lemma's role in forcing R7. The counts used here are
P3:A3=2, Q8:A1=8, R7:A0=1 and A2=3. In particular R7's center and
three leaves cancel in A2-3A0; omitting the center would invalidate
the new balancing argument.

## The two global errors

For nu=11 the new lemma gives p<=10; for nu>=12,
p<=2nu-4<=10(nu-10). Consequently D>=-30E is valid without an
assumption on finite/finite edges. The earlier potential proof also
establishes this error bound; it is not newly discovered here.

I reread the full coefficient-two successor
`right_angled_mixed_nu12_carrier_error_sharpening_research_20261019.md`.
Its nu12 argument excludes three all-ideal-edge corners by K3,3 or
the ideal-edge capacity plus the bad-triangle attachment lemma. Its
nu13 argument excludes five by the same cusp capacity. For nu>=14
the elementary A3<=(nu-2)/2 bound suffices. With the separately reviewed
nu<=11 P3/D3 recognition input, this gives H<=2E. The D3 contribution
is two against one unit of E, and is not accidentally discarded.

## The balanced inequality and elimination

Put U=C(11-k,3), V=C(k,2)(11-k), W=C(k,3). At a finite corner,
X>=U+V, D=V-3U, H=W-t. In the claimed inequality the terms
in t cancel exactly. Its remaining check is

    40U+28V-16(165-W)>=0.

The complete twelve values in the source are correct. Summing uses
sum t=2B, and 3D-16H>=-122E. Thus

    31X>=16(165f-2B)-122E,
    680E+124R>=2640f-32B.

I independently eliminated a2,a3 to get

    R=-5B/4-165f/8+5760c-3E/4,
    360c<B-75f/2-E.

The first substitution gives 587E-123B-10395f/2+714240c>=0.
The second is strict, and 714240/360=1984 is positive. It therefore
gives exactly

    3722B>159195f+2794E.

Both strictness and the sign of the cusp substitution are correct.
The two cross-product comparisons are 257580 and 656. The unbalanced
new row has the larger excess coefficient2254/1379, so the two new
rows should be retained together.

An upper bound 2B<=84f+(3/2)E would contradict the balanced row:
159195/1861-84=2871/1861>0 and
2794/1861-3/2=5/3722>0. That proposed upper bound remains unproved.
The current unconditional mixed ceiling remains11.

This is internal written review of the stated sources and interfaces,
not an external review, a formal proof of their geometry or a claim of
worldwide priority.
