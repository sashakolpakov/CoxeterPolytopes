# Independent local transfer for the nine dropped D16 pairs

The nine pairs between the two fork triples of the17-node affine D16 diagram can be removed in a literal shared-L25 H17 quadrilateral star if the stated finite half-payments are added. The transfer survives arbitrary additional supporting whites: either they give an immediate length-two payment, or they are orthogonal to the entire19-root core and preserve its finite component types.

The [independent exact checker](../scripts/fresh_h17_half_coupon_transfer_20261017_review.py) and [certificate](fresh_h17_half_coupon_transfer_20261017_review.json) verify the core signature, radical, affine positive kernel and all27 finite corners in the nine cases. They use Fraction elimination and literal graphs, without importing the proposed global producer's arithmetic. The actual small chamber is the frozen [known finite-volume H17 pyramid](fresh_nonincident_carrier_incidence_20261017.md); the transfer below concerns any actual quadrilateral with this same literal core and role assignment.

## Local hypotheses and finite coupons

Let the19-node core be two copies of L25 sharing node9, with the numbering of the [D16 bottleneck audit](fresh_h17_d16_payment_bottleneck_20261017.md). Its opposite black pairs are

`{a,3}` and `{b,13}`, with a in{0,1,2} and b in{10,11,12}.

All other core roots are supporting whites. Further supporting whites E are allowed, and their Gram is positive definite because the face has an interior hyperbolic point. The four displayed blacks are the actual four sides. The Gram on the displayed19 roots is assumed to be the literal core, including its zero contacts; obtaining that literal core from a general classification is a separate obligation.

Retain a unit finite payment when the two marks have distance at most9 or lie in a connected component of order at most17. Add weight1/2 at either of the following kinds of finite corner:

* the marks lie in distinct A1 and D16 components, with the D16 mark at its fork center or either short fork leaf;
* the marks lie in distinct A2 and D14 components, with either A2 mark allowed and the D14 mark at its fork center or either short fork leaf.

The component names here are finite types: D16 has16 roots and D14 has14 roots. These half-payments are disjoint from the retained finite rule, since their marks lie in different components. This note proves the required local transfer, not a universal bound on the sum of all such coupons.

## Extra supporting whites cannot change the required components

Write F for the17-node affine D16 subdiagram obtained by omitting3 and13 from the core. The actual corner retaining a,b contains F. It is necessarily ideal, since a finite corner cannot contain a singular principal diagram. Its positive affine kernel has coefficients1 at the four leaves and2 elsewhere. Every supporting white e in E lies at that corner. The positive kernel and nonpositive Coxeter contacts force e to be orthogonal to every root of F.

The exact19-root Gram has signature(17,1,1) and radical

`(-6,-3,-4,-2,-5,-4,-3,-2,-1,0,6,3,4,2,5,4,3,2,1)`.

This radical is an actual linear relation in the ambient Lorentz realization. Indeed, its represented vector is null and orthogonal to the indefinite core, hence to a timelike vector; a timelike vector's orthogonal complement is positive definite, so the represented vector must be zero.

Since e is already orthogonal to F, pairing that relation with e leaves

`-2(e·n3)+2(e·n13)=0.`

If these equal contacts are nonzero, the actual corner with marks3,13 contains the white-interior path3–e–13. That corner has a unit payment under the finite distance9 or ideal distance7 rules. At an ideal corner the actual factor-fiber lemma places this path in the varying affine factor. Thus this case already pays the quadrilateral.

Otherwise every e in E is orthogonal to all19 core roots. The finite corner stars then split as the displayed finite core corner plus the positive-definite diagram E. They remain finite, and the marked core components cannot be enlarged or joined by E. Consequently their literal A1/D16 or A2/D14 coupon types, and their existing connected-component payments, persist unchanged.

This argument concerns the supporting whites in the actual face star. It does not assume a simple ideal vertex or claim that all ambient facets incident to the ideal corner are among them. Those extra ideal facets are irrelevant to the finite coupon corners.

## All nine literal cases

In F, the fork triples are{0,1,2} and{10,11,12}. Their nine cross-pairs are exactly the eight pairs at distance at least13 plus the pair of fork centers at distance12. The checker verifies that F has corank1 and the other three corner stars are positive definite of order17, for every case.

| Dropped ideal pair | Old finite payment total | Added coupon total |
|---|---:|---:|
|0,10|0|1|
|0,11|1|1/2|
|0,12|0|1|
|1,10|1|1/2|
|1,11|3|0|
|1,12|1|1/2|
|2,10|0|1|
|2,11|1|1/2|
|2,12|0|1|

Only four cases need coupons. Their two half-paid marked pairs are:

| Dropped pair | First finite corner | Second finite corner |
|---|---|---|
|0,10|0 in D14 fork center;13 in A2|3 in A2;10 in D14 fork center|
|0,12|0 in D16 fork center;13 in A1|3 in A2;12 in D14 short leaf|
|2,10|2 in D14 short leaf;13 in A2|3 in A1;10 in D16 fork center|
|2,12|2 in D16 short leaf;13 in A1|3 in A1;12 in D16 short leaf|

Each required pair of half-payments sums to1. Together with the sealing argument, this proves the local transfer for every actual quadrilateral satisfying the explicit core hypotheses, with arbitrary additional supporting whites.

## Verification and limits

Generation and fresh saved replay both terminated successfully. Syntax, all7 direct source hashes and all4 local report links pass. The independent finite arithmetic covers one19-root signature/radical, one17-node affine positive kernel, nine15-node white supports and all27 finite corner matrices. It also checks the literal marked component roles and each face payment. A separate internal written-proof review passed the actual-corner affine-kernel argument, actual radical relation, length-two alternative and finite-component sealing; that review did not duplicate this arithmetic or any global embedding census.

This package does not reproduce the full old108 embedding census, prove that these are the only affected cases in that census, establish the global finite coupon budget or the all-rank ideal budget, or compose a new dimension theorem. Those are distinct obligations in the proposed parent argument. No frozen predecessor was edited. Source, data and report are final/frozen.
