# Context-dependent cusp payments and an exact allocation criterion

The known finite-volume H17 chamber admits a complete context-dependent payment certificate with cost109 at its affine D16 cusp and614 at its E8×E8 cusp. It satisfies every triangle and quadrilateral demand while keeping the retained finite payments unchanged. The total ideal payment723 is minimal for those fixed finite payments. This bypasses the136 D16 cost forced by intrinsic-diagram invariance in the [preceding bottleneck audit](fresh_h17_d16_payment_bottleneck_20261017.md).

The [checker](../scripts/fresh_h17_contextual_cusp_payments_20261017.py) and [certificate](fresh_h17_contextual_cusp_payments_20261017.json) cover all3,456 actual two-faces. The result is a certificate in this fixed chamber, together with a general exact criterion for allocating face deficits among cusps. It supplies no universal bound on those deficits and no new ambient dimension bound.

## An allocation lemma for actual face incidences

Let F range over the triangular and quadrilateral faces of a polytope. Set d(F)=2 for a triangle and1 for a quadrilateral. Fix nonnegative payments at its finite corners, and let

`delta(F)=max(0, d(F) - sum of its fixed finite payments)`.

Choose a set I(F) of eligible actual ideal corners, identified with their cusp vertices. A convex polygon has at most one corner at a given vertex. Payments at the other ideal corners are zero. We seek weights w(F,v) in[0,1] such that

`sum_{v in I(F)} w(F,v) >= delta(F)`

and, for each cusp v, their total is at most a prescribed budget B(v). No diagram invariance or transitivity at nonsimple cusps is assumed.

**Allocation criterion.** Such payments exist if and only if, for every subset U of the cusp vertices,

`sum_F max(0, delta(F) - |I(F) outside U|) <= sum_{v in U} B(v).`  (1)

This includes U empty, which requires delta(F)<=|I(F)| for every face. If all deficits and budgets are integral, integral weights0 or1 suffice.

To prove the criterion, build a finite network with edges from the source to each face of capacity delta(F), from that face to each eligible cusp of capacity1, and from each cusp to the sink of capacity B(v). A complete flow supplies the payments. It may be assumed that each face receives exactly its deficit, since any excess payment can be removed.

Fix the cusp vertices U on the source side of a cut. The face F should be placed on that side precisely when delta(F)>|I(F) outside U|. Minimizing over the face choices gives the cut capacity

`sum_F delta(F) - sum_F max(0,delta(F)-|I(F) outside U|) + sum_{v in U}B(v).`

Thus all cuts admit the required total flow exactly when(1) holds. For completeness, the flow/cut implication needs no geometric hypothesis: a maximum flow exists by compactness of the bounded finite flow polytope; if its residual graph has a source-to-sink path it can be increased, while otherwise the reachable set gives a cut whose capacity equals its value. For integral capacities, the same augmentations can preserve integrality. This proves the assertion, including real nonnegative budgets.

The criterion is not by itself a curvature estimate. Proving(1) for useful geometric budgets is the remaining substantive obligation. Its left side measures actual unpaid face demand concentrated on a cusp subset, rather than all faces having a nonincident affine carrier.

## Complete H17 reconstruction

The chamber and its complete cone are the frozen [known H17 pyramid](fresh_nonincident_carrier_incidence_20261017.md), with80 finite vertices and two proper ideal vertices. The apex has affine factors E8×E8; the other cusp has affine D16, a17-node diagram. This package does not repeat the original chamber or volume proof.

At each finite corner retain weight1 if its marks have distance at most9 **or** lie in one connected component of order at most17; otherwise retain weight0. Every actual ideal corner is eligible. The source reconstructs the marked pair and its component in the actual face-corner star, not merely the full diagram of a nonsimple cusp. It verifies every face demand and independently accumulates the finite-vertex pair costs. Each of the80 finite vertices is incident to136 two-faces, and its paid cost equals the sum of binomial(component order,2), at most136<=9·17.

The complete deficit ledger is:

| Face type | Ideal corners | Required ideal payment | Number of faces |
|---|---:|---:|---:|
|Triangle|0|0|1,456|
|Triangle|1|0|44|
|Triangle|1|1|644|
|Triangle|2|1|2|
|Triangle|2|2|14|
|Quadrilateral|0|0|1,232|
|Quadrilateral|1|0|15|
|Quadrilateral|1|1|49|

The single-cusp demands split as598 at the E8×E8 apex and95 at D16. Fourteen triangles through both cusps require both unit payments. Two further triangles require only one payment, which can be assigned to either cusp. Therefore the exact feasible budget region is

`B_D16 >= 109,    B_E8xE8 >= 612,    B_D16+B_E8xE8 >= 723.`  (2)

The checker saves every nonzero term of all four cusp-subset inequalities, including the empty set. It also supplies three explicit integral payment certificates:

| D16 cost | E8×E8 cost | Total ideal payment |
|---:|---:|---:|
|109|614|723|
|110|613|723|
|111|612|723|

For each certificate every corner weight is0 or1 and every triangle/quad inequality is checked separately. The first allocation pays exactly the109 D16 pairs from the full-baseline bottleneck audit, but it drops34 unnecessary E8 apex payments as well. All723 ideal units are necessary in total because they equal the sum of the fixed face deficits.

For comparison with factor-normalized budgets, the D16 cusp has sV=17 and the E8×E8 cusp has sV=2·81=162, where V counts link vertices and s counts positive-dimensional simplex factors. In this fixed chamber the least common coefficient c allowing budgets B(v)=c·s(v)V(v) is exactly109/17: it is the largest of the three nonempty-subset ratios in(1), and the first explicit certificate fits those budgets. This number is **not** a universal coefficient or dimension result.

## What this resolves, and what it does not

The earlier136 obstruction concerned weights invariant under the intrinsic D16 diagram automorphisms, with the other cusp payments fixed. The present weights depend on the actual face and the payments available at its other corners. They preserve the finite rule and satisfy all actual faces simultaneously. Hence the H17 example does not force cost136 for such contextual rules.

There are only two flexible cross-cusp triangles here; the14 two-cusp triangles with zero finite payment need both unit payments under the[0,1] restriction. Transferring more mass to the E8 apex would require allowing some individual ideal corner weights above1 and a correspondingly different capacity criterion. No such change is used in the certificate.

The frozen nonincident-carrier tests remain binding cautions. The same chamber has1,048 unpaid quadrilaterals under the much finer Vinberg rule, defeating the proposed link-edge bound and the coefficient6 small-factor budget. A separate finite-volume H3 family has unbounded numbers of carried but fully paid compact quadrilaterals at a fixed cusp-link type. These are different payment rules and different quantities; none is contradicted by the present certificate. We do not charge every face carrying E8 to that cusp, or infer a multiplicity bound merely from a carrier's presence.

A future universal argument must bound the exact left side of(1), for all relevant cusp subsets, in terms of their allowed geometric budgets. The present package provides a tested allocation formulation and a complete positive example beyond intrinsic-diagram invariance; that global inequality remains open.

## Verification scope

Generation and fresh saved replay both terminated successfully. The finite computations reconstruct all3,456 face deficits,80 finite-vertex costs, three complete allocations and the four exact Hall cuts. Syntax, all7 direct source hashes and all4 local report links pass. Separate internal reviews passed the full source/mathematical argument and the allocation theorem/budget derivation respectively; neither review duplicated the full face-ledger arithmetic or the original volume proof. Source bindings preserve the original cone/face ledger and the full-baseline bottleneck audit. No new geometric existence statement, optimizer or ambient dimension claim is part of this calculation. Source, data and report are final/frozen.
