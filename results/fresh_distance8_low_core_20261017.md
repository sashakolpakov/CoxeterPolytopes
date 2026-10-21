# Complete small-core domain at distance eight

**Final, internally checked.** Exact generation and a fresh saved replay passed, and a separate narrow mathematical/source review found no gap in the stated scope. This is a local face-payment lemma, with no new ambient dimension claim.

Let H=W+U+V be the star of an actual quadrilateral in a finite-volume hyperbolic Coxeter polytope. W is its set of supporting white roots; U,V are the opposite black pairs. Normalize every root to square two. Suppose there is no path of length at most eight from U to V with white internal nodes. If a minimal indefinite core in W+U or W+V has order at most five, the configuration is impossible or has an adjacent marked corner whose **entire connected component has order at most12**, even after all additional supporting whites are restored.

The [source](../scripts/fresh_distance8_low_core_20261017.py) and [certificate](fresh_distance8_low_core_20261017.json) use the unchanged complete order-four/five graphs in the [distance-nine low-core certificate](fresh_distance9_low_core_20261017.md), and the high-core catalogue and actual-face/cusp conventions of the [completed distance-nine proof](fresh_selective_distance9_dimension108_20261017.md). The new ports, joins, identifications, corner domains and extension conclusions are recomputed. There is no assumption of arithmeticity or simple ideal vertices.

## Domain and the unbounded order-three family

A connected core of order k containing two black roots has marked radius at most k−2. The complete high-order catalogue has maximum radius seven. If both cores have order at most five, every possible cross edge gives a path of length at most seven, and every overlap gives a still shorter path. With no cross edge their two indefinite blocks already give two negative directions. Thus only a small/high pairing needs attention.

A minimal order-two core has only its two black nodes. Every contact with a high core gives a path of length at most eight; otherwise the indefinite blocks are orthogonal. It cannot occur in an unpaid configuration.

An order-three core has a single white w, at distance one from its black pair. An overlap is impossible, because1+7=8. Any disjoint cross edge must meet a high port at distance seven. The strict actual-edge inverse test leaves precisely L25 node9, with high black pairs {1,2}, {1,3}, {2,3}. Its full inverse diagonal is exactly

    (G_L25 inverse)99 = 0.

Write the small Gram as

    [[2,-x,-a],[-x,2,-b],[-a,-b,2]].

The black-black magnitude x can have an unbounded finite label or be parallel. A truly minimal order-three core has x≤2: a dotted black-black pair would already be an indefinite proper order-two subdiagram. Its two black deletions are positive definite. If one also applies the low port inverse test, a²,b²<3, so the white-incident labels are right,3,4 or5. No finite bound on the black-black label follows or is assumed.

There is nevertheless a uniform exclusion. Every cross contact lands at high node9. Eliminating the nonsingular Lorentzian high core leaves the small Gram unchanged, since the Schur correction is identically zero. The high core supplies one negative direction and the small indefinite core supplies another. This holds for **every real x,a,b and cross magnitude for which the small core is indefinite**, including the physical unbounded-label family. It also covers zero cross contact. The proof is a symbolic identity, not a sampled angle limit.

The same zero-correction identity permits an arbitrary column of low contacts when all high endpoints are node9; the finite computations below retain every such contact pattern rather than rely on a single-edge approximation.

## New finite order-four/five domain

The frozen low-core census proved, before enumeration, that order-four finite labels are at most six and order-five labels at most five. It enumerated all weighted graphs, not just trees, and found32 and14 canonical minimal-indefinite types. The current checker recomputes every black pair with positive-definite deletions and every white port of distance at least one.

For disjoint joins, low radius one requires high radius seven; radius two requires high radius at least six; radius three requires high radius at least five. The high inverse filter gives18 relevant ports, belonging to L18,L22,L25. Every such port has maximum normalized black-deletion inverse diagonal at least6/7. Actual-edge positive definiteness therefore forces any cross magnitude beta to satisfy beta²<7/3. Since the label-five squared magnitude is (3+sqrt(5))/2>7/3, only labels3 and4 remain. This includes parallel and dotted exclusions and is not an imposed label cutoff.

The checker applies both endpoint inverse tests and the supporting-white Schur condition

    beta² (G_Wlow inverse)xx (G_Whigh inverse)yy < 1.

It then eliminates the full high core and computes the exact inertia of the remaining small Schur matrix. The726 low-first joining-label cases split as follows:

| Low order | Simple two-negative joins | Simple white-PD rejections | Label-four two-negative joins | Label-four white-PD rejections |
|---|---:|---:|---:|---:|
|4|294|11|17|6|
|5|342|22|15|19|

Thus668 joins have two negative directions and58 violate white positive definiteness. The131 underlying small Schur types are recorded exactly. Unlike the distance-nine case, some doubled joins pass the initial white tests; their full two-negative inertia is essential.

There are417 marked context pairs with at least one surviving edge choice:249 have one possible cross-edge position,117 have two and51 have three. Every one of the321 subsets containing at least two positions already violates actual-edge positive definiteness when all selected bonds are simple. Explicit nonnegative vectors of nonpositive square certify these rejections. Larger native magnitudes only decrease the squares, so all mixtures of simple and doubled edges are covered.

## Complete overlaps, including two shared nodes

An overlap needs marked-distance sum at least nine. Order four therefore contributes radius two against high radius seven. Order five contributes radius two against seven, or radius three against six or seven. All raw eligible high ports are retained before the inverse test: with two shared nodes an endpoint need not have an exclusive new neighbor, so filtering every shared high node in advance by its port inverse would be invalid.

The checker tests all injective partial identifications and every contact on their shared roots. Among630 raw candidates,523 are consistent. Conditional inverse tests reject422; actual-edge tests reject35. The66 remaining low-first contexts have14 literal core types, all involving L25. They comprise12 order-four and54 order-five contexts. Reversing the cores covers the other orientation.

All four adjacent corner Grams are positive semidefinite in every survivor. An additional bond between exclusive parts of a core either creates a path of length at most eight (2,184 tests) or has an actual-edge obstruction (12 tests). Those nonnegative witness certificates remain valid after any extra acute contacts. Thus each full surviving core Gram is fixed before its inverse or radical is used.

## Complete supporting-white extensions

At a selected corner, a further supporting white is right to every affine component by its positive kernel. On a finite component a possible contact must have normalized inverse diagonal at most two. The omitted black roots are retained in the possible-contact set. These tests cover contacts anywhere on the core.

The66 cases divide into three complete mechanisms.

* **39 radical cases.** The core has one negative direction and a one-dimensional radical. In an actual Lorentzian realization the radical is a linear root relation: its null root combination lies in the positive-definite perpendicular space to a timelike vector and hence is zero. The selected corner contact domain avoids every zero radical coordinate. Acuteness and the relation force any nonzero contact vector to meet both signs, giving a path of length at most eight. Every extra supporting white is therefore orthogonal to the core.
* **24 cases with no relevant contact.** The nonsingular core is Lorentzian. The contact domain splits into two positive inverse end blocks. A root meeting both ends gives a short path; in either one end its projection cost is at least any one contacted inverse diagonal. Every possible contact in the marked component has full-core inverse diagonal greater than two. Thus no outside root can touch that component. Projection-one contacts elsewhere may allow arbitrarily long outside chains; no bound on those irrelevant chains is claimed or needed. There is no route into the marked component without a first contact to it.
* **Three sealed simple extensions.** The only potentially relevant core-contact node has inverse diagonal p=(5+sqrt(5))/4. All other contacts are excluded by the same end-domain tests. A doubled bond already costs2p>2, so the one allowed contact is simple. Its positive-space residual has square rho=(3−sqrt(5))/4<1/2. Two roots with that same contact would have actual product at least2p−2>0, so at most one occurs. Every remaining core-orthogonal norm-two white has product magnitude at most sqrt(2rho)<1 with it. The normalized acute Coxeter gap forces that product to be zero. The new root therefore has no indirect tail or bridge, and adds at most one node.

The source retains every corner decomposition and then chooses a certified corner. The first two mechanisms give component order at most12; the third gives at most12 after its one possible new root. Arbitrarily many additional supporting whites are allowed in these arguments. The conclusion is not a bounded one-root search.

At an ideal corner, the inherited actual factor-fiber lemma places the whole varying affine factor inside W plus the marked pair. Other incident cusp walls cannot enlarge the component outside this calculation. Hence the component bound has the same meaning for finite and ideal corners, without simple-cusp assumptions.

## Exact certificates and scope

Finite/affine positive-semidefinite recognition uses the frozen exact tree/cycle criterion. Every failed strict actual-edge test now includes an explicit nonnegative nonpositive-square vector. A shortest cycle uses an all-ones vector; for a forest, positive leaf Schur eliminations construct the witness by positive coefficients until the first nonpositive pivot. The source verifies each stored vector in the native normalized Gram. Squared tree bonds lie in Q(sqrt(5)); the remaining small matrices and inverse/radical calculations use exact SymPy arithmetic. No floating-point test is used.

Generation15769 and fresh saved replay59083 both terminated with exit0 and `PASS_COMPLETE_DISTANCE8_SMALL_CORE_DOMAIN`. The replay recomputed the exact certificate and verified its input bindings. The parent agent independently read the complete new source and proof, checking the radius domain, analytic triangle Schur argument, doubled-join and multiple-edge coverage, conditional overlap inverse tests, nonnegative obstruction certificates, all three extension mechanisms and ideal-factor transfer; this was a narrow mathematical/source review, not a second independent enumeration.

The result completes the small-core portion of a distance-eight argument only. The remaining high-core disjoint/overlap and extension geometry, and any resulting dimension theorem, are separate obligations. No external refereeing, formalization, construction or publication-priority claim is made.
