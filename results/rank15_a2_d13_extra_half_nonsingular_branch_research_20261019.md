# The four D15 routes do not need the extra A2/D13 half

19 October 2026. Written application of the reviewed nonsingular full-seal
rank lemma, with saved-record extraction only. No checker, matrix calculation,
search, forest program, Lean or build was executed. All predecessor bytes are
unchanged. This establishes a local replacement for four historical routes;
it does **not** yet authorize deleting the extra half from the complete
current rule.

## Exact route and core

The rank15 sparse exchange sets the two whole-D15 pairs(1,13),(1,14)
from half to zero. Its four whole-component contexts are1162,1175,1910,2303.
At a distinct corner it upgrades the whole A2 endpoint / D13 short-leaf
payment from half to unit. There are four physical marked pairs per
A2/D13 component pair, so this upgrade costs2 per component pair.

The saved B15/D15 theorem explicitly supplies a **nonsingular18-root**
Lorentz core of signature(17,1) in every whole case, not the17-root core
used by its proper D14 cases. Its core key is
`(25, 25, ((8, 9), (9, 8)))`. There are14 displayed common white roots
and four black side roots. Its selected finite corner has whole D15
plus a displayed A1, hence16 incident roots.

| Context | D15 selected marks | A2/D13 coupon marks | Extra A1 at selected corner |
|---:|---|---|---|
|1162|(0,11)|(2,11)|3|
|1175|(0,12)|(2,12)|3|
|1910|(1,10)|(1,12)|13|
|2303|(2,10)|(2,12)|13|

In1162/1175 the remaining allowed contact positions after imposing the
whole D15 condition are2,3,13; in1910/2303 they are3,12,13. Their inverse
diagonals are4,4,14. The inherited complete arbitrary-contact theorem
partitions these positions into the two end sets. A contact reaching both
ends gives an independently paid short8 path. Every nonzero single-end
contact has projection square at least4, larger than its norm2. This
argument allows all permitted contacts and arbitrary native magnitudes;
it is not an assumption that omitted spectator contacts vanish.
Consequently every further actual common white is orthogonal to all18
core roots in the no-short8 branch. This is exactly the whole-core seal
proved in the frozen B15/D15 report and checked by its source `proof()`.

## Actual-face rank contradiction

Let the ambient dimension be n and assume an actual quadrilateral realizes
one of those whole-D15 routes without the independent short8 escape.
At a point in its relative interior the common-white normals span a
positive subspace of dimension n−2. They are independent: if a dependence
is split into disjoint positive sums P=M, acuteness gives
P²=P·M≤0; positivity forces P=M=0, while strict pairing with an interior
point of the polytope rules out a nonzero positive sum equal to zero.
There are therefore exactly n−2 common-white normals.

Fourteen are in the displayed core, leaving n−16 further common whites.
The full seal makes their span E orthogonal to the core span H. Since
the18-root Gram matrix is nonsingular, dimH=18 and H is nondegenerate,
so E intersects H trivially. Independence gives dimE=n−16. It follows
that the ambient Lorentz space would have dimension at least

    dimH+dimE=18+(n−16)=n+2,

contradicting its actual dimension n+1. If n<16, the displayed16-root
finite corner already cannot occur. Thus these four no-short8 branches
are impossible in every ambient dimension.

In the remaining branch the independent short8 route already pays the
face. The historical half+half or zero+unit pairing is unnecessary for
these four actual faces. This application uses the rank lemma only on a
nonsingular core; it does not exclude a19-root corank-one core or assert
an invalid capacity restriction for a triangular prism.

## Payment consequence and remaining interface gap

Locally the D15 zeros may therefore be retained while restoring the
A2-end/D13-short coupon from unit to its older half. Its old fork-center
half is retained as well. Relative to the current component charge, the
candidate decreases the finite cost by2*a2*d13 and changes no ideal weight.
No new positive family would result from that monotone decrease.

This is not yet a complete current payment theorem. The earlier broad
small-component/D13 theorem already needs the retained half. Every later
physical output asking for **more** than that half must still be audited;
if a later route uses the upgraded unit without one of these four complete
whole-D15 models, the above local contradiction does not cover it.
Subcore inclusion or a source hash alone cannot establish the required
whole-factor guard. The current D16/D13 worker targets A4-end/D13-short
marginal alpha and does not perform this separate A2/D13 audit.

For scale only, the old formal forest3D13+4A4+A2 would save6, but the
current physical A4/D13 capacity proposal already lowers its allowance
below364. The unchanged370-family forests with no A2/D13 pair are not
affected. No new global finite maximum or dimension bound follows here.

## Pinned inherited inputs

* `results/finite_quad_nonsingular_full_seal_rank_obstruction_20261019.md`,
  SHA256 `9d16f3bcdcefe85d3138aae6967a484907f3978badf02eef6f5fde3cc878a2c4`.
* `scripts/fresh_finite_b15_d15_existing_payment_transfer_20261018.py`,
  SHA256 `b562201a79642540bd208c716860e0294f1bf0707fbf7daee4547edcac5199a1`;
  `proof()` asserts signature(n−1,1), whole-case n=18 and the inverse
  diagonals4,4,14 before recording all four corner roles.
* Its saved JSON, SHA256
  `88490ce0b2401edf37f4ec45321bfa3e4fd50afac958345f76c4006567e8e290`;
  its report, SHA256
  `ddabd354250b2b7915e3662549136966e2093224c1e5fc22b3d30c5f49c40205`.
* `results/fresh_rank15_sparse_payment_exchange_20261018.json`, SHA256
  `50d7949f65aa549a90048b5cd1badf57572a7065746ff3f93770d22aec4d2770`;
  its report, SHA256
  `ae8f89387c6e34ee9a2dec54c6cb47e0c78625d4c340a61d07e20e33d703cc02`.
  These identify the precise half-to-unit upgrade being reconsidered.

The complete native contact domain, independent short8 routing and exact
inertias are inherited written/computational premises. This note does not
repeat their historical enumeration or claim kernel verification of their
geometric completeness.
