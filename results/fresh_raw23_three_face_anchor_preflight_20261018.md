# Raw23 three-dimensional faces: exact native-contact preflight

No norm-two normal satisfying Coxeter contact conditions with the known walls can have its Lorentz projection in the relative interior of any of the552 entirely spacelike three-dimensional faces found in the three original raw23 cones. The proposed normal need not have ground field Q(√5): the exclusion applies also to larger-field or nonarithmetic extensions of these fixed known cones. Arbitrary components orthogonal to the known Lorentz span are allowed. This does not exclude arbitrary acute completions or any raw23 pair, and it adds no capacity, model column or arithmetic bound.

The new [producer](../scripts/fresh_raw23_three_face_anchor_preflight_20261018.py), [certificate](fresh_raw23_three_face_anchor_preflight_20261018.json) and [saved checker](../scripts/check_raw23_three_face_saved_20261018.py) bind the frozen [extreme-ray polar packet](fresh_raw23_union_polar_completion_preflight_20261018.md). The separate frozen [two-ray packet](fresh_raw23_two_ray_completion_preflight_20261018.md) covered336 adjacent spacelike two-dimensional faces. This packet addresses a larger non-extreme domain without a native125-choice search.

## Complete stated face domain

Each of the three original cones, of kinds(1,1), (1,6) and(6,6), is pointed and full-dimensional in its ten-dimensional known-wall Lorentz span. The frozen polar certificate gives all50 projective extreme rays:32 spacelike and18 timelike. Completeness and the exact Q(√5) geometry of these known inputs remain inherited; no field restriction on an additional normal is used below.

For each cone the producer examines all4,960 triples of distinct spacelike rays. It retains a triple exactly when the rays are linearly independent and the known-wall rows common to their zero sets have rank7. Because every ray has nonpositive known-wall products, a positive sum of the three has exactly those common active walls. Their rank7 cuts out a three-dimensional face, and the three independent rays already span that face. Conversely, every three-dimensional face all of whose extreme rays are spacelike has an independent triple in this domain.

For each retained triple the producer gathers every one of the50 original extreme rays satisfying its active equations, intersects their zero sets again, and checks that their span has dimension3. Deduplication uses the complete active-wall set. Thus an extra fourth ray or a timelike ray cannot be silently discarded by selecting one convenient triple.

|Per original cone|Count|
|---|---:|
|Spacelike triples inspected|4,960|
|Independent triples defining a3-face|544|
|Distinct faces|184|
|Simplicial faces with3 rays|64|
|Nonsimplicial faces with4 rays|120|
|Faces with a positive-definite3-plane|184|

All552 distinct faces found across the three cones are entirely spacelike and have positive-definite span. There is no non-PD or timelike-containing exception within the found domain. This does not claim enumeration of three-dimensional faces that have no independent triple of spacelike rays.

## Exact maximum on each positive triple cone

Let x1,x2,x3 be an independent ray triple, H its positive-definite Gram matrix, and let p_i=−(x_i·r)≥0 for a known scaled root r of norm2s. For a projected vector x=Σu_i x_i with u_i≥0 and x²≤2, the squared physical contact with that known norm-two wall is

    (p·u)² / s.

If p is nonzero, its maximum on the nonnegative ellipsoid occurs on a nonempty coefficient support J. At an interior stationary point of that support, u is a positive multiple of H_J^−1 p_J. The support is therefore admissible precisely when every entry of this vector is strictly positive. Its candidate maximum is

    2 p_Jᵀ H_J^−1 p_J / s.

The producer considers all seven nonempty supports. Boundary maximizers occur on smaller supports and are consequently included. Taking the largest admissible value gives the sharp maximum over the whole triple cone. A zero contact form has maximum0. The saved checker verifies every stationary vector through H_J z=p_J and independently checks all support/maximizer decisions.

For the192 simplicial faces, the packet additionally verifies three independent known-wall contact forms whose relative-interior squared bounds are below4. In this auxiliary finite-anchor check, a maximum4 is allowed only if it is attained on the actual simplicial boundary. Under the earlier Q(√5) candidate-field restriction this would permit a finite125-choice anchor reduction. That enumeration and field restriction are unnecessary here: every such face already has a stricter obstruction below1.

## Four-ray faces: complete conic cover, not one chosen triple

Each nonsimplicial face has four extreme rays in a positive-definite three-plane. Conic Carathéodory, or directly eliminating a linear dependence from a positive representation, expresses every vector in this cone with at most three independent generators. A representation with fewer generators can be extended to an independent triple with zero coefficients. The producer enumerates every independent triple of the full four-ray set; all four triples are independent in these cases.

It computes the seven-support sharp bound on each triple cone and takes the maximum over that complete cover, separately for every nonzero known-wall contact form. The resulting value bounds the contact on the entire nonsimplicial face. The certificate records all1,632 independent subcones across the552 faces, rather than assuming that a single simplicial piece controls a four-ray face.

No strictness is inferred from an artificial triangulation boundary. Such a boundary can pass through the relative interior of the actual four-ray face. The exclusion below uses only a strict maximum less than1 over every subcone, so this issue does not arise.

## Why the strict bound excludes arbitrary passive lifts

For every one of the552 faces, at least one known-wall contact form is nonnegative and nonzero on the full ray set, while its sharp full-cone squared maximum is strictly less than1. A nonzero nonnegative linear functional is strictly positive on the relative interior of the face: otherwise its zero set would be a proper face containing an interior point. Hence a projected normal in that relative interior has a known-wall squared contact strictly between0 and1.

For any real Coxeter angle between distinct norm-two wall normals, the squared contact is0 at a right angle, while a nonright finite angle π/m has

    4*cos²(π/m) >= 1,    m >= 3.

Parallel and ultraparallel walls have squared contacts at least4. A contact in(0,1) is therefore impossible for any Coxeter extension, independently of a ground field or arithmeticity. The bound above ranges over all real nonnegative ray coefficients, so it has not silently restricted candidate normals to Q(√5). That field is used only to compute and verify the exact known-input coefficients. No conjugate or additional unknown-wall gate is needed for this rejection.

Write a full proposed norm-two normal as n=x+z, where x is its projection to the known Lorentz span. That span contains a timelike vector, so its orthogonal complement is positive definite. Thus x²+z²=2 implies x²≤2. Every known-wall product equals the product with x. The maximization already allowed every projected vector in the face with squared norm at most2; allowing arbitrary passive z cannot evade the bound. No product with an unknown passive wall is assumed zero.

The native125-choice search count is zero. There is no viable additional normal in these relative interiors, but this is a restricted domain result. Boundaries, timelike-generator faces, higher-dimensional faces, projections outside these faces, and faces arising after a previously added pole require their own arguments. The frozen two-ray and extreme-direction packets cover some boundary domains; they are not silently extended here. An arbitrary completion, pair exclusion, shared vertex capacity or arithmetic improvement remains unproved.

## Verification scope

The initial classification-only generation53764/2c89ae and separate replay16791/b1247f passed before the complete four-ray extension. They established the face domain but left nonsimplicial faces unresolved. Extended generation47557/fc1a68 and separate producer replay40841/06d9cc terminated exit0/PASS with all552 native exclusions. The first independent saved check43904/7b5e71 and final saved check76247/684ee4, after explicit headline/no-claim guards, also terminated exit0/PASS. Metadata preflight9d4375 checked5 direct bindings,2 Python ASTs and5 local links.

Construction completed a full mathematical producer/report/standalone-helper read (fda88c/8057e8), with a saved-headline comparison(d1d941): PASS. It checked the complete stated face domain, full-ray expansion, all independent-triple covers, sharp seven-support maxima, strict interior contact argument and arbitrary passive-norm relaxation. It performed no execution and did not reprove the inherited cone, native field or catalogue. The requested explicit552-face and no-claim guards were added to the helper before its final run; the final guard delta read(acf599) also passed.

The parent completed a full mathematical producer/report/standalone-helper read (a07b96/332de9/b9084e): PASS, covering all14,880 triples, full-face expansion, positive-definite spans, conic Carathéodory covers, sharp seven-support maxima, nonzero interior functionals and passive projected-norm relaxation. It performed no execution and did not replay the inherited ray catalogue or field proof. That review identified the stronger field-independent conclusion: the universal Coxeter minimum1 suffices, so the report now imposes no candidate-field restriction. This is a proof-scope clarification; the source, data and checker calculations are unchanged. All computations and both stated reviews are complete, and the packet is final and frozen.

The producer uses exact Q(√5) arithmetic and freshly reconstructs all14,880 triple decisions, full ray memberships, plane Grams and complete support-wise bounds. The standalone checker imports no producer; it uses the earlier saved-checker arithmetic design, reverse-column rank elimination, direct Gram products and linear-equation checks of every stationary certificate. These are independent checking operations, not a claim of unrelated field-arithmetic designs. Both inherit the frozen complete ray lists, actual-star interpretation, native catalogue and field description of the known cones. Neither performs native125 enumeration, a full chamber computation, optimization or a model replay.
