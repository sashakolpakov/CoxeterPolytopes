# A finite local target for the compact dimension-29 problem

Any hypothetical compact hyperbolic Coxeter polytope in dimension 29 must have a vertex in an explicit list of **2,126 spherical Coxeter forests**, with only nine possible component-rank partitions. This is an unconditional necessary condition: no arithmetic, facet-count, or nerve assumption is imposed. It is **not an exclusion of dimension 29** and is not a classification of full polytopes.

The [exact checker](../scripts/fresh_compact_dimension_bound_20261015.py) and [saved certificate](fresh_compact_dimension_bound_20261015.json) enumerate the complete local list. They check graph weights, integer partitions and the combinatorial identity below. A [separate independent review](fresh_compact_dimension_bound_20261015_review.md) reproduces the h-vector arithmetic and9/2,126/1,354 counts, and checks the occurrence argument. It inherits the connected graph weights after source review instead of repeating the66-graph census. The spherical-diagram classification and the global planar-face inequality remain published theorem inputs.

The unrestricted compact upper bound remains 29 in the checked primary survey, [Lee–Marquis, Theorem 3.10 and Remark 3.12](https://arxiv.org/pdf/2109.06758). Alexandrov's stronger bound 12 concerns the additional hypothesis that every missing face has size 2; it does not exclude general compact dimension 29. [Alexandrov, Theorem C and §4](https://arxiv.org/pdf/2203.07248).

## The exact combinatorial threshold is 84

Use Vinberg's planar-angle weights. For a vertex v write σ(v) for its total weight, and for a k-gonal two-face F write σ(F) for the sum at its corners. The published construction has σ(v)≤3n and σ(F)≥5−k. Its componentwise weights and face proof are in the original [Vinberg paper, §§1 and 6, printed pp.70–72 and 80–91](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

For n=29, every vertex lies on 406 two-faces. The exact odd-dimensional average-side bound is 116/28=29/7. Consequently f₂>98f₀, and

\[
\sum_v\sigma(v)\ge5f_2-406f_0>84f_0.
\]

Thus a universally valid weighting with vertex bound 84 would suffice. This is a less demanding numerical target than substituting c=23/8 into the coarser formula n<8c+6. Existing Vinberg weights give 87, so the outstanding gap is three units per vertex in a uniform bound, or an appropriate global averaged inequality.

There is an exact version of this comparison. Let h₀,…,h₂₉ be the h-vector of the dual simplicial polytope; hᵢ=h₂₉₋ᵢ, hᵢ≥0 and h₀=1. The standard face-count identities give

\[
 f_2-98f_0=\sum_{i=0}^{14}(14-i)(15-i)h_i.
\]

Define s(F)=σ(F)−(5−k_F)≥0. Double-counting planar angles gives the exact balance

\[
\boxed{\quad
\sum_v(\sigma(v)-84)
=5\sum_{i=0}^{14}(14-i)(15-i)h_i+\sum_Fs(F).
\quad}
\]

Its right side is at least 1050 using h₀ alone. Because σ(v)≤87, at least 350 vertices must have σ(v)>84. This modest lower count uses no additional h-vector theorem. The certificate checks all 15 coefficients of the identity.

There is a necessary occurrence refinement. Set every angle weight on a face with five or more sides to zero, retaining the original weights on triangles and quadrilaterals. Denote the new vertex total by σₛ(v). All face inequalities remain valid: on a larger face the new total is0≥5−k. Also σₛ(v)≤σ(v)≤87. The same balance therefore forces at least 350 vertices with σₛ(v)>84. **At least one listed extremal corner has more than 84 units of weight on its incident triangles and quadrilaterals.** This supplies the small-face occurrence premise; extremal vertex existence alone would not supply it.

## Complete enumeration of the extremal vertex types

Let w(r) be the contribution of an irreducible spherical component of rank r, and δ(r)=3r−w(r). The result of applying the published weighting is independent of the irreducible type at fixed rank:

| Rank | w(r) | δ(r) |
|---|---:|---:|
|1|0|3|
|2|2|4|
|3|6|3|
|4|12|0|
|5|10|5|
|6|15|3|
|7|21|0|
|8|24|0|
|9|27|0|
|r≥10|r+19|2r−19|

For the classical A/B/D graphs the checker constructs every relevant end segment of order at most 6 and counts the union of its node pairs with adjacent pairs. It directly checks 66 graphs of ranks 8 through 29. At ranks 8 through 12 it also checks the end-segment list against a separate exhaustive connected-subset definition. For D, invariance under swapping the two fork leaves is enforced. The E₈ contribution 28−4=24 is an explicit published input; the checker does not independently reconstruct the four pictured exceptional pairs or repeat Vinberg's face proof.

Since a rank 29 forest has σ=87−Σδ, the strict inequality σ>84 requires Σδ≤2. All components therefore have ranks 4,7,8,9 or 10; there can be at most two rank 10 components. Exhausting all 4,565 partitions of 29 leaves exactly the following list.

| Component ranks | σ | Typed forests | Without H₄ |
|---|---:|---:|---:|
|4+4+4+4+4+9|87|378|168|
|4+4+4+7+10|86|420|240|
|4+4+4+8+9|87|420|240|
|4+4+7+7+7|87|300|200|
|4+7+8+10|86|240|192|
|4+7+9+9|87|120|96|
|4+8+8+9|87|150|120|
|7+7+7+8|87|80|80|
|9+10+10|85|18|18|
|**Total**||**2,126**|**1,354**|

The available irreducible types are A₄/B₄/D₄/F₄/H₄; A₇/B₇/D₇/E₇; A₈/B₈/D₈/E₈; and A/B/D in ranks 9 and 10. Repeated equal-rank components are unordered. The certificate saves every typed forest, rather than only these counts. B and C are the same unoriented Coxeter type here.

No angle-label cutoff was used. Every finite I₂(m), including arbitrarily large m, has deficit 4 and therefore cannot occur in an extremal forest. In fact all extremal vertices use only labels 3,4,5. Of the 2,126 types, 385 contain E₇, 245 contain E₈ and 772 contain H₄.

## A known obstruction to the most obvious optimization

Vinberg explicitly warns in §6.12, printed p.91, that a uniform improvement c<3 is obstructed when nonnegative weights depend only on the planar-angle diagram and vanish across its components. E₇ supplies the obstruction through quadrilateral stars. This is a published limitation of that local architecture, not a new theorem proved by this checker. It must not be ignored when setting up a reweighting program. [Original remark](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

A rank 29-specific weighting depending on the full forest, weights across components, or a global charge argument requires separate analysis; the remark alone is not being promoted into an impossibility theorem for all such strategies. The exact balance above provides a test: any proposed global charge must control its positive vertex excess while retaining the nonnegative face slacks.

Our [Lorentz-inertia ledger](fresh_h9_lorentz_inertia_20261015.md) supplies valid forbidden configurations in arbitrary Lorentz dimension only under their literal spherical-face and contact premises. The strongest high-label configurations are not present inside the extremal vertices, whose components have no label ≥7. They cannot simply be declared unavoidable by increasing dimension.

## Concrete next verification target

The 2,126 forests are finite entry points for an exhaustive local analysis. They are not by themselves a finite algorithm deciding every chamber completion. For each forest T and each pair of its nodes, delete those nodes to form a rank 27 white diagram W. A triangular or quadrilateral two-face through that angle introduces respectively one or two further black nodes.

The occurrence refinement gives a particularly small certificate target. If the original vertex weight is W∈{85,86,87}, a vertex with σₛ>84 can lose at most W−85 units to larger faces. Thus it suffices, for **each** listed type, to prove that some marked pairs of total original weight at least W−84 cannot occur on triangles or quadrilaterals. The required blocked weight is respectively 1, 2 or3. A pair of weight2 counts twice. This is an exhaustive marked-angle coverage requirement, rather than a claim that one arbitrary failed extension eliminates a corner type.

The necessary conditions can be stated exactly before any numerical search:

1. Each adjacent pair of black nodes, together with W, gives a positive-definite rank 29 vertex diagram. Each single black node together with W is also positive definite.
2. After taking the Schur complement of W, a triangular star has a 3×3 residual of Lorentz signature (2,1). A quadrilateral star has a 4×4 residual of rank 3 and signature (2,1,0); the opposite black pairs give indefinite residual pairs. All full-star principal subdiagrams must have negative index at most one, and compactness excludes parabolic subdiagrams.
3. Every discrete label and continuous dotted distance must remain covered. New rank-two spherical factors can arise after deletion, even though the original extremal forest has none; such I₂(m) tails cannot be truncated without proof.

These conditions give the appropriate interface for the universal inertia obstructions. Rejecting all small-face extensions for every type would suffice; more economically, the blocked-weight coverage just described would exclude compact dimension 29. Alternatively, an exact charge or modified-weight certificate satisfying every retained triangle/quadrilateral constraint could contradict the balance identity. **Neither rejection nor that global certificate has yet been produced.** Merely verifying a local corner Gram, or rejecting a bounded selection of extension labels, would not achieve the goal.

## Separate arithmetic refinement

For compact arithmetic groups in dimension at least 22, the relevant ground-field candidates are Q(√2), Q(√5), and Q(cos(2π/7)); this is a separate hypothesis from the unconditional reduction. [Vinberg, Theorem 2 and §7](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

An H₄ component contains a label 5 edge. Its length-two cyclic product gives 4cos²(π/5)=(3+√5)/2, so its ground field must contain √5. Q(√2) does not contain √5, and the degree-three field Q(cos(2π/7)) has no quadratic subfield. Therefore those two arithmetic cases retain only the 1,354 H₄-free extremal forests. Q(√5) retains all 2,126 at this stage. Labels 4 in B/F components have rational squared cosine and are not removed by this argument. No arithmetic lattice or chamber classification is claimed.

Generation and saved replay pass:

```sh
python scripts/fresh_compact_dimension_bound_20261015.py --verify
```

The compact-star search remains frozen at its 60,015-visit cursor. This package does not change any prior construction, obstruction or search certificate.
