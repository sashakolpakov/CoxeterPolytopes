# Hostile scope audit of the arithmetic dimension-23 candidate

**No mathematical gap was found in the transitions audited here.** This is a fresh mathematical and source audit of the frozen [dimension-23 checkpoint](arithmetic_dimension23_and_complete_face_coverage_2026_10_15.md), with a separate addendum for the new dimension-23 weight rule and its bound-22 composition. It is not another full catalogue or queue computation, and it does not establish external validation or priority.

The [standalone checker](../scripts/fresh_arithmetic_dimension23_scope_audit_20261016.py) and [certificate](fresh_arithmetic_dimension23_scope_audit_20261016.json) use only the Python standard library and exact `Fraction` arithmetic. Generation and a fresh saved replay pass. They independently check the h-vector scalars, finite-type and hub-context counts, all **112,764** saved four-corner payment sums in the two binders, completeness of the spectator-count quantifiers on every row, and the two field/dimension Cartesian products. They do **not** reconstruct the saved corner weights or the forest maxima from diagrams; those are supplied by the named independent weight reviews.

## 1. Where compactness is needed

Use the Lorentz convention with spacelike outward normals, nonpositive distinct-wall products, and a strict future timelike interior point H.

If an affine principal diagram existed, take a positive null-kernel vector on one irreducible affine component and form its corresponding positive linear combination delta of wall normals. Then delta has square zero. It is nonzero because delta.H is strictly negative. Acuteness gives delta.n <= 0 for every remaining wall, so delta represents an ideal point of the chamber. This contradicts compactness. The argument concerns the **whole compact chamber**, not merely a compact two-face in a noncompact chamber.

A positive-definite principal Z-Gram B does give an actual face. Let b = B^{-1}(H.n_i). Since B^{-1} is entrywise nonnegative, b <= 0. The projection H_0 = H - sum b_i n_i is timelike, is orthogonal to the chosen walls, and satisfies every other wall strictly: subtracting sum b_i(n_i.n_j), a nonnegative quantity, only decreases H.n_j. Thus a purported opposite-side subset of a quadrilateral cannot have positive-definite Gram. It cannot be positive-semidefinite singular either, by the preceding affine argument. It is indefinite.

Minimal indefinite subsets in this compact setting therefore have all proper subsets positive definite: they are Lanner diagrams. For a quadrilateral with white set W and opposite black nodes a,c, a minimal indefinite subset of W+{a,c} contains both a and c, because W plus either one is contained in an actual finite corner. Similarly, a minimal indefinite subset of a triangular star contains all three black nodes. This justifies the marked-piece inputs; it does not assume that a minimal union is already the complete face star.

At a compact Coxeter vertex the incident Gram is positive definite and the chamber is simple. A quadrilateral two-face in dimension n has n-2 white and four black walls, hence n+2 nodes. Its rank is n+1: after quotienting the positive white span, its four side normals span the Lorentz three-space of the compact polygon. A smaller span would leave an unbounded direction or a common pencil, incompatible with compactness. Consequently a nonsingular 31-node principal star cannot occur when n <= 29. The cutoff is a dimension bound, not an empirical rank cutoff.

These facts also agree with the star conditions reproduced in [Burcroff, Proposition 9.3](https://arxiv.org/pdf/2201.03437). They do not use the invalid inference that every compact face of a finite-volume chamber has no affine subdiagram in its ambient star.

## 2. Arithmetic restrictions and dotted parameters

The cyclic-product normalization is consistent. For normals of squares d_i,d_j and product g_ij, the length-two cyclic invariant is

    t = 4 g_ij^2 / (d_i d_j).

For an arithmetic reflection chamber it is an algebraic integer in the ground field. The other real embedding has positive-semidefinite principal Grams, giving 0 <= t' <= 4. This is the normalization obtained by putting every root at square 2; no assumption that every unscaled Gram entry is integral is made. [Nikulin, Section 2, formulas (1)–(3)](https://arxiv.org/pdf/0708.3991) states the cyclic-integrality criterion and also explains that an indefinite principal subdiagram has the same arithmetic ground field. In particular, no maximality assumption is needed in this reduction.

For each of the 24 surviving minimal dotted cores, the reviewed common submatrix B obtained by deleting the two dotted endpoints is nonsingular Lorentzian. In any larger actual chamber the endpoint Schur complement must be positive semidefinite: B already accounts for the unique negative direction. With residual diagonal a,b and off-diagonal -x-c,

    (x+c)^2 <= ab,  hence x^2 <= 2ab+2c^2.

This gives the recorded upper bound t <= U, independently of all additional white nodes. In Q(sqrt5), write t=A+B phi. Since t>4 and 0<=t'<=4,

    B=(t-t')/sqrt5 > 0,       A=t'-B phi' > 0.

Also A<=U and B<=U, so the integer search box is complete. The endpoints t'=0 and t'=4 cannot occur for a nonrational t>4, by injectivity of the conjugation map. Thus the 532 checked lattice points and 32 surviving parameter seeds do not represent a distance sample or finite-label cutoff. Their arithmetic remains inherited from the [independent dotted-parameter review](fresh_sqrt5_dotted_parameters_20261015_review.md).

For a Lanner determinant, the ground-field assertion can also be seen directly. The normalized determinant lies in the cyclic field, is negative in the physical embedding, and is positive under every nonidentity ground-field embedding. If it generated a proper subfield, another embedding fixing that subfield would retain its negative value, a contradiction. This is the exact role of [Vinberg, Proposition 17, printed page 91](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt); rational equivalence of unrelated integral lattices is not being used.

## 3. Why the attachment queues cover arbitrary white extensions

The initial no-short-path reduction is legitimate only because the retained baseline already pays any short corner; it is not a claim that all quadrilateral stars have no short corner. The independently regenerated 5,440 literal minimal unions and 1,180 marked classes cover the unpaid case, with all cross contacts and shared-white identifications included. See the [incidence review](fresh_sqrt5_quad_exhaustion_20261015_incidence_review.md).

Given an actual larger unpaid star, order the white nodes outside an initial core by connected growth. Any outside component touching the core has a first directly adjacent node. Its contacts occur in the complete one-white domain. Positive-definite corner tests, field conditions, kernel compatibility, and physical/conjugate inertia failures are necessary principal-submatrix conditions; adding more nodes cannot repair them. An existing short path remains short, and a paid entire H4 corner component cannot acquire another vertex in a finite Coxeter forest. These two baseline pruning rules survive further attachments because all subsequent bonuses are nonnegative.

If a finite core has no unpaid accepted touching profile, every remaining white node is orthogonal to it. There is no path that evades a first contact. The remaining white Gram is positive definite. Hence a sealed nonsingular core plus such spectators stays nonsingular and cannot be a quadrilateral star, whose corank is one. For a singular core, its zero-further-touch case is a genuine payment obligation even if other touching extensions are also possible. The complete queues keep **all** such singular cases, including expanded and boundary states. A traversal that retained only final sealed states would be insufficient; that is not what these certificates do.

The finite queue's 5,172 states and 2,070 obligations, together with the dotted queue's 438 states and 264 obligations, are complete only in the explicit dimension-at-most-29 range. That range is justified before the search by the published compact dimension bound. The [finite queue review](fresh_sqrt5_nonsingular_rank29_fast_20261015_review.md) reconstructs every state and domain independently; the [dotted queue review](fresh_sqrt5_dotted_rank29_exhaustion_20261015_review.md) does the corresponding independent job for all 32 parameter seeds. This audit reads their completeness arguments and binds their results; it does not repeat those expensive traversals.

Old paid-core shortcuts require care when weights change. The complete binders do not merely reuse an old numerical payment: they reevaluate all 2,661 zero-touch rows with the new rule, including the earlier sealed, singular-extension and dotted branches. The 471 literal maps and the source reviews remain explicit dependencies. The frozen [contact review](fresh_sqrt5_quad_contacts_20261015_review.md) correctly qualifies its 112 singular and 44 nonsingular seals as **unpaid** seals, allowing touching extensions already discharged by the unchanged baseline.

## 4. Whole-forest weights and common spectators

At a finite vertex let a count entire A1 components and b count entire I2(5) components. A1-to-root cross weights divide by the actual a; I2(5)-to-root weights divide by the actual 2b. For fixed a,b, each eligible destination root occurs in exactly a or 2b such pairs. The multiplicities cancel these denominators in the vertex total. The second family excludes A1 and I2(5) destinations, so it neither double-counts A1-to-I2(5) pairs nor adds I2(5)-to-I2(5) pairs. A capped b may select a coefficient context; it must not replace actual b in the denominator. Both reviewed implementations make this distinction.

Once exact a,b are fixed, component costs are additive. The available Q(sqrt5) finite types through rank n are A_1,...,A_n; B_2,...,B_n; D_4,...,D_n; E6,E7,E8,F4,H3,H4; and I2(5),I2(6),I2(10). The aliases A2/I2(3) and B2/I2(4) are counted once. Thus the counts are 74,77,92 at ranks 23,24,29. Allowing finite types that never occur in an actual arithmetic chamber only enlarges the upper-bound domain; excluding a possible type would not be safe.

At a zero-touch face core, the external positive-definite white block is the **same** at all four corners. Its exact A1 and I2(5) counts e,f are therefore common. The internal corner counts may differ and are combined with the same e,f separately. For a k-node core in dimension n the remaining rank is n+2-k. After removing e+2f, possible residual rank is exactly zero or any integer at least two; A2 and A3 realize every such total. There is no remaining rank-one component because A1 was counted separately.

The checker independently enumerates these quantifiers, checks no duplicate/missing case, and sums each saved four-corner payment. It finds:

| Excluded dimension | Exact vertex contexts | Zero-touch rows | Rows fitting dimension | Common spectator cases | Maximum vertex weight |
|---:|---:|---:|---:|---:|---:|
| 24 | 157 | 2,661 | 1,897 | 60,735 | 267527/5000 < 54 |
| 23, new addendum | 144 | 2,661 | 1,764 | 52,029 | 485211/10000 < 99/2 |

The diagram-dependent weight values and forest dynamic programs remain supplied by the [rank-24 independent review](fresh_sqrt5_dimension24_i25_cross_20261015_review.md) and the [new independent hub-count review](fresh_sqrt5_hub_count_weights_20261016_review.md). The new shared CLI correctly synchronizes dimension, cap, target, input name and preceding complete coverage. Its bound-22 binder excludes dimension 23 only after checking the complete predecessor that excludes 24 through 29.

## 5. The strict global counting step

For a simple n-polytope, each unordered pair of incident facets at a vertex identifies one two-face angle, so the total number of such angles is binom(n,2) f_0. If nonnegative weights pay at least 2 on every triangle and at least 1 on every quadrilateral, then

    total weight >= 5 f_2 - binom(n,2) f_0.

This remains valid on k-gons with k>=5 because their required value 5-k is nonpositive. No condition on their Coxeter diagrams is needed.

Writing f_0=sum h_i and f_2=sum binom(n-i,2)h_i, Dehn–Sommerville symmetry pairs the coefficients into

    ((n-i)(n-i-1)+i(i-1))/4.

The h_i are nonnegative, h_0=h_n=1, and the extreme coefficient is strictly greater than the minimum. The resulting average lower bound is therefore **strict**, including when the saved vertex maximum equals the threshold. The exact thresholds for n=24,...,29 are 54,60,65,143/2,77,84; for n=23 it is 99/2. This explains why equality at dimension 29 still contradicts existence. This argument uses simple-polytope h-vectors, not an assumption that a facet or two-face is itself a Coxeter chamber. The checker recomputes all coefficients rather than importing the producer's closed formula.

## 6. Other fields and the all-arithmetic scope

For compact dimensions under consideration, Vinberg's original Theorems 1–2 give the bound 29 and, in dimension at least 22, the three possible fields Q(sqrt2), Q(sqrt5), Q(cos(2pi/7)). [Original paper, printed page 69](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt).

The Q(sqrt2) proof needs its complete external-attachment closure, not just 115 union rows or nine finite principal cores. The [reviewed closure](fresh_sqrt2_quad_exhaustion_20261015_review.md) treats the seven sealed singular cores, both nonsingular completions, and dotted connecting-component versus pendant cases. Its weights are separately reviewed. For the cubic field, the [independent review](fresh_compact_global_bound_20261015_review.md) uses both the degree-three obstruction to a nonzero rational quadratic in cos^2(pi/7) and the full-star rank condition; a local Lorentzian six-node diagram alone is not enough. The surviving D5 half-weight configuration is proved to be an entire finite corner component. These essential qualifications are present. I found no source-level omission, but did not rerun those field censuses in this audit.

The checker independently confirms all 18 old and 21 new field/dimension branches. The noncompact arithmetic part is a separate published bound 21, with dimension 20 excluded; it is not deduced from compact h-vectors. Arithmetic lattices generated by reflections are of simplest type, so there is no additional reflection-lattice type omitted by the composition. See [Belolipetsky, introduction and Lemma 2.1](https://arxiv.org/pdf/1506.03111), which cites the original [Esselmann theorem](https://doi.org/10.1006/jnth.1996.0141) and Vinberg's lemma. This does not bound arbitrary arithmetic lattices without the reflection condition.

## 7. Literature status and provenance correction

A bounded current-source check located the official 2026 Mühlherr–Weiss survey. Its Section 7, pages 9–10, still reports the published compact/arithmetic exclusion at dimension 30, Prokhorov's general exclusion at 996, known examples through 19 and in dimension 21, and no known examples at or above 22 or compact examples above 8. [Official primary author publication](https://londmathsoc.onlinelibrary.wiley.com/doi/pdf/10.1112/jlms.70379). This corroborates the comparison used in the checkpoint; it does **not** establish that no other publication proves an equal or stronger bound.

The inspected newer restricted result of Burcroff concerns bounded facet excess, with the consequence that a compact 29-polytope would need at least 40 facets. It is not an unrestricted arithmetic dimension improvement. [Paper, abstract and Section 10](https://arxiv.org/pdf/2201.03437). No novelty/current-best claim follows from this limited search. Nikulin's older announced unrestricted <300 estimate and the repository's nonarithmetic candidates are separate from the arithmetic composition audited here.

The old complete finite-queue independent review did pass generation and a fresh saved replay: session 23308 exited 0 with `PASS_INDEPENDENT_COMPLETE_FINITE_ATTACHMENT_QUEUE_DIMENSION29`, 5,172 states and 2,070 obligations. Its final JSON hash is `d18ea4b16a962595c76bb735494d81466fc70f1c996e5bc0758d3fa0b7b02088`.

One frozen prose statement requires clarification. The dimension-27 coverage note says the fast primary export and independent review both passed fresh saved replay. More precisely, the fast primary had generation plus terminal structural database/map checks; the **independent full reviewer** had generation and fresh saved replay. The final checkpoint describes this correctly. This note records the distinction without modifying the frozen artifact. It does not change the mathematical coverage premise.

The remaining validation limit is explicit: this is an internally reviewed computer-assisted research proof candidate. The present audit adds independent proofs of the listed transitions and small exact quantifier/scalar checks. It inherits the reviewed large catalogues and arithmetic queues, does not formally verify their implementations, and does not claim external peer validation.

```sh
python3 scripts/fresh_arithmetic_dimension23_scope_audit_20261016.py
python3 scripts/fresh_arithmetic_dimension23_scope_audit_20261016.py --verify
```
