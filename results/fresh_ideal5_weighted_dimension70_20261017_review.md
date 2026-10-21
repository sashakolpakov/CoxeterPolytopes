# Independent scalar and composition review of the dimension70 proof

The [dimension70 proof](fresh_ideal5_weighted_dimension70_20261017.md) passes the independent exact scalar comparison and the mathematical composition review described here. No remaining logical gap was identified in the reviewed scope. The [review source](../scripts/fresh_ideal5_weighted_dimension70_20261017_review.py) uses only the Python standard library and imports no producer arithmetic. Final hash-bound generation and a separate fresh saved replay both terminated with exit0/PASS; the [saved review certificate](fresh_ideal5_weighted_dimension70_20261017_review.json) binds the frozen primary and all direct inputs.

## Independent exact scalar calculation

At ambient dimension71, the cusp rank is70. The reviewer selects the optional high rank, optional low rank and optional rank12 first, then fills the remaining rank with arbitrary copies of11,13,14,15,16,17. The producer uses its inherited repeated-parts-first enumeration. Exact set equality confirms that both enumerate the same1,432 profiles without duplicates.

For each profile the reviewer multiplies the factors `1+t+...+t^r` by sliding-window convolution. It independently checks the central coefficient sum using the cumulative bounded-composition identity

    G = floor(d/2)*V - 2*sum_shift a_shift*binom(floor(d/2)-shift+s,s+1),

where `sum a_shift*t^shift = product(1-t^(r+1))`. It also verifies the absolute-deviation formula

    curvature/V = 3d/4 - (5/2)*E|sum U_i-d/2|.

All polynomial coefficients, vertex counts, G values, curvature, costs and margins agree exactly with all1,432 primary rows. The minimum margin is **23641/338688**, attained at `(11,11,11,11,13,13)`. The finite-vertex slack is **19/2**, and the inherited global additive slack is350. At ambient70, the boundary profile `(11,11,11,11,11,14)` has margin **−1637/51840**. This negative formal profile establishes neither realizability nor optimality of the theorem.

Only the rank15 coefficient changes, from111/16 to107/16. The unchanged concentration moves preserve or increase cost: low ranks concentrate to11; two12 ranks become11 and13, increasing cost by16/91; and two ranks at least18 concentrate until one is17, preserving the constant tail cost. Integer squared-part sums strictly increase and are bounded. Thus the inherited complete terminal domain remains valid. Repeated rank15 parts are retained. The reviewer also reconstructs the unchanged fifteen piecewise affine cost branches and their endpoint/limit maxima; it does not infer an all-rank maximum from a finite sample.

## Changed local interface and literal checks

The [order16 local transfer](fresh_d16_partial_weight_transfer_20261017.md) changes ten still-paid affine-D15 pairs from1 to1/2. The reviewer reconstructs the canonical16-node graph, its distances and eight automorphisms. The ten pairs are exactly the previously paid pairs at distance greater than9 and form an invariant set. Its complete120-entry weight vector sums to106. The nine inherited zero pairs stay zero, while all94 elementary triangle/short8 pairs retain weight1. No finite coupon is added.

From the frozen [small-factor map certificate](fresh_distance8_ideal5_small_factor_20261017.md), the reviewer independently selects all54 affected cases and compares the complete328 saved maps. Each literal map and marked image is checked against the actual core's labelled Gram. Their66 grouped completion outcomes are:

|Outcome|Literal completions|Maps|
|---|---:|---:|
|Multiple required core-contact roots|20|64|
|Required contact projection square exceeds2|14|112|
|Completed Gram has two negative directions|22|72|
|Already-whole sealed factor|10|80|

Independent Fraction Gaussian elimination computes the saved core inverses. Exact congruence elimination, including two-by-two hyperbolic pivots when needed, recomputes the signatures of the nine distinct full Grams underlying the22 index exclusions and ten sealed completions. It also checks the14 projection exclusions from inverse diagonals. The20 multiple-root outcomes use the independently checked required-contact counts and the local theorem's reviewed at-most-one-contact-root lemma; this review does not claim a new exhaustive classification of supporting-root profiles.

For all ten sealed completions, the reviewer independently verifies the positive affine kernel, the white and four edge principals, all40 corner signatures, their entire finite components and the marked coupon roles. The finite halves are whole A4 against whole D12 at a fork center or short leaf, or whole D5 against whole D11 at a short leaf. Each occurs at a different corner from the selected ideal half. Every total is at least1; four cases have additional finite unit payment. These checks use actual corner occurrences rather than arbitrary finite subdiagrams.

The mathematical/source read checks why the completed Gram's unspecified entries really are zero. The inherited contact and no-bridge lemmas give at most one relevant core-contact root. In the relevant-component cases, the possible extra exterior root would bound the whole marked component by14, so it cannot coexist with the actual order16 factor. Once the forced factor is already whole, no additional supporting white can contact its core. The displayed finite components are therefore whole even with arbitrarily many spectators. This geometric justification is necessary before an exact inertia calculation can exclude a completion.

The original complete core/overlap and induced-map censuses remain explicit frozen inputs. This review compares their saved affected support and the literal completions; it does not rerun those original censuses.

## Complete proof composition

The local modification concerns only distance10 or11 pairs. It leaves short8 units, elementary triangle units and all inherited nine zero-pair replacement proofs intact. The complete prior split into low, disjoint, overlapping, protected/basic and immediate-factor alternatives is retained. Only the54 affected saved contexts change, and each possible completion is either excluded or receives the ideal half plus an existing finite half.

The [finite budget](fresh_finite8_d12_b10_coupon_budget_20261017.md) remains8n through dimension75. Its13,034 endpoint checks and proof were independently reconstructed in the [dimension71 review](fresh_ideal5_weighted_dimension71_20261017_review.md); they are bound and inherited here without repeating that unchanged enumeration. The order16 maximum is now107, attained by affine B, while D has cost106. All other small-order costs and the all-rank tail are unchanged.

For an ideal vertex with V link vertices, a factor of order r+1 contributes its weighted pair sum multiplied by V/(r+1). Fractional weights preserve this counting formula, including at nonsimple cusps. The [shape inequality](fresh_cusp_link_shape_curvature_20261017.md), the positive ideal margins and finite slack then contradict the total face-payment upper bound at dimension71. The completed [dimension71 theorem](fresh_ideal5_weighted_dimension71_20261017.md) excludes all higher dimensions. Together these prove the unrestricted ceiling70 under the stated finite-sided, finite-volume real-hyperbolic Coxeter hypotheses.

This argument adds no arithmeticity, simple-cusp, cusp-factor or facet-count restriction. It is separate from the compact-arithmetic partial models. No new construction, external mathematical refereeing or proof-assistant formalization is asserted.

## Verification status

The independent scalar preparation, complete primary-row/local-interface comparison, final hash-bound generation and separate fresh saved replay all terminated with exit0/PASS. The final commands completed directly, without outstanding process handles. All21 direct source bindings, Python syntax and nine local report links pass. The saved primary hashes are source `83d1d3a33f86a25153f381497071e96d8a885a694306d6eb1bcf0045496e6beb`, data `7f683b9d1a9ce11cec039335d9008b2da50cda380a2f0e0a57a27f6b192795d7`, and report `ebf0d1b19f7f1ef4a0a4486c64cf93cb3fc00895964b99543b142617d4431257`.

The local source and written proof, and the full primary source and composition, received a complete mathematical read in the scopes above. The parent also independently read this reviewer's entire source and its mathematical reconstruction, including scalar enumeration, inverse/inertia routines, all saved maps and corner roles, and the global interface; no gap was identified. That source review does not claim to repeat this reviewer's1,432-profile calculation. No upstream exhaustive census is claimed to have been independently repeated by this reviewer. The source, data and report are final/frozen; all validation is internal.
