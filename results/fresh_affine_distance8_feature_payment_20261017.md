# Affine distance8 plus feature-near distance9 pair budget

Let D be a connected affine Coxeter diagram with m nodes. Mark as features every vertex of degree at least3 and both endpoints of every nonsimple bond (label at least4, including the parallel bond). Give an unordered pair weight1 if its graph distance is at most8, or if its distance is9 and at least one endpoint is within14 graph steps of the feature set. Optionally give every pair weight1 when m≤17. Then the number P(D) of paid pairs satisfies

\[
P(D)\le\min\{8m,\tbinom m2\}.
\]

The rule is intrinsic under every weighted-diagram automorphism. The [source](../scripts/fresh_affine_distance8_feature_payment_20261017.py) and [certificate](fresh_affine_distance8_feature_payment_20261017.json) give exact graph counts and verify the formulas below. This is an affine graph-cost lemma only. It does not prove that the weights pay every triangle or quadrilateral, and it does not assert dimension97 or any other new dimension bound.

## Small diagrams and completeness

For m≤17, even paying every pair gives binom(m,2)≤8m. Thus the optional whole-component rule is safe without any further graph assumption. The certificate additionally checks all63 standard connected affine catalogue entries of order at most17: affine A1; the available A/B/C/D families; E6,E7,E8,F4,G2. It records each literal labelled graph and exact pair-distance histogram.

The standard affine Coxeter classification is an inherited input, as in the frozen [classical graph-count package](fresh_simple_cusp_distance_bound_20261017.md) and [affine endpoint catalogue](fresh_h20_parabolic_endpoint_profiles_20261017.md). In orders m≥18 only the classical A/B/C/D families remain. The finite-rank checks corroborate the arbitrary-rank counting proof below; they are not a cutoff argument for completeness.

## Arbitrary-rank counts

For m≥18 the exact counts are:

| Affine family | Distance≤8 | Extra feature-near distance9 | Total |
|---|---:|---:|---:|
|A cycle|8m|0|8m|
|C path|8m−36|min(m−9,32)|min(9m−45,8m−4)|
|B single fork and doubled far end|8m−35|2+min(m−11,31)|min(9m−44,8m−2)|
|D double fork|8m−34|4+min(m−13,30)|min(9m−43,8m)|

Every total is at most8m. Since m≥18 implies8m≤binom(m,2), this proves the asserted budget in these ranks.

For an affine A cycle there are no features. Each distance k=1,...,8 contributes exactly m unordered pairs because2k<m. No distance9 pair is added.

For the tree families, number a spine with M vertices by0,...,M−1. Its short-pair count is sum(k=1,...,8)(M−k)=8M−36. A spine pair at distance9 is uniquely(i,i+9), with0≤i≤M−10, giving M−9 positions. The first min(a,L) positions and last min(b,L) positions in a list of length L have union size min(L,a+b): if disjoint their sizes add, and otherwise their union is the entire list. This identity also covers a or b larger than L.

For C the whole graph is the spine, M=m. The doubled terminal bonds make0,1,M−2,M−1 features. The near-feature vertices satisfy i≤15 or i≥M−16. Therefore a distance9 pair(i,i+9) is added exactly when i≤15 or i≥M−25: the first16 and last16 positions. This gives min(m−9,32) extras.

For B, M=m−2 and two simple leaves are attached to the left spine endpoint0. The far spine bond(M−2,M−1) is doubled, so the features are0,M−2,M−1. Near-feature spine vertices satisfy i≤14 or i≥M−16. A distance9 spine pair is added when i≤14 or i≥M−25: the first15 and last16 positions, giving min(m−11,31). Each of the two left leaves adds its pair with spine vertex8, for two more distance9 pairs. At distances≤8 the leaves contribute16 leaf-spine pairs and one leaf-leaf pair, hence8M−36+16+1=8m−35.

For D, M=m−4 and two simple leaves attach to each spine endpoint. The features are0 and M−1. Near-feature spine vertices satisfy i≤14 or i≥M−15. A distance9 spine pair is added when i≤14 or i≥M−24: the first15 and last15 positions, giving min(m−13,30). The four leaves each have one distance9 pair along the spine, adding4. Opposite-end leaf pairs have distance M+1≥15 and contribute neither short nor extra payment. At distances≤8 the leaves contribute32 leaf-spine pairs and two same-end leaf-leaf pairs, hence8M−36+32+2=8m−34.

The deficits from8m are max(45−m,4), max(44−m,2) and max(43−m,0), respectively. Affine D therefore attains equality for every m≥43. The radius14 is sharp for this uniform rule: radius15 on affine D of order100 pays802 pairs, exceeding800. This diagnostic is also checked exactly.

Dividing the bound by m gives the per-factor form P(D)/m≤min(8,(m−1)/2). Applying it to actual cusp incidence multiplicities requires the separate geometric counting argument; no transfer is inferred merely from this graph estimate.

## Exact checks and review scope

Generation and fresh saved replay both terminated exit0 with `PASS_ALL_RANK_AFFINE_DISTANCE8_FEATURE_PAIR_BUDGET`:444 breadth-first classical checks (orders18–128 in all four families),63 small affine catalogue checks and450 endpoint-interval checks. Every extra pair is counted in its literal spine-spine, leaf-spine or leaf-leaf category and compared with the formulas. Frozen labelled graph constructors are reused without alteration.

All six direct source hashes, Python syntax and four local report links pass. Frontier_literature independently read the mathematical proof and source and found no gap in the endpoint windows, leaf contributions, arbitrary-rank formulas, small-order reduction or sharpness diagnostic. Root separately read the complete source/report and inherited labelled constructors and confirmed the same counts and per-factor division. These are narrow mathematical/source reviews, not duplicate444/63/450 arithmetic replays. No optimizer is invoked. This is a complete affine-factor cost bound; face payment and ambient dimension conclusions remain outside its scope.
