# Physical compatibility of previously unowned covering pairs

This successor tests a complete, small new physical-pair domain at the unchanged face masses of the [current historical weighted certificate](fresh_current_historical_weighted_corners_20261017.md). In five selected finite vertex forests, 2,250 of 2,286 previously unowned pairs of literal prism-corner stars are incompatible. Five exact new weighted columns lower the repaired full-witness ratio from 12.335925889381436 to **12.332859482381691**. The named-Q restricted ratio falls from 12.058423241708674 to **12.055550455551511**. Both remain above the dimension-13 target 12. These are exact repairs of specified witnesses in a partial model, not a feasibility result or dimension theorem.

The [source](../scripts/fresh_unowned_cover_pair_geometry_20261017.py) and [saved certificate](fresh_unowned_cover_pair_geometry_20261017.json) preserve every frozen predecessor. No global optimization is performed.

## Domain and ownership

The preceding complete weighted calculation exposed 275 positive covering terms in 117 forests whose chosen literal representatives use at least one never-owned pair. We selected the five leading forests by their total such term mass. The domain here is **every** never-owned pair in each selected forest, including pairs absent from those particular representatives. Thus selection of a forest does not truncate its new pair domain.

| Vertex forest | New pair tests | Incompatible | Retained unknown |
|---|---:|---:|---:|
| 982 | 1,752 | 1,734 | 18 |
| 2524 | 272 | 268 | 4 |
| 2190 | 216 | 204 | 12 |
| 3324 | 6 | 6 | 0 |
| 3165 | 40 | 38 | 2 |
| Total | **2,286** | **2,250** | **36** |

The five forests contain 255 distinct literal stars. Each is referenced by its frozen source and literal index. The source rebuilds the predecessor union, checks equality with the saved union, and revalidates current explicit or unique provenance, corner multiplicity, raw-core order and spectator rank through the frozen preparation routine. A prism of dimension d has 2d actual corners; its face-obligation value is d(d−1)/2. These two quantities are kept separate.

A historical source owns a pair precisely when it contains both complete native star keys. Source order is applied before inspecting the outcome. The new domain consists only of pairs contained in none of those source domains. Inherited unknown outcomes therefore cannot be overridden. A future consumer must resolve the frozen predecessor's first-owner domains first, and then consult this package's explicit `newly_owned_pair_domain`. All pairs outside these five domains retain their predecessor outcome.

## Actual wall products

At a finite simple vertex, the thirteen incident norm-2 wall normals form a positive-definite Gram basis for the spacelike complement of the vertex vector. The saved complete native contact covector determines an external wall's spacelike projection s. With the future unit timelike vertex vector v, the outward actual wall is

\[
n=s+\sqrt{s^2-2}\,v.
\]

The positive square root is forced by the strict side of the finite vertex. Consequently two actual external walls with projection product p and squared timelike heights h₁,h₂ have product

\[
n_1\mathbin{\cdot}n_2=p-\sqrt{h_1h_2}.
\]

The source reconstructs each finite forest's exact active-component inverse Gram matrices. It recomputes both positive heights for all 255 stars and checks agreement with the frozen values. Each star's two external covectors also reproduce its native omitted-wall product exactly. Passive components have no contact support, as directly asserted. Thus both stars are interpreted in one complete labelled vertex frame, with their common future signs fixed; their orientations are not chosen independently after evaluating a product.

The exact pair-field routine selects at most two nonsquare native-contact classes over K=Q(√5), with positive K coefficients for proportional squareclasses. Every projection is recomputed in that immutable selected field. Cached keys are native contact vectors, never coefficient vectors interpreted in another radical basis. Identical external covectors identify the same outward wall and are skipped as an incompatibility witness.

The necessary exclusions are:

| Reason | Pair exclusions |
|---|---:|
| Distinct complete stars occupy the same actual quadrilateral port pair | 508 |
| Positive product of distinct actual walls | 792 |
| Product strictly between −1 and 0 | 134 |
| Nonintegral squared product in K | 374 |
| Squared product outside K | 442 |
| Additional forbidden finite Coxeter product | 0 |

The first reason uses uniqueness of the actual two-face through the marked pair of ports and its complete quadrilateral wall star. The next two use acuteness and the norm-2 Coxeter angle gap. The arithmetic reasons use the fixed ambient arithmetic ground field K and its necessary two-cycle integrality. The general exact gate is inherited from the reviewed [raw110 package](fresh_raw110_mixed_actual_corners_20261017.md); it treats both square roots already in the selected field and the independent-root case, testing physical signs before squaring. The certificate stores 138 distinct exact wall expressions underlying the exclusions. Failure to find an exclusion leaves a pair compatible in the conservative graph; it does not assert physical realizability.

## Exact weighted separation and repair

A simultaneous set of actual corner stars is a clique in the conservative compatibility graph. Complete maximal-clique enumeration in the five revised graphs gives 268 literal cliques and 53 nondominated orbit-count vectors. The occurrence coefficients are inherited actual prism-corner multiplicities, not the number of orientations in an orbit.

For each of the full and restricted demands, a nonnegative orbit charge w has capacity 1 when w·c≤1 for every maximal count vector c. Its exact load is the sum of the corresponding actual face masses with their corner multiplicities. The saved covering mixtures dominate every orbit marginal and have total mass equal to the charge load. Thus they certify the exact optimum over all nonnegative orbit weights for each of these five conservative graphs at the specified masses. There are 48 full-demand and 40 restricted-demand mixture terms. The replay reconstructs the graphs and all maximal count vectors, then checks these saved rational certificates without another optimization.

Only five distinct charge columns violate one of the two prior vertex budgets. Four full-witness vertex masses and three restricted-witness masses increase; all face masses remain unchanged. Forest 2190 has a larger exact graph-cone optimum but its prior vertex budget already covers it. The unchanged 545 forests retain their same-mass exact graph-cone certificates. All 439,410 inherited finite-model columns, including 96,176 signed zero-budget columns and 52 signed positive-vertex-cost columns, and all 78 preceding weighted columns remain valid: the variables are nonnegative, so their dual constraints are inequalities; holding all face masses fixed and increasing nonnegative vertex budgets only relaxes them.

The exact full ratio is

\[
\frac{5064933545222915}{410686066151853}
=12.332859482381691\ldots,
\]

with vertex normalizer 136895355383951/136500000000000. The named-Q restricted ratio is

\[
\frac{18670750401222348}{1548726494908789}
=12.055550455551511\ldots,
\]

with normalizer 1548726494908789/1544400000000000. “Named-Q” means the same frozen 248-class named domain, of which 209 classes are positive. The historical weighted domain has 248 positive classes, including 39 outside that named set; these are different sets. This successor adds no face class and no new mass restriction.

## Verification and limits

Generation session 18439 and a separate fresh saved replay, session 84267, both terminated with exit 0 and PASS. The fresh replay rebuilds the native projections, all new pair outcomes and the revised graph domains; it verifies saved exact charges and mixtures without running a global or local optimizer.

A separate standard-library Fraction saved-data check, with no producer imports, passed the complete ownership-domain enumeration, all exclusion references, all 53 charge inequalities, all 88 mixture terms, all five face-coefficient loads, both monotone repairs, exact ratios, nine source hashes and Python syntax. This check does not independently reconstruct physical projections, field gates or maximal cliques; those belong to the producer replay.

The parent independently read the complete new source and report and found no gap in this delta, inheriting the previously reviewed field gate and provenance preparation. A separate construction agent read the complete new source/provenance/mathematical delta and also found no gap. These reviews cover ownership, common-frame heights and future signs, native-key caching, conservative capacities and monotone repair; they do not claim an independent physical-pair census, field-arithmetic reconstruction, graph enumeration or local optimization replay.

This is complete weighted separation for the five revised conservative graphs, together with the 545 unchanged predecessor graphs. It is not a complete physical compatibility classification, a new-alpha simplex/joint-family certificate, or a complete global arithmetic model. The 36 new unknown pairs and every untouched historical or unowned pair remain conservatively allowed. No dimension claim follows from these ratios.
