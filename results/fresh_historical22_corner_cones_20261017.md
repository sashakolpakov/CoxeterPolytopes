# Complete capacities after adding the22 historical corner classes

The unchanged5,053-face witness now satisfies the complete nonnegative orbit-weight cone of the enlarged conservative corner graph in all553 represented forests. The75 changed forests receive new exact certificates;478 inherit their frozen same-demand certificates. There are no open forests. The new17 capacity inequalities require17 vertex-mass increases and lower the full fixed-witness ratio from12.3328594824 to **12.315398448566656**. The named restricted witness remains **12.055550455551511**. Both exceed12; this is a fixed-mass repair, not an optimized finite-model result or a dimension theorem.

The [source](../scripts/fresh_historical22_corner_cones_20261017.py) and [saved data](fresh_historical22_corner_cones_20261017.json) extend the frozen [22-class provenance/domain diagnostic](fresh_historical22_corner_domain_20261017.md). That diagnostic supplies six actual raw25 classes and16 actual raw39 classes from their historical certificates, increasing the positive represented domain from248 to270 classes. All22 added classes lie outside the named restricted Q domain. This package introduces no new physical wall-pair tests.

## Frozen ownership and unchanged demands

The new source reconstructs each of the75 changed domains and compares the entire result against the frozen diagnostic. Actual corner multiplicities, complete native-star orbits, current face IDs and both demand vectors must agree exactly. Source membership determines the first owner of a pair, before its outcome is consulted: the original eight complete historical domains come first, the explicit tested-pair successor comes next, and the two newly added historical domains come last. Earlier unknown outcomes cannot be overridden by later exclusions. All unowned pairs remain compatible.

The [unowned-pair successor](fresh_unowned_cover_pair_geometry_20261017.md) is the starting vertex-mass vector for both repairs. Its changed forest2190 is outside the75 newly affected forests and retains its explicit-pair certificate. The other477 unchanged forests inherit the [eight-source weighted certificates](fresh_current_historical_weighted_corners_20261017.md). Exact equality of both full and restricted face-mass lists is checked before using any inherited certificate. No previous cover at different face masses is reused.

The75 changed domains contain8,699 literal stars and601 corner orbits, representing601 distinct actual corner occurrences. Three of the forests,207/374/433, are new;72 were already represented. Together with the478 inherited forests this gives553 complete forest domains.

## Complete count vectors and exact weighted capacities

Every simultaneous actual set of corners is a clique in its conservative compatibility graph. The source enumerates all maximal literal cliques, counts their members in each complete corner orbit, and removes only coordinatewise dominated count vectors. Nonnegative orbit weights ensure that this last reduction preserves every relevant capacity inequality.

Across the75 forests, the exhaustive enumeration yields **2,985,031 maximal literal cliques**,613 distinct orbit-count vectors and572 nondominated vectors. Forest143 alone has2,594,484 maximal literal cliques,12 distinct count vectors and8 nondominated vectors. Its explicit guard is three million cliques and30seconds; all other forests retain the inherited one-million-clique/30second guard. A guard never licenses a truncated certificate: it records an open forest. The final generation completes every enumeration, with zero open forests.

For each forest and each of the full/restricted demands d, the saved certificate contains nonnegative charges w and nonnegative mixture coefficients lambda satisfying

    w · v <= 1                         for every maximal count vector v,
    sum(lambda_v * v) >= d             coordinatewise,
    sum(lambda_v) = w · d.

The first inequality proves a valid normalized capacity on actual corners. The second is a covering certificate for the conservative graph cone. Equality proves the exact optimum for all nonnegative orbit charges at these fixed demands. Symmetry averaging uses the complete actual corner orbits, and face coefficients count actual corner occurrences rather than literal orientation multiplicities. The source also checks each saved representative against the compatibility graph and its stated orbit-count vector.

The changed forests use390 full-demand mixture terms and357 restricted-demand terms. Generation uses bounded local LP proposals only; all final inequalities, coverage and equality are checked with rational arithmetic. Fresh `--verify` reconstructs the domains and complete count vectors, then checks the saved rational certificates without invoking an optimizer. No global LP is run.

## Exact repair and preservation

For each forest, the new vertex mass is the maximum of its frozen starting mass and the exact new covering optimum. Only genuinely violating charge columns are exported;17 distinct columns suffice here. The largest increases occur at forests534,975 and500. There are no restricted-domain increases.

The full repaired normalizer is

    4021290490212741 / 4004000000000000,

and the exact full ratio is

    445714151979616520 / 36191614411914669
    = 12.315398448566656… .

The restricted normalizer and ratio remain

    1548726494908789 / 1544400000000000,
    18670750401222348 / 1548726494908789
    = 12.055550455551511… .

Every face mass is fixed and every vertex mass is nondecreasing. Consequently all439,410 inherited finite-model inequalities remain valid, including96,176 zero-budget signed columns and52 signed face columns with positive vertex costs. Their primal variables are nonnegative; signed face coefficients do not turn the corresponding dual inequalities into equalities. Every previously exported corner capacity and weighted constraint is preserved by the same monotonicity. This package does not claim implicit simplex or joint-port coverage at the current face masses beyond the inherited statements that actually apply.

## Validation and scope

Final generation82499 and a separate fresh saved replay42002 both terminated with exit0/PASS. A separate no-producer Fraction checker also terminated with exit0/PASS, independently checking both charge systems against all572 count vectors, both marginal-cover systems and equal objectives, the747 mixture terms, all17 exported loads and repairs, the478 inherited-certificate references, and all17 direct source hashes. It does not re-enumerate the compatibility graphs or reconstruct upstream physical pair arithmetic.

A separate construction-agent read of the complete source and report passed the domain comparison, inheritance, charge/cover accounting, occurrence coefficients, export and monotone-preservation arguments, including the explicit larger guard for forest143. It did not repeat graph enumeration, local LPs or certificate arithmetic. A root mathematical/source read also passed the complete source and report. Separately, the root checked all75 saved charge/cover certificates,747 mixture terms, both charge systems against572 vectors,17 source hashes and the exact repairs using standard-library Fraction arithmetic, without re-enumerating graphs or reconstructing physical pair data. The larger guard changes only the permitted exhaustive enumeration size; it does not alter any graph or certificate criterion.

Completeness here concerns the stated270 positive historical classes, their553 conservative graph domains and all nonnegative orbit weights at the fixed full/restricted demands. Compatible graph cliques need not be physically realizable. Untested wall pairs, additional face classes, other geometric compatibility conditions and other ground fields are outside this conclusion. The remaining ratios above12 are concrete obstructions to these fixed-demand graph constraints alone; no full-model feasibility or new dimension ceiling follows.
