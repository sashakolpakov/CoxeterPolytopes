# Independent mathematical review of all-rank4 leaf refunds

18 October 2026. **PASS for the written theorem and its stated scope: no flaw
found.** I read the full new report, the preceding B4/F4/H4 refund report,
and the E7 facewise-refund report. This review involved no proof execution,
census, numerical calculation, optimizer, Lean process or TeX build.

Reviewed new source:
`results/compact29_all_rank4_leaf_refund_20261018.md`, SHA256
`9cf22ed38a0f8b57e5fe0c28ef146f171ad6aa8dff70bb0b5499a6ae73d3cb7f`.

## The potential E8 loophole is closed geometrically

The old A4/D4 is a whole component, so after deleting its leaf b every
remaining A3 node has no old outside contact. In a finite E8 continuation,
the new root c has only one A3 contact. Contact at the middle would force
two terminal arms of length one, incompatible with E8. Endpoint contact
makes the surviving A3 the terminal three vertices of E8's unique arm of
length four. Therefore c is that arm's fourth vertex, and the other four
old vertices constitute a whole A4. This proves the outside placement;
it does not assume those contacts vanish or select an A4 inside a larger
old component.

If the second replacement corner is unpaid, c is orthogonal to b by
adjacency positivity. Together with the original whole-component products,
this gives exactly the diagrams listed in the report. I checked the
neighbor-sum identity for every displayed coefficient: the A4 far-end
case has arms (1,2,5) and coefficients (6;3;4,2;5,4,3,2,1); the near-end
case has two forks at distance two with path coefficients two and leaf
coefficients one; D4 has two forks at distance three with the same rule.
The latter two are actual principal subdiagrams, with the explicitly
omitted vertices removed, rather than claims that the full nine-node
diagram is positive semidefinite.

The positive null combination cannot vanish because its product with a
strict interior point is negative. Its products with all other walls
are nonpositive, so it gives a feasible ideal ray and contradicts
compactness. Thus the proof needs no identification of Vinberg's omitted
E8 pairs at all. It would be incorrect to carry this compactness argument
unchanged to a noncompact chamber.

## Completeness and consumers

I also checked the other continuation branches. Classical endpoint contact
gives a terminal four-node chain; middle contact gives the D fork segment
containing both whiskers. Outside continuation can only pass through c,
so the stated end segments retain the required invariance. The remaining
small whole exceptional components have positive original weight. If c
does not meet A3, two unpaid replacement corners would make all old roots
plus c positive definite, contradicting positive index n. All outside
whole components are retained in this argument.

Every selected A4 pair and every D4 pair contains an eligible leaf. The
simultaneous rule preserves all old zero weights and leaves every old
positive weight at least one. Hence the distinct positive partner supplied
by the lemma remains sufficient even when it too is modified. A changed
quad keeps a full unit, and longer faces remain nonnegative. This covers
all actual face consumers without an unproved quad classification.

The selector sigma0>84 controls every rank4-containing heavy forest since
the minimum component saving is three and sigma0<=87. It leaves the E7
source forests of original charges 80,84,78 unchanged, including the
boundary case 84. The essential-quad implication survives because the
zero set is unchanged; its rank-two source at total 78 is protected.
Therefore the active bound 81 and inactive bound 171/2 still compose.

## Counts, dependencies, and limits

Conditional on the inherited complete 2126-type domain, the remaining
rank partitions give 80+18=98 types: 79 at allowance 87, the single
E7^3+E8 type at 171/2, and 18 at 85. The reduction controls charges;
it does not exclude geometric realizations of the other types. I checked
the stated necessary excess inequality against the rank29 balance:
pairing h_i with h_(29-i) gives minimum coefficient 168, and h_0=1
contributes a further 1050 above 84 times the vertex count.

The original Vinberg rule/face theorem, finite Coxeter classification,
compact acute actual-corner geometry, prior E7 contact theorem and
essential-quad proof, and inherited complete typed domain remain explicit
dependencies. I did not freshly re-establish the primary weighting or
rerun that domain. This is an independent written review, not an
end-to-end formal or external verification. The compact dimension ceiling
remains 29; the 98 remaining types and their incidences are still open.
