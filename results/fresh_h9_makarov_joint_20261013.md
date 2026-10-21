# Joint mixed-contact caps on the compact H5 face seed

The exact result is a **27-facet infinite-volume Coxeter chamber in H9**. Its 16 new caps jointly turn all 40 previously open edges from the ten compact-face vertices into finite edges. All ten original vertices survive, and the 40 first endpoints are distinct actual finite vertices. An exact strict spacelike escape still proves infinite volume. This is a changed joint construction, not a compact H9 example or a new dimension record.

The chamber retains the eleven walls of the [compact H5 face plus H4 seed](fresh_h9_makarov_face_extension_20261013.md). It does **not** retain the previous twelfth, changed-center cap. The final selection uses mixed contacts to the H4 walls and therefore goes beyond the previously excluded forty vertex-centered, single-weight families.

## Exact checkpoint

- [Exact normals, pair gates, forty endpoints and escape](fresh_h9_makarov_joint_20261013_exact.json).
- [Exact producer and saved replay](../scripts/fresh_h9_makarov_joint_20261013_exact.py).
- [Source/proof review](fresh_h9_makarov_joint_20261013_exact_review.md).
- [Separate arithmetic audit](fresh_h9_makarov_joint_20261013_exact_independent.md).

The form has signature (9,1) over \(K=\mathbb Q(\sqrt5)\). Each new normal is stored as \(A+B\sqrt D\), with both coefficient vectors and positive radicand in \(K\); two normals already lie in \(K\). The normal recipes prescribe nine independent pairings and solve the norm-two quadratic. Branch proposals come from the numerical search, but all subsequent assertions use exact arithmetic.

All 351 pairs pass:

| Relation | Count |
|---|---:|
| Right angle | 161 |
| Angle pi/3 | 31 |
| Angle pi/5 | 12 |
| Strictly ultraparallel | 147 |

The inherited center remains strictly inside every halfspace. Every wall is an actual facet: its orthogonal projection of this timelike center lies on that wall, remains timelike, and is strictly inside every other inequality because all off-diagonal pairings are nonpositive. Thus the count is 27 actual facets, not just 27 inequalities.

For each old finite vertex \(v\) and each released H4 wall \(j\), let \(d_j\) be its H4 dual vector. The old edge is \(v-t d_j\), with the other eight seed walls identically zero. Each cap gives a linear intersection parameter. The producer selects a proposed first positive parameter and then checks the resulting point exactly:

1. Its norm is strictly negative, with the correct time orientation.
2. Its zero set is exactly the eight old edge walls and one new cap.
3. Every other one of the 27 inequalities is strictly negative.
4. The old rank-eight Gram block is positive definite and its extension has a positive Schur complement.

These checks certify firstness independently of the numerical ordering used to select the proposal. The 40 incident sets are distinct. This gives at least 50 actual finite vertices, not a complete vertex census.

The integer escape is

    (140532998, 178861908, 130353897, 81870848, 33432956,
     12307295, -284558239, -345319008, -230212893, -115106557).

Its square is

\[
276679291647808184-123399368968947096\sqrt5>0.
\]

Its pairing with the strict center and all 27 wall normals is strictly negative. Consequently an open set of spacelike directions survives; the chamber has infinite volume despite closure of the original edge frontier.

## Bounded proposal search

All material in this section is numerical scouting. Its counts are not exact distinct-mirror counts or exhaustive cap classifications.

The first bank prescribed five old-wall pairings, at most two nonright among them, and one H4 weight direction. The proposal labels were 3,4,5,6,10. It produced 11,961 rounded entries. A whole-bank LP scout found a spacelike point even when incompatible entries were imposed together, so no compatible-selection search was justified for that bank alone.

The changed bank used the [complete old-edge spherical profile census](fresh_h9_makarov_frontier_20261013.json). It took the 1,193 finite profiles, sampled the 56 isolated dihedral tails, and imposed a ninth finite anchor chosen from the remaining seed walls. Sampled anchor labels were 2,3,4,5,6,8,10,12,16,24. The resulting 13,140 rounded proposals include contacts to two or three H4 nodes.

| Numerical stage | Saved outcome |
|---|---|
| Combined banks | 23,841 rounded entries after numerical deduplication |
| Whole combined bank, 100 LP objectives | No spacelike proposal found; this does not certify boundedness |
| Greedy compatibility search, 300 selections | Reached 39 of the 40 original null ends |
| Direct endpoint-cover search | Found a coherent cover; exactification distinguished ideal from finite contacts |
| Counterexample-guided selection, 30 rounds | Every proposed coherent selection had a remaining spacelike gap |
| Final restricted-field selection | Sixteen caps with strictly positive old-null cuts; all 40 first endpoints verified finite |

The final exact selection restricts prescribed pairings to labels 2,3,5, and any prescribed contact with the unequal-norm polar wall is right. This keeps every linear recipe over \(K\); it is a practical arithmetic choice, not a completeness claim.

Relevant proposal artifacts:

- [First bank](fresh_h9_makarov_joint_20261013.json) and [whole-bank gap](fresh_h9_makarov_joint_20261013_volume.json).
- [Mixed endpoint-profile bank](fresh_h9_makarov_joint_20261013_profiles.json) and [combined-bank scout](fresh_h9_makarov_joint_20261013_profiles_volume.json).
- [Exact-friendly selection recipe](fresh_h9_makarov_joint_20261013_cover_simple.json).
- [Thirty-round joint selection trace](fresh_h9_makarov_joint_20261013_cegar.json).

## Reproduction and next obligation

The exact saved replay is independent of rerunning the numerical optimization:

    python3 scripts/fresh_h9_makarov_joint_20261013_exact.py --verify
    python3 scripts/fresh_h9_makarov_joint_20261013_exact_independent.py --verify

The first command reconstitutes the selected exact recipes, all pairs, center, forty endpoint points and escape. It binds the frozen seed and numerical recipe files by hash. The second uses separate arithmetic for the normals and pair geometry; the linked reviews describe their precise scopes.

The next construction task is to control the new chamber's remaining projective rays and outgoing edges. Covering the old forty null endpoints was a necessary frontier test, but the explicit new escape shows why it cannot serve as a global compactness criterion.
