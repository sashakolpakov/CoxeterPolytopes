# Global H5 closure test: both finite banks are excluded

The search tested whether the saved cap banks could close the entire distinguished H5 face, rather than adding more isolated endpoints. The result is an exact finite-bank obstruction, including optional old-wall deletion. No new walls were appended to the authoritative [125-wall checkpoint](fresh_h9_h5_closure_next_20261013.md), whose H5 face still has 26 facets, 29 finite vertices, 46 finite edges and 53 open edges.

There are two certificates:

- [Original bank](../data/enumerations/fresh_h9_global_h5_20261013_verify.json): all 38,209 literal formula caps share one strict future spacelike escape with the 125 old walls. Compatibility is irrelevant to this first obstruction.
- [Enlarged bank](../data/enumerations/fresh_h9_global_h5_20261013_enlarged_verify.json): adding 119 opposite quadratic branches still permits no finite-volume pairwise Coxeter-compatible subset. Four further escapes and 1,198 exact incompatibilities prove this, even if optional old walls are deleted.

These are literal entry counts, including repetitions, not counts of distinct mirrors. The theorem does not exclude other normals, other profiles or a different seed.

## The original global relaxation

Twelve saved roof, truncation, conic and endpoint banks contribute:

| Bank category | Entries |
| --- | ---: |
| Original roof, changed roof and roof profile proposals | 2,100 |
| Earlier finite truncation and truncation proposals | 1,010 |
| Three four-contact conic banks | 23,229 |
| Earlier endpoint bank | 244 |
| Recent 46-edge endpoint bank | 11,522 |
| Two subsequent coupling banks | 104 |
| **Total** | **38,209** |

The [numerical global scout](../scripts/fresh_h9_global_h5_20261013_scout.py) intersects all these halfspaces simultaneously, ignoring every conflict. It immediately found a surviving spacelike point. The [exact verifier](../scripts/fresh_h9_global_h5_20261013_verify.py) reconstructs the exact formula representatives and proves a common strict escape with norm

`21120203161130719809369943401 / 4000000000000000000000000 - (5221189072472564933543152633 / 1600000000000000000000000) phi`,

approximately 0.0124271888374. Every cap pairing is strictly negative, with the largest certified upper bound below -0.00204116. Projection into the fixed H5 face remains spacelike, with norm above 0.01242718848. Therefore no subset of the original bank can close either chamber.

Each representative is defined by nine exact prescribed pairings, norm two and a specified branch of the residual quadratic. The nine-contact spans have only 121 distinct patterns. The [source schema](../scripts/fresh_h9_global_h5_20261013_schema.py) and unique-row guards cover every row of every bank. The rational decimal snapshot selects the sign of the one-dimensional complement; the checker certifies that sign. Rounded numerical normals are not themselves treated as exact roots. The three additional anchor formulas are also entries of the original bank.

The checker verifies every exact formula branch, without assuming that every proposal has already passed all Coxeter pair gates. Including unusable candidates only strengthens a common-escape obstruction. Rational cosine bounds cover the high-label source proposals; the norm and finite-contact equalities follow from the defining projection equations. A [source review, full saved replay and separate small Fraction audit](fresh_h9_global_h5_20261013_review.md) pass.

## Opposite branches and coherent replacements

Reversing the quadratic branch produces 1,502 numerical cutters of that common escape. Of these, 119 pass the numerical gate against the retained original seed and fixed face walls. Each requires removal of optional old walls. Testing all 119 single-branch replacements showed severe global deterioration; none was adopted as a new chamber.

The enlarged *conflict-ignoring* relaxation was more encouraging: a [complete numerical H5 vertex census](../data/enumerations/fresh_h9_global_h5_20261013_enlarged_geometry.json), including a dual-hull boundedness test, finds 122 timelike vertices and maximum norm approximately -0.00786453. This is a geometric relaxation containing incompatible walls, not a Coxeter polytope or an exact compactness result.

The next step therefore conditioned on each opposite branch and allowed every wall judged compatible with it, while ignoring all remaining mutual conflicts and permitting old-wall deletions. Every one of the 119 resulting numerical supersets still had an escape. These failures compress to four witnesses, which the [supplemental verifier](../scripts/fresh_h9_global_h5_20261013_enlarged_verify.py) checks exactly:

| Witness | Positive norm, approximately | Complete cutter set | Opposite branches assigned |
| --- | ---: | ---: | ---: |
| 0 | 23.01987783 | 15 entries | 68 |
| 1 | 3.14523764 | 54 entries | 2 |
| 2 | 94.99597177 | 1 entry | 28 |
| 3 | 344.84776485 | 2 entries | 21 |

For every witness the checker classifies all 38,453 entries: 125 optional old walls, 38,209 original proposals and 119 opposite branches. All 153,812 signs are strict. For each assigned branch, every possible cutter conflicts with it. The 1,198 conflicts comprise 449 positive pairings and 749 pairings strictly between consecutive finite Coxeter-angle values. These intervals exclude **all integer labels**; the argument uses no finite label cutoff. The largest adjacent labels needed are 28 and 29.

The covering argument is short. If a selection uses no opposite branch, the original common escape survives. Otherwise choose any selected opposite branch and its assigned witness. Every wall that could cut that witness is incompatible with the chosen branch. All selected halfspaces are therefore strict at a future spacelike point. Deleting old walls cannot remove this escape. Keeping the four named fixed face walls also gives an escape within the H5 slice, since the projection only decreases pairings with the other compatible normals.

## Replay and consequence

```sh
python scripts/fresh_h9_global_h5_20261013_verify.py --verify
python scripts/fresh_h9_global_h5_20261013_enlarged_verify.py --verify
```

Both saved replays pass. The [enlarged-bank review](fresh_h9_global_h5_20261013_enlarged_review.md) also passes a separate Fraction check of the witness projections and cosine-gap bounds. It distinguishes those independent calculations from the full replay of the primary arithmetic engine. Both primary checkers use rational interval arithmetic and frozen source hashes, with no optimizer in the proof.

This changes the next search decision: neither selecting more carefully from the original banks nor admitting these 119 opposite branches with coherent old-wall deletion can yield finite volume. A successful continuation needs normals outside this enlarged literal dictionary. The result is not an all-profile completion obstruction and does not settle compact H9 existence.
