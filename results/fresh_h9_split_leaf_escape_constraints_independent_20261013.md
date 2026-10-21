# Independent exact escape-clause review

The current 1,041-clause snapshot passes an independent integer/Fraction replay of **all 3,164,640 cap–escape pairing signs**. It contains 496,882 actual cutter incidences, with no omitted cutters and no extra uncertain cutters. This is a necessary-clause certificate, **not an unsatisfiability or finite-volume nonexistence result by itself**.

The replay uses neither the search producer nor its interval-sign code. It reads the literal exact normals and independently reconstructs the seed Gram matrix over `Q(phi)`. Its inertia is `(9,1,6)`. The sum of the sixteen seed normals has norm `−64−128phi<0`. Every one of the 3,040 saved cap normals has norm two and negative pairing with that common center. For each cap, the replay checks all sixteen covector entries directly against its saved normal and the Gram matrix.

For each saved rational escape, it checks positive norm, negative center pairing, and strict satisfaction of every seed inequality. If the cap covector has entries `0,−1,−phi,−d` and the escape has integer coefficients `w/D`, its cap pairing is

`−(u+v phi+t d)/D`,

where `u,v,t` are the corresponding sums of coefficients. Write `d=d0+d1 sqrt(Delta)`. The independent replay decides the sign in this expression exactly by comparing `A²` with `B² Delta` in `Q(phi)` whenever `A` and `B sqrt(Delta)` have opposite signs. It includes explicit zero-radical and square-radical cancellation tests; no independence assumption about the radical is needed.

The producer's lower-interval test is sound: if the certified lower bound for `u+v phi+t d` is nonnegative, the cap pairing is nonpositive, so that cap cannot cut the escape. Otherwise retaining the cap as a possible cutter is conservative. The independent replay verifies that every actual positive pairing occurs in its archived clause. In this snapshot the conservative and exact sets coincide.

Why a cutter is necessary: if a selected collection cuts no saved escape, that future spacelike vector satisfies all selected inequalities. Join it to the common strict timelike interior. The segment crosses the future light cone at a point strict on every inequality, giving an open set of ideal directions and therefore infinite hyperbolic volume. Hence every finite-volume selection must contain at least one member of each certified cutter clause. A separate certified compatibility graph and Boolean contradiction would be needed to turn these necessary conditions into a finite-bank obstruction.

The audit verifies the exact-normal archive's upstream source hash and records hashes of both the normal archive and **the complete clause archive**. It also rejects a mutation deleting one true cutter. A future clause addition or change requires a fresh audit, even when its normal bank is unchanged. The clause count above describes the audited snapshot, not any later continuation.

Files:

- [Independent replay](../scripts/fresh_h9_split_leaf_escape_constraints_independent_20261013.py)
- [Audited snapshot and source hashes](../data/enumerations/fresh_h9_split_leaf_escape_constraints_independent_20261013.json)
- [Necessary clauses](../data/enumerations/fresh_h9_split_leaf_escape_constraints_20261013.json)

```sh
python scripts/fresh_h9_split_leaf_escape_constraints_independent_20261013.py --verify
```

If the necessary-clause archive has intentionally changed, run the same command without `--verify` to generate a new complete exact audit, then replay it with `--verify`.
