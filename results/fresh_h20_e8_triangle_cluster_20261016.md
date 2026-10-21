# A numerical E8 attachment with a compact triangular two-face

This package gives an **exact-promotable four-wall recipe**, numerically compatible with all158 retained H20 walls. It creates a compact triangle on the new E8-supported two-face and predicts four new ambient finite corners, three lying on the chosen H3 face. It does not close that H3 face or establish finite volume. Exact reconstruction and interval signs for the four new walls remain outstanding.

The [bounded search](../scripts/fresh_h20_e8_triangle_cluster_20261016.py), [search data](fresh_h20_e8_triangle_cluster_20261016.json), [rational-contact recovery](../scripts/fresh_h20_e8_triangle_cluster_20261016_recipe.py), and [recipe data](fresh_h20_e8_triangle_cluster_20261016_recipe.json) all pass generation and fresh saved replay. They are final for this numerical package. The input [finite attachment census](fresh_h20_spherical_attachment_profiles_20261016.md) is frozen separately.

## Prescribed geometry

Write T for the fixed support `A10+A4+A2+A1`, with paths `[1,...,10]`, `[12,13,14,15]`, `[61,87]`, `[151]`. New root r, index158, has square2 and label3 contacts to12,61,151, and is right to every other root in T. Its T-projection has square59/30; its residual square is1/30. Thus T+r is `A10+E8`.

The three partners, indices159,160,161, also have square2 and are right to r. Their T contacts are respectively only root1 at label3, only root10 at label3, and no nonright contacts. Their mutual products are

`159.160=0`, `159.161=-sqrt(2)`, `160.161=-1`.

Their full Gram is B3. Removing the positive T projections gives the Lorentzian matrix

```
 [ 12/11   -1/11   -sqrt(2) ]
 [ -1/11   12/11   -1       ]
 [ -sqrt(2) -1       2       ].
```

Every two-by-two principal block is positive definite. This gives the compact triangular two-face on support T+r when all ambient gates and coherence are certified. The three T-face corners have types `E8+A12`, `E8+B12`, `E8+A12`.

The full21-node diagram T+four-new splits as E8 and a13-cycle having exactly one label4 edge. All13 vertex deletions of that cycle are tested by exact rational leaf elimination. Exactly four are positive definite: omission of oldroot1 or newroot159,160,161. The first gives a fourth ambient corner, type `E8+B12`, outside the fixed T-face. Removing an E8 node leaves the indefinite cycle, and therefore cannot give a finite corner. Since every old wall outside T is numerically dotted to every new wall, there are no other candidate old/new finite-corner supports. This last ambient conclusion remains conditional on exact verification of those dotted gates.

## Recoverable algebraic recipe

The saved recipe uses only exact finite-contact targets and six rational free old pairings. Every step solves twenty independent linear conditions in the21-dimensional ambient space followed by the square2 quadratic. The numerical branch reference selects the intended real root; no numerical SVD basis is part of the algebraic specification.

| New root | Nonright T contacts | Further exact constraints |
|---|---|---|
|158|12,61,151 at label3|95 product `-75089/10000`;142 `-116551/10000`;145 `-36367/5000`|
|159|1 at label3|right158;95 `-26659/2000`;142 `-209307/10000`|
|160|10 at label3|right158,159;95 `-142043/10000`|
|161|none|right158;159 product `-sqrt(2)`;160 product `-1`|

Unlisted T products are zero. The minimum normalized cosine against old walls outside T is approximately1.37576,3.08452,3.84493,3.08719 for these four roots. Thus the numerical proposal has substantial strict dotted margins. The source records all linear-frame Grams, particular solutions, quadratic radicals, branch choices, and residual coordinates for independent reconstruction.

## Bounded search and remaining boundary

The search tests36 initial E8 rays, eight boosts, and1440 rotations of the prescribed triangular partner geometry:414,720 rotations total. Of these400,964 pass all158 native ambient gates. Every accepted grid index is retained, with one literal maximum-margin representative for each of288 families. The rational recipe recovers the profile0, direction5, boost1.25 representative.

In the same fixed future-center chart used by the input scout, the old H3 maximum squared radius is3.06383257946. The selected four-wall cluster reduces it to2.22840934760. Its numerical cone has24 rays:10 finite,14 exterior, none in the `1e-6` sphere band. There are three new finite vertices on this H3, but exterior boundary remains. The counts concern this face only; they are not a complete ambient H20 vertex census.

All finite gate lists and rotation/boost grids are bounded discovery choices. The source has no optimizer-based or sampled assertion that an untested family is impossible. This package is a qualitatively different attachment experiment because its E8 root and two A10-contact partners are not right to S16. No exact162-wall chamber, new dimension, or finite-volume construction is claimed here.
