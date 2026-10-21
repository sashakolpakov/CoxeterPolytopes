# One finite anchor and one varying dotted contact

The [source](../scripts/fresh_h20_one_anchor_dotted_families_20261016.py) and [data](fresh_h20_one_anchor_dotted_families_20261016.json) enlarge the preceding one-parameter cap families. Only wall 151 and one actual old side of the tracked H3 face have prescribed finite angles. A different old-wall contact has a continuous dotted cosine c; fixing c leaves a conic parameter t. The saved base projection, slope and kernel permit later joint variation of (c,t).

The pilot uses 151 labels 2 or 3; genuine old anchors `{95,142,145,150,152,156,157}` with labels `{2,3,4,5,6,8,10}`; and dotted-cosine slices `{1.01,1.1,1.5,2,3}`. Among 3,920 patterns, 1,357 numerical interval families survive. Every surviving interval is split at every pairing zero against all 61 saved test points, giving 18,569 checked representatives. This retains narrow coverage cells within each tested slice. It does not cover all c values and is not a rigorous interval certificate. One numerically unstable norm sample was rejected explicitly.

Of these representatives, 15,044 have one actual old finite neighbor and 3,525 have none. Saved points 16, 35 and 39 have no cutter in this sliced bank. Consequently the bank cannot close the sampled face by itself.

The [mixed acute cover source](../scripts/fresh_h20_one_anchor_dotted_families_20261016_cover.py) and [data](fresh_h20_one_anchor_dotted_families_20261016_cover.json) combine these representatives with the preceding 42,203 roots. They start from all 121 saved acute test points and stop at the explicit 30-iteration guard. The chosen clusters still have degree-one facets and non-Coxeter mutual contacts. No new exact cap, closed face or finite-volume chamber is claimed.

Both source/data pairs passed generation and a fresh saved replay; direct source hashes and syntax are checked separately at freeze. The next finite target is a coordinated edge-endpoint construction, allowing a new wall to meet two new walls without a prescribed finite old anchor.
