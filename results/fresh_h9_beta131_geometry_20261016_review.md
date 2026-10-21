# Independent complete cone census for the beta131 pilot

The [checker](../scripts/fresh_h9_beta131_geometry_20261016_review.py)
enumerates all **220 nine-wall supports** of the twelve-wall height-20
prefix, using its own exact Gaussian elimination. It imports no cone
enumeration code and therefore independently checks the producer's
incremental double-description census.

The [certificate](fresh_h9_beta131_geometry_20261016_review.json) confirms
**22 finite vertices, zero ideal rays, and 28 hyperideal rays**, with no
opposite-sheet causal rays. Every incident-facet set agrees with the
[primary geometry](fresh_h9_beta131_pilot_20261016_geometry.json).

Completeness follows because the wall covectors span rank ten and there
is a strict interior point. The cone is full-dimensional and pointed,
and each extreme ray consequently has a rank-nine active support. Every
such support is included, including supports of nonsimple extreme rays;
duplicates are identified by their full active sets and checked for
projective equality.

An exact spacelike point lies strictly behind every wall and has negative
pairing with the timelike interior point. The segment joining these
points crosses the future null cone in the strict interior of every wall
halfspace, giving an open ideal boundary patch. Thus the **bounded prefix
has infinite volume**. Later mirrors may cut this patch; no
nonreflectivity statement about the whole lattice follows.

Generation and saved replay both pass:

```sh
python3 scripts/fresh_h9_beta131_geometry_20261016_review.py
python3 scripts/fresh_h9_beta131_geometry_20261016_review.py --verify
```

The [separate independent prefix review](fresh_h9_beta131_pilot_20261016_review.md)
checks root enumeration and the norm alphabet's mathematical premises.
