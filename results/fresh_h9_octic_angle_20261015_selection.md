# Numerical selection from the 152 saved escape cutters

The [diagnostic](../scripts/fresh_h9_octic_angle_20261015_selection.py) and [saved replay](fresh_h9_octic_angle_20261015_selection.json) examine exactly the 152 cutters in the frozen [29-wall rescore](fresh_h9_after_fourth_rescore_20261015.md). They find **no compatible pair**, so the largest compatible subset in this numerical graph has size one. A single renewed strict escape also survives the union of all 152 proposed walls with the 29 retained walls, even when their mutual incompatibility is ignored.

All 11,476 unordered pairings are evaluated from a 100-decimal physical frame, using the saved floating-point proposal coordinates. Of these, 11,362 are positive and 114 negative; none is within `1e-7` of zero. Their range is approximately `[-1.786888970686, 2.550169731077]`. The closest finite-angle fit is still `0.0007855558603466` away, at a positive near-right pair. There are no dotted or parallel pairs.

The finite-label scan checks labels 2 through 256, but a separate numerical monotonicity margin excludes the entire larger-label tail: every mutual pairing exceeds `-2cos(pi/7)` by at least `0.01504876511922`. Thus only labels 2 through 6 could occur in these numerical values. This is not an exact all-label certificate for the saved approximate proposals.

The literal renewed escape is embedded in the script and JSON, with no temporary-file dependency. Its numerical square is `0.02007686339762>0`; its center pairing is negative. Its largest pairing with the old normalized walls is `-1/1024`, and its largest pairing with the 152 cutters is approximately `-0.2051007955668`. The same point is [independently certified exactly](fresh_h9_octic_escape_caps_20261015_review.md) for the two selected octic caps and the 29-wall base; the remaining 150 signs here remain numerical. The reconstructed proposal norms differ from two by at most `1.20e-12`.

Generation and saved replay pass. Source hashes bind the frozen chamber, rescore and compressed proposal bank. The result concerns this bounded bank only. It does not exclude caps on new supports, other parameter values, different contact profiles or changed retained walls, and it supplies no new wall to a chamber.

```sh
python scripts/fresh_h9_octic_angle_20261015_selection.py
python scripts/fresh_h9_octic_angle_20261015_selection.py --verify
```
