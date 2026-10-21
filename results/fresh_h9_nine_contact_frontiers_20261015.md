# The new frontiers still miss the latest escape in bounded searches

The exact 30-wall chamber has eight open finite-corner germs at its one new vertex. Their [certificate](fresh_h9_nine_contact_germs_20261015.json) and [checker](../scripts/fresh_h9_nine_contact_germs_20261015.py) cover 117 finite spherical profiles and five symbolic dihedral tails. The older 293 open germs are not repeated. All eight limiting null points have unbounded links.

The [complete-line census](fresh_h9_nine_contact_lines_20261015.json), produced by its [checker](../scripts/fresh_h9_nine_contact_lines_20261015.py), gives 75 complete geodesic one-faces with no finite endpoint. All previous 72 remain. The three new supports are

* `[0,1,3,5,6,8,10,29]`;
* `[0,2,3,5,6,8,10,29]`;
* `[1,2,3,5,6,8,10,29]`.

Of 3,604 relevant exterior-ray pairs, 1,918 fail the clique test and 1,611 fail positive definiteness. The complete ray input is the independently reviewed 30-wall cone. These derived frontiers have primary exact checks; no separate independent frontier implementation has been run at this checkpoint.

The [new-germ scout](fresh_h9_nine_contact_germ_scout_20261015.json) produces 2,312 distinct numerical proposals and no cutter of the latest escape. The [complete-line scout](fresh_h9_nine_contact_line_scout_20261015.json) produces 71,375 proposals, also with no escape cutter. Both bind their complete compressed banks. Their anchor labels and sampled dihedral tails are bounded, and their global pair gates are numerical. These zero counts do not exclude every label or parameter.

This is evidence to change the construction search, not a proof that the 30-wall chamber admits no finite-volume completion. Joint additions whose vertices require multiple new walls, changed retained walls, other sources and other incidence types remain outside these single-endpoint proposal searches.

Generation and saved replay pass for both exact frontiers and both bounded scouts.

```sh
python scripts/fresh_h9_nine_contact_germs_20261015.py --verify
python scripts/fresh_h9_nine_contact_lines_20261015.py --verify
python scripts/fresh_h9_nine_contact_germ_scout_20261015.py --verify
python scripts/fresh_h9_nine_contact_line_scout_20261015.py --verify
```
