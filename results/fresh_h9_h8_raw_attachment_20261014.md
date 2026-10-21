# The fixed perpendicular H8 attachment is excluded in every cap direction

**No real all-dotted cap direction completes the frozen mixed twelve-wall H9 chamber to finite volume by attaching this particular compact H8 source with its perpendicular sidewalls alone.** The source is Bugaenko's eleven-facet H8 chamber. The result covers every isometric placement, every cap direction and all finite Coxeter labels. It uses no search cutoff.

The [composition checker](../scripts/fresh_h9_h8_raw_attachment_20261014.py) and [certificate](fresh_h9_h8_raw_attachment_20261014.json) bind an exhaustive division into four cases. Normalize old wall normals and the cap to square two, and call a projected squared norm *low* when it is below `B=sqrt(5)-1`.

| Number of low centers | Obstruction | Independent review |
|---|---|---|
|Two or more|[The two-anchor theorem](fresh_h9_h8_corner_anchors_20261014.md) excludes every pair of low centers.|[Separate exact and interval check](fresh_h9_h8_corner_anchors_20261014_pair_review.md)|
|One, at an old wall 0–10|[The radial theorem](fresh_h9_h8_global_cap_gap_20261014_radial.md) leaves too few balls to cover the other two exceptional source corners.|[Separate physical-coordinate review](fresh_h9_h8_global_cap_gap_20261014_radial_review.md)|
|One, at old polar wall 11|[The four-root projection bound](fresh_h9_h8_corner_anchors_20261014_polar_gate.md) forces cap square greater than `2651/918>2`.|[Separate polynomial and interval review](fresh_h9_h8_corner_anchors_20261014_polar_review.md)|
|None|[The source-diameter argument](fresh_h9_h8_all_high_20261014.md) reduces coverage to three old dotted pairs, each contradicting the cap norm.|[Separate algebraic-field and interval review](fresh_h9_h8_all_high_20261014_review.md)|

The earlier [all-real source contact-gap theorem](fresh_h9_h8_facet_cover_20261014.md) remains an explicit dependency: every low center lies on one of three specified source-vertex rays. Its full gap census is not repeated by the new reviews.

For the raw attachment, finite volume requires the twelve old-wall balls to cover the placed source. The proof contradicts that necessary condition. Additional walls not perpendicular to the cap can supply further covering balls and are outside the theorem. Changed retained walls and other source chambers are also outside its scope. This is neither a general obstruction to compact H9 Coxeter polytopes nor a proof that arbitrary larger completions of the old chamber are impossible.

The earlier [exact cap proposal](fresh_h9_h8_global_cap_gap_20261014.md) passed scalar size tests; its [specific radial rejection](fresh_h9_h8_corner_anchors_20261014_radial_gate.md) already showed why those tests were insufficient. The present theorem subsumes that rejection and the old continuous-line exclusion. Further raw placements of this same source cannot improve the outcome.

```sh
python3 scripts/fresh_h9_h8_raw_attachment_20261014.py --verify
```
