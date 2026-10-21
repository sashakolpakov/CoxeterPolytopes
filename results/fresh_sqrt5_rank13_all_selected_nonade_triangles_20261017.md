# All non-ADE selected endpoint triangles in rank 13

The [source](../scripts/fresh_sqrt5_rank13_all_selected_nonade_triangles_20261017.py) applies the frozen exact endpoint method to every one of the 363 non-ADE members of the selected set. That set contains all 578 rank-13 ADE/H3/H4/I2(5) forests with the old baseline at least 12. The [complete literal catalogue](fresh_sqrt5_rank13_all_selected_nonade_triangles_20261017.json) contains:

- 2,090,517 endpoint contact tests and 151,857 accepted finite attachments;
- 911,142 framing tests, of which 237,388 have the required two-place Schur signs;
- 2,848,656 third-corner tests and 173,434 accepted triangle occurrences;
- 147,234 triangles whose three vertex forests all belong to the selected set.

The accepted framing-label counts are 64,504 for label 2, 61,210 for 3, 8,292 for 4, 33,404 for 5, 1,394 for 6, and 4,630 for 10. No sampled angle or dotted-distance cutoff is used: a triangular face has finite black-pair contacts, and these six labels exhaust squared contact products in Q(sqrt5).

Endpoint completeness uses the finite classification. The new root has degree at most three, only 3-bonds at degree three, and at most one 5-bond at degree two. Every allowed contact subset and assignment is tested, followed by full finite-forest recognition. The original finite Gram and its conjugate are positive definite. A fixed nonzero path between the old port and new root forces the framing product itself into Q(sqrt5); without that path, its square suffices. The exact physical negative and conjugate positive Schur residual is checked for every retained framing. Every possible second deletion is tested against the full finite classification, so the third vertex may lie outside the selected set.

Together with the earlier [ADE endpoint triangle census](fresh_sqrt5_rank13_high_ade_triangles_20261017.md), this covers every triangle with two endpoints in the selected set: if either endpoint is non-ADE, orient that endpoint as the new starting type; otherwise the earlier ADE census applies. This does not cover triangles with fewer than two selected endpoints, general quadrilaterals, or every rank-13 vertex forest. It is a necessary local domain, not a simultaneous compact realization or a dimension bound.

Generation and a separate fresh complete saved reconstruction both pass. Source, exact data, and this report are now frozen.
