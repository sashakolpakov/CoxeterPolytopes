# Next tree-incidence targets: literal gate eligibility only

**Case 0 is the recommended different-structure target.** It is the remaining star incidence, with 378 vertices and 132 distinct conditional gates from the three new theorem families: 24 dense/high, 72 strong-contact and 36 two-dense conditions. Its old necessary-label survivor was never a complete metric; subsequent exclusions of specified assignments do not exclude this whole incidence. No new label search was run here.

For maximum gate count alone, the T-type cases 6 and 20 each have 155 conditions; case 11 has 154. These counts measure available pruning predicates, not predicted runtimes or likelihood of a metric survivor.

The [checker](../scripts/fresh_h9_dotted_central_20261014_next_targets.py) and [certificate](fresh_h9_dotted_central_20261014_next_targets.json) examine all 33 frozen representatives. The ten closed cases are **1, 2, 3, 4, 5, 8, 17, 22, 25 and 26**; 23 remain eligible. Their statuses are inherited from the bound completed proofs. No closed case is recommended for restart.

Each placement uses an ordered pair A,B of the four declared two-facet wedge blocks and distinct outside leaves x,y. Both own triples A+x and B+y must be literal minimal nonfaces. Every required finite support is checked against a containing **saved actual vertex**, and every thresholded pair has its own such witness. This does not use an incidence isomorphism or merely assume that the case5 supports transfer.

The three inherited all-label theorems are the [dense/high gate](fresh_h9_dotted_central_20261014_gate.json), [strong-contact gate](fresh_h9_dotted_central_20261014_strong_contact_gate.json) and [two-dense gate](fresh_h9_dotted_central_20261014_dense_pair_gate.json), with their separate proof reviews. Their threshold hypotheses remain conditional. A placement does not force those labels to occur. All unlisted contacts and dotted distances retain the freedom stated in the original theorems.

The table gives distinct threshold-condition counts. The dense/high column also shows its literal placement count in parentheses, since different unspecified leaves can produce the same condition. For the other two families the placement and distinct-condition counts agree in this catalogue. All witnesses and conditions are saved.

| Case | Tree | Vertices | Status | Dense/high (placements) | Strong | Two dense | Total |
|---:|:---|---:|:---|---:|---:|---:|---:|
| 0 | star | 378 | eligible | 24 (36) | 72 | 36 | 132 |
| 1 | star | 432 | closed | 21 (30) | 60 | 30 | 111 |
| 2 | star | 432 | closed | 19 (28) | 56 | 28 | 103 |
| 3 | star | 432 | closed | 18 (26) | 52 | 26 | 96 |
| 4 | star | 432 | closed | 17 (24) | 48 | 24 | 89 |
| 5 | T | 444 | closed | 24 (44) | 88 | 44 | 156 |
| 6 | T | 488 | eligible | 23 (44) | 88 | 44 | 155 |
| 7 | T | 488 | eligible | 21 (42) | 84 | 42 | 147 |
| 8 | T | 462 | closed | 20 (38) | 76 | 38 | 134 |
| 9 | T | 480 | eligible | 20 (36) | 72 | 36 | 128 |
| 10 | T | 444 | eligible | 22 (42) | 84 | 42 | 148 |
| 11 | T | 488 | eligible | 22 (44) | 88 | 44 | 154 |
| 12 | T | 488 | eligible | 20 (42) | 84 | 42 | 146 |
| 13 | T | 480 | eligible | 21 (38) | 76 | 38 | 135 |
| 14 | T | 420 | eligible | 24 (40) | 80 | 40 | 144 |
| 15 | T | 456 | eligible | 23 (36) | 72 | 36 | 131 |
| 16 | T | 464 | eligible | 20 (36) | 72 | 36 | 128 |
| 17 | T | 438 | closed | 21 (36) | 72 | 36 | 129 |
| 18 | T | 456 | eligible | 20 (32) | 64 | 32 | 116 |
| 19 | T | 444 | eligible | 22 (42) | 84 | 42 | 148 |
| 20 | T | 488 | eligible | 23 (44) | 88 | 44 | 155 |
| 21 | T | 488 | eligible | 20 (40) | 80 | 40 | 140 |
| 22 | T | 462 | closed | 21 (40) | 80 | 40 | 141 |
| 23 | T | 480 | eligible | 22 (38) | 76 | 38 | 136 |
| 24 | T | 396 | eligible | 24 (36) | 72 | 36 | 132 |
| 25 | T | 432 | closed | 20 (30) | 60 | 30 | 110 |
| 26 | T | 432 | closed | 19 (28) | 56 | 28 | 103 |
| 27 | T | 414 | eligible | 23 (34) | 68 | 34 | 125 |
| 28 | T | 432 | eligible | 19 (28) | 56 | 28 | 103 |
| 29 | T | 488 | eligible | 22 (42) | 84 | 42 | 148 |
| 30 | T | 488 | eligible | 21 (40) | 80 | 40 | 141 |
| 31 | T | 462 | eligible | 23 (42) | 84 | 42 | 149 |
| 32 | T | 480 | eligible | 23 (38) | 76 | 38 | 137 |

Case0 uses star edges (0,1), (0,2), (0,3), (0,4), cut in that order, with negative coordinate0 unwedged. This is structurally different from the recently closed T-type case5. The recommendation is preparation for a future search; it is not a new survivor or realization.

As a regression, the new literal enumeration reproduces all case5 maps and condition sets exactly: 44/24 dense-high placements/conditions, 88 strong-contact conditions and 44 two-dense conditions.

Generation and saved replay pass. The bounded table exhausts only the declared-block placement scheme; it does not claim to enumerate arbitrary embeddings using other facet pairs. It performs no label propagation, backtracking, optimization or metric solving.
