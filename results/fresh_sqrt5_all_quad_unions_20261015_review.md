# Independent audit of the all-quadrilateral traversal and resume protocol

**PASS for the stable source, transactional protocol, and two complete left-piece blocks. The full enumeration remains pending.** The [independent checker](../scripts/fresh_sqrt5_all_quad_unions_20261015_review.py) and [certificate](fresh_sqrt5_all_quad_unions_20261015_review.json) do not open, modify, or bind the active producer database. They bind the stable traversal source and its frozen catalogue/arithmetic dependencies. Generation and saved replay pass.

The independent comparison rebuilds all 153 marked minimal pieces from the reviewed Lannér catalogue, using ordered black choices and an independent white-permutation quotient. For every pair of white counts from zero through three, a separate overlap construction chooses shared right slots and injective left targets, then assigns canonical fresh names. Its complete sets agree with the producer's identification generator.

The source-level coverage argument is valid. In a compact quadrilateral star, each opposite black pair lies in a minimal indefinite piece. All such pieces come from the bound catalogue; interchanging the two opposite pairs allows the unordered pair of piece indices. Permuting white names within a piece loses no attachment possibilities, since every shared identification and every unspecified cross-contact is subsequently enumerated. Every unspecified pair belongs to an actual adjacent-corner Gram and therefore has a finite label from the complete arithmetic alphabet. Opposite black pairs are already specified by their pieces, with label zero representing an unrestricted dotted parameter. There is no residual no-short restriction and adjacent black labels are free.

Partial corner failure is a safe prune: adding nonpositive off-diagonal entries to a symmetric Z-matrix cannot restore positive definiteness, by the nonnegative-vector Rayleigh/Perron criterion. The separate finite-Dynkin classifier used here does not use the producer's Schur leaf elimination. Finite-edge cycle pruning is also necessary: a cycle involving only assigned finite edges has a fixed squareclass, and later contacts cannot alter it. Dotted edges are ignored in this test, which retains a necessary relaxation rather than imposing an unproved dotted-distance restriction.

At a completed finite-label union, the physical and conjugate inertia tests are necessary for an arithmetic Lorentzian realization. The independent comparison uses a different tree-product rescaling and full Fraction Q(sqrt5) inertia calculation, rather than the producer's color rescaling and two-by-two Schur test. Dotted rows remain unfiltered except for fixed finite cycles, finite corners and the orthogonal-Lannér obstruction. Rejecting a graph at its first orthogonal-Lannér witness is safe even if it has another piece decomposition: the two orthogonal negative subspaces are present in the full graph itself. All other deduplicated terminal tests depend only on the complete labeled graph.

The independent recursive traversal visits every completion in left blocks 0 and 1, using reverse label order. It agrees exactly with the producer:

| Quantity | Count |
|---|---:|
| Base identifications | 623 |
| Passing base identifications | 305 |
| DFS nodes | 36,820 |
| Distinct completed graphs | 10,237 |
| Dotted-parameter survivors | 9,465 |
| Finite principal-test survivors | 109 |
| Orthogonal-Lannér rejections | 305 |
| Two-negative-direction rejections | 307 |
| Conjugate-PSD rejections | 51 |

Within those completed blocks, all 128 earlier no-short rows are accounted for: 116 dotted survivors, two finite survivors, and ten exact two-negative-direction rejections. This is a bounded containment check, not yet the requested containment check for the whole previous census.

Three fresh temporary SQLite runs test the actual producer protocol. An uninterrupted run, a sequence of closes/reopens at node budgets and a support-storage guard, and a run interrupted immediately after an INSERT following a committed checkpoint produce identical ordered key/status/payload rows, counters and final cursor. The multipause test includes nonempty stacks of depths one, three, four, five, seven and nine. The injected interruption rolls back the uncommitted insertion and resumes from the last atomic result-plus-cursor checkpoint, including a valid empty-stack boundary state.

The frame invariant is explicit: each ancestor's next-label cursor identifies its live child's assigned label; the last frame stores the next unvisited sibling. A terminal row is inserted and its frame advanced before the next commit. Result rows, stack, counters and outer indices are committed in one transaction. A crash between commits therefore loses only replayable work, and the key uniqueness constraint prevents duplicate result rows. Source hashes are checked on reopen. The exported JSON is not the transaction authority; if interrupted during export it may be stale and must be regenerated from the database by resume before saved-output verification.

The producer source is used only for these protocol tests on new temporary databases. The comparison traversal, overlap construction, spherical classification, and full field inertia arithmetic are separate. The catalogue's all-label theorem and previously reviewed arithmetic-field hypotheses remain explicit inherited inputs.

Before claiming a complete all-union enumeration, a later final audit must verify completion of all 153 left blocks, the final cursor/stack and database counters, all saved rows and terminal statuses, and containment or exact rejection of every earlier no-short row/class. Arbitrary external white attachments, rank completion, actual quadrilateral realization of every survivor, and any global dimension theorem are outside this enumeration's scope.

```sh
python3 scripts/fresh_sqrt5_all_quad_unions_20261015_review.py --verify
```
