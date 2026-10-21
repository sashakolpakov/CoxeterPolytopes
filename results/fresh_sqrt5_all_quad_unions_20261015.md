# All quadrilateral minimal-piece unions over Q(sqrt5)

The complete necessary-union enumeration has finished: **38,755 retained literal supports form 8,298 marked isomorphism classes**, comprising 5,904 finite-contact classes and 2,394 classes with undetermined dotted distances. It removes the no-distance-two restriction from the previous enumeration, supplying the previously omitted quadrilateral constraints needed when a free-weight rule decreases the old distance-two payments. These are necessary union supports, not certified realizable full quadrilateral stars or a dimension bound.

The [terminal certificate](fresh_sqrt5_all_quad_unions_20261015.json) records all 153 completed starting-piece blocks and all 127,875 shared-white identifications. Its terminal cursor is `(left,right,mapping,stack)=(153,153,0,null)`, with `complete=true`. The [SQLite database](../data/enumerations/fresh_sqrt5_all_quad_unions_20261015_running.sqlite) is now immutable; `running` in its filename is historical.

| Quantity | Count |
|---|---:|
| Entered DFS states | 977,040 |
| Label trials | 5,257,464 |
| Corner-PD prunes | 4,252,699 |
| Partial finite-cycle-field prunes | 41,782 |
| Completed labeled leaves, including duplicates | 100,796 |
| Distinct completed labeled graphs | 98,963 |
| Orthogonal-Lannér-pair rejections | 11,781 |
| Two-negative-direction rejections | 25,154 |
| Conjugate-PSD rejections | 23,273 |
| Retained finite-contact supports | 29,290 |
| Retained symbolic-dotted supports | 9,465 |

The [producer](../scripts/fresh_sqrt5_all_quad_unions_20261015.py) retains the frozen 153 marked Lannér-piece representatives and all 127,875 possible shared-white identifications across unordered piece pairs. It inherits the independently reviewed finite alphabet 2,3,4,5,6,10, symbolic dotted pairs, and exact Q(phi) corner-positive-definiteness test. All adjacent black contacts are now free to take every finite label, and white paths of length two between adjacent black nodes are allowed. No old distance-two payment is used as a pruning rule.

Every unspecified pair belongs to a spherical corner, so it receives a finite label. Physical corner failures are monotone when further acute contacts are added. At a partial assignment, the checker also rejects a finite-edge cycle whose squareclass product is outside Q(sqrt5). A dotted edge has unknown squareclass and is ignored by this particular test; cycles involving it are not silently excluded. At a complete finite-contact leaf, the separately reviewed radical gauge and principal Schur tests reject two negative directions or a nonpositive-semidefinite conjugate. Those conjugate tests are not applied to an incomplete matrix with unknown contacts artificially set right. Dotted leaves retain their undetermined distances. Disjoint orthogonal Lannér pieces are excluded by Lorentz index one.

Results are deduplicated by the full literal label array in a SQLite database. The finite-diagram cache is bounded, and each batch has node and retained-support limits. SQLite commits the new rows and the exact unfinished DFS stack, counters and source hashes in a single transaction. Each stack frame records its depth and next label index; its preceding label is the branch leading to the child. Consequently a pause can retain an entered but unprocessed leaf, and all later siblings remain reachable. Only fully processed outer pieces count as completed.

An actual mid-DFS pause at 36,857 visits was saved and validated with seven live frames; its immutable [historical snapshot](fresh_sqrt5_all_quad_unions_20261015_pause_test.json) records that state. Subsequent batches resumed the same database without restarting completed pieces. Thirteen exported batches reached the terminal cursor. The historical snapshot is not a claim that its old database hash remains current.

`--verify` opens the database read-only. It checks SQLite integrity, source bindings, stored counters, every retained literal piece embedding, all corner and finite-cycle conditions, the full finite-contact principal tests, and the live stack. This is a complete saved-row/frontier replay; it does not independently re-enumerate previously rejected branches.

The frozen [independent protocol review](fresh_sqrt5_all_quad_unions_20261015_review.md) passes source-level coverage, separate reconstruction of all 153 pieces and shared-white maps, and a full independent traversal of left blocks 0 and 1. Its separate finite-Dynkin and Fraction-inertia implementation matches all 36,820 visits and 10,237 distinct completed graphs in those blocks. Fresh temporary databases also agree after multiple mid-stack pauses, a storage-limit stop, and an injected insertion followed by rollback. That review explicitly certifies the protocol and two-block traversal; a separate final audit checks the terminal cursor and all earlier no-short supports.

The [independent completion audit](fresh_sqrt5_all_quad_unions_20261015_review_complete.md) checks the full terminal accounting, 155,020 retained corner diagrams with a separate finite-Dynkin classifier, and every quotient map. All 5,440 earlier no-short literal supports are accounted for: 116 dotted and 1,685 finite supports are retained; 2,847 have two negative directions, 350 fail conjugate PSD, and 442 fail the fixed finite-cycle field condition. It independently rechecks those earlier arithmetic decisions. It does not duplicate the complete 977,040-state traversal or independently recompute every new retained principal inertia; the source/protocol proof and primary saved arithmetic replay have their separate stated scopes.

The [marked quotient](fresh_sqrt5_all_quad_unions_20261015_classes.json) is computed by a [separate script](../scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py). It examines all eight black-label permutations preserving the two opposite pairs, and all white permutations within invariant adjacency-signature groups. Each merged member has an explicit full-label bijection to its representative; 1,015,927 pair equalities are checked. The 8,298 class counts by node count 4 through 10 are respectively 143, 792, 1,281, 1,961, 1,903, 1,565 and 653. These are marked graph classes, not polytope incidence classes.

Generation and saved replay pass for the complete retained-row arithmetic and the marked quotient. Syntax, all directly bound source hashes, and report links also pass. No running jobs or unfinished search branches remain in the producer.

```sh
python3 scripts/fresh_sqrt5_all_quad_unions_20261015.py --verify
python3 scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py --verify
```

This classifies the chosen minimal-piece union with necessary principal tests. Dotted-distance feasibility, arbitrary external white attachments, full-star rank and actual face geometry remain separate obligations. For a sufficient weighting proof, imposing a face-weight inequality on every retained union is conservative, but the pair weights must also remain valid under all possible finite-component extensions. No new dimension bound is claimed by this enumeration.
