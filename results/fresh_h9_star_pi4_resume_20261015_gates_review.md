# Independent scope and transport check of the pi4-third binder

The [separate checker](../scripts/fresh_h9_star_pi4_resume_20261015_gates_review.py) and [certificate](fresh_h9_star_pi4_resume_20261015_gates_review.json) independently reconstruct the literal transport performed by the [new binder](fresh_h9_star_pi4_resume_20261015_gates.json). Generation and saved replay pass. No issue was found.

Each of the twelve saved maps preserves the full induced maximal-face complex on the eleven-node core, and every core pair occurs in an actual vertex. The checker transports the 114 path, kernel and budget support premises through every map and verifies all **1,368 actual-face occurrences** directly against the saved incidence data. Transporting the eight third4 orientations produces exactly **96 distinct predicates**, disjoint from the inherited 397, giving 493 in total. The saved survivor at cursor 60,015 satisfies three of the new forbidden predicates and none of the old ones.

The third4 mathematical exclusion and prior 397 predicates remain dependencies on their independently reviewed certificates. The twelve maps are induced-core maps; no full-nerve automorphism is assumed. This check does not replay the resumed DFS, prove that its pending frontier is empty, or supply a metric realization. A search stop or a newly excluded saved survivor must not be described as a whole-incidence exclusion.

The search source's use of the new gate path and `old_family_count=397` is consistent with this transport. Its search outcome remains separately reported.

```sh
python scripts/fresh_h9_star_pi4_resume_20261015_gates_review.py
python scripts/fresh_h9_star_pi4_resume_20261015_gates_review.py --verify
```
