# Independent source and geometric review of the five-port domain

**PASS: no gap identified in the stated necessary-domain argument.** This is a separate source/proof read by exact_hub_counts, not a duplicate enumeration, facet discovery, arithmetic replay, or saved-dual test.

Reviewed frozen artifacts:

| Artifact | SHA256 |
|---|---|
| [Source](../scripts/fresh_sqrt5_joint_five_port_domain_20261017.py) | `348e80947555b1a840d91091b8c1a73b4603a53b601395f5e0d25841e63cd77f` |
| [Certificate](fresh_sqrt5_joint_five_port_domain_20261017.json) | `e8ce83355fe29e41311ab0dae9fd96da70b6707c45affc82738d9ed57ac05f01` |
| [Proof/report](fresh_sqrt5_joint_five_port_domain_20261017.md) | `df64ff2917a2aee6f826e3336a9bd461206296b76c2999ae4e3cefe2d32f0c67` |

For a fixed pair of distinct maximal-simplex blocks, any two incident selected square edges have their other endpoints in a common triangle block. Triangle-plus-two-squares completion therefore identifies their canonical aggregate. Connectivity propagates this identification through each bipartite component. A selected dotted aggregate has a star footprint, so a component must be a star. Conversely its nonempty footprint is connected, its leaves are in one simplex block, and its hub is in one other block. Thus it belongs to exactly one block pair and one component.

This argument is correctly restricted to canonical dotted-prism aggregates. A general non-dotted product, such as a product of two triangles, need not satisfy the star condition. Distinct block pairs may reuse a vertex; treating their star choices independently enlarges the necessary domain and does not remove an actual configuration.

The labelled enumeration covers all partitions and cross-pair subsets. The separate component recursion is exhaustive: a distinguished vertex is isolated, is a star center with at least one opposite leaf, or is a leaf of an opposite center with at least two leaves. The one-edge case is already included in the center alternative. Removing the chosen component yields the smaller pair problem. The recursion tracks invariant vectors, not labelled multiplicities; that distinction matches the source's set comparison.

The saved checker uses rational affine ranks and integer slacks on the complete proposed invariant set. Its claims are properly limited to the fifteen displayed valid facets, not completeness of the facet list or realizability of every retained vector. Numerical convex-hull discovery is not trusted by the saved verifier.

The five-star fractional witness is also correctly scoped. Its stars partition the ten square pairs, it has zero simplex mass, and the source directly checks every pair, triple and four-subset inequality claimed. It need not obey the global coface flags or arithmetic-form conditions, which the report explicitly excludes from the separation claim.

For a subsequent scan, the simplex coefficients must be evaluated on the **exact intersection sizes**. Some rows are nonmonotone in block size, so replacing them by threshold indicators would be invalid. Exact-size histogram averaging and one common actual quad provenance supply the stated valid charge normalization.
