# Independent check of the surviving marked-facet flow

**PASS:** the inherited marked-facet baseline has maximum 38, above the
required threshold 35. The marked `A1+D20` type has a genuine narrow
same-type, same-port local edge at each of its 20 in-facet ports. The
resulting balanced port-flow relaxation therefore still permits cost 38.
These separately realized edges need not fit into one simultaneous star
or one compact facet.

The [reviewer](../scripts/fresh_sqrt5_marked_flow_20261016_review.py) and
[certificate](fresh_sqrt5_marked_flow_20261016_review.json) independently
check the baseline and the 20 selected witnesses in the
[primary transition certificate](fresh_sqrt5_marked_maximum_transitions_20261016.json).
They do not repeat the exhaustive enumeration of all 54,040 attachments or
the primary exclusions at five `D21` ports. Those broader enumerations
retain their primary-only verification scope.

The baseline is reconstructed using breadth-first graph distances and the
two stated H4 adjustments, without calling the producer's weight function.
Every connected-component total and marked-root debit agrees. A two-state
finite-forest optimization, recording whether exactly one component has
already been marked, gives the exact rank-21 maximum 38. The only maximizing
marked components are the isolated A1 with a residual D20, and the long-arm
endpoint of D21. Independent paired h-vector coefficients give the facet
dimension-20 threshold `5*45-binomial(20,2)=35`, with strictness from
`h0=h20=1`.

For each selected edge, the reviewer reconstructs the literal 22-root Gram
from the old `A1+D20` graph, the saved reattachment, and its framing product.
An independent tree code colors the facet root and edge port differently;
it checks that both endpoints have the same marked diagram and port orbit.
The old rank-21 Gram is positive definite. Its last Schur complement is
negative at the physical embedding and positive at the other embedding,
giving signatures `(21,1)` and `(22,0)`. All normalized Gram entries are
integral, and every new contact has an allowed Coxeter value with framing
product strictly between minus 14 and zero.

The 40 saved endpoint rays are checked by direct full Gram multiplication.
Each has exactly the required zero pairings, one strictly negative omitted
wall pairing, and negative square. The two rays have negative mutual product,
so are on the same future sheet. They lie in the marked facet. These checks
verify actual local finite edges, not only algebraically admissible endpoint
diagrams. Giving every same-port edge mass one half supplies incidence one
at each of the 20 ports.

The review uses the project's independent Golden arithmetic and frozen
review matrix operations, but imports neither transition producer nor its
quadratic arithmetic. Generation, fresh saved replay, source hashes and
local report links pass. No global polytope realization or new dimension
bound follows from the surviving relaxation.

```sh
python3 scripts/fresh_sqrt5_marked_flow_20261016_review.py
python3 scripts/fresh_sqrt5_marked_flow_20261016_review.py --verify
```
