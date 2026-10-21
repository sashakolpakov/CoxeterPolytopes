# Five-port constraints coupling simplex blocks and dotted prisms

For five edge ports at a vertex of a compact simple Coxeter polytope, the
selected simplex and dotted-prism incidences satisfy fifteen nonnegative
inequalities certified below. Thirteen were found to separate local
relaxations containing the earlier four-port constraints during discovery;
the saved verifier retains one particularly simple exact separation.
This is a new local constraint family, not an ambient dimension bound.

Let `s_i` count maximal simplex port blocks meeting the five ports in exactly
`i` directions, for `i=2,...,5`. Let `q_j` count selected canonical dotted
full-star prism aggregates whose restricted footprint contains exactly `j`
square pairs, for `j=1,...,4`. One useful inequality is

\[
\boxed{s_3+2s_4+3s_5+q_2+q_3+q_4\le3.}
\]

Another is

\[
s_2+3s_3+4s_4+5s_5+q_2+q_3+q_4\le5.
\]

All coefficients are nonnegative. Consequently, retaining only the simplex
or prism classes represented in a partial model gives valid weaker charges.
For a fixed five-port set, average the whole inequality under the finite
vertex forest's automorphisms. Use the same actual prism provenance in every
coefficient. No symmetry of the ambient chamber is assumed.

## Why the necessary domain is finite and exhaustive

The [maximal-simplex theorem](fresh_simplex_disjoint_ports_20261017.md)
partitions the ports into triangle blocks. The selected dotted aggregates
have star footprints whose leaves lie in one triangle block, as proved and
checked on all frozen templates in the
[four-port package](fresh_sqrt5_joint_four_port_domain_20261017.md).

Consider the selected square pairs between two different triangle blocks.
Two incident square edges have their other endpoints in the same triangle
block. The [triangle-plus-two-squares completion](fresh_prism_simplex_shared_capacity_20261017.md)
forces these squares into the same canonical aggregate. Therefore each
connected component of this bipartite graph belongs to one aggregate and
must be a star. Conversely, the nonempty restriction of one selected
aggregate is connected, so it is precisely one such component. The entire
selected square graph between a pair of blocks is thus a disjoint union
of stars. This assertion concerns each pair of blocks separately; different
block pairs may reuse a port.

Enumerating all 52 labelled partitions and every subset of cross-block
pairs gives 11,985 patterns. The star-component test retains 10,445 and
rejects 1,540. The retained patterns give 74 distinct vectors in the eight
coordinates above. Every actual selected configuration is included.
Retained patterns are not asserted to be geometrically realizable. General
non-dotted products, including a product of two triangles, are outside this
selected-prism scope.

A second enumeration independently reconstructs the same 74 vectors. For
a pair of triangle blocks with sizes `a,b`, distinguish one vertex on the
first side. It is isolated, is the center of a star with `j` leaves on the
other side, or is one of `j>=2` leaves of a star centered on the other side.
Removing that component gives an exact recursion on smaller block sizes.
Combine these profiles over all pairs of blocks and the seven integer
partitions of five. The exact match checks the labelled enumeration against
a different counting method.

## Exact certified inequalities

Each row lists coefficients in `s2,s3,s4,s5,q1,q2,q3,q4` order, followed by
the upper bound. Every row is checked with integer arithmetic on all 74
vectors. The domain has exact affine dimension eight; each row has tight
points of exact affine dimension seven. Thus every listed inequality is a
facet of this necessary invariant domain. Completeness of this facet list
is not claimed.

| Simplex coefficients | Prism coefficients | Bound |
|---|---|---:|
| 0,1,1,1 | 0,0,0,0 | 1 |
| 1,0,1,2 | 0,0,1,1 | 2 |
| 1,1,2,2 | 0,0,0,0 | 2 |
| 0,1,2,3 | 0,1,1,1 | 3 |
| 1,3,4,5 | 0,1,1,1 | 5 |
| 1,7,9,10 | 1,1,1,1 | 10 |
| 3,5,9,10 | 1,1,1,1 | 10 |
| 3,5,8,10 | 1,1,2,2 | 10 |
| 3,6,10,12 | 1,2,2,2 | 12 |
| 4,7,11,14 | 1,2,3,3 | 14 |
| 3,8,12,15 | 1,3,3,3 | 15 |
| 5,10,17,20 | 2,3,3,3 | 20 |
| 5,10,16,20 | 2,3,4,4 | 20 |
| 6,12,19,24 | 2,4,5,5 | 24 |
| 7,22,33,40 | 4,7,7,7 | 40 |

## An exact separation from four-port capacities

Put zero mass on every simplex block, and mass one on each of these five
two-edge prism footprints:

```
01,12    02,24    03,34    13,23    04,14
```

Together they partition the ten pairs of five ports. Every pair capacity
holds with equality. On any triple at most three distinct footprints occur,
so every strengthened `2s+t+q<=3` constraint holds; at most one footprint
occupies two of its pairs, so the older wedge bound also holds. Every
four-port set contains at most three full two-edge footprints, satisfying
`r+q_ge2<=3`. Zero simplex mass has the empty packing cover, satisfying the
entire simplex-only packing cone.

But on all five ports the boxed inequality has load five against capacity
three. This is an exact local relaxation witness, not a realizable chamber.
In particular, it is not claimed to satisfy the global coface flags or a
common arithmetic-form condition. It proves that the new local family adds
information beyond the listed pair, simplex, three-port and four-port
constraints.

## Reproduction and scope

Floating convex-hull calculations proposed the fifteen facets, and fifteen
small local linear programs explored separation. The saved
[standard-library checker](../scripts/fresh_sqrt5_joint_five_port_domain_20261017.py)
uses neither method: it re-enumerates the complete necessary domain twice,
checks all facet inequalities and ranks exactly, and directly verifies the
displayed witness. No global optimization or saved-dual scan occurs here.
Generation and a fresh saved replay pass. The
[certificate](fresh_sqrt5_joint_five_port_domain_20261017.json) stores every
labelled pattern, all invariant vectors, tight points and witness checks.

```sh
python3 scripts/fresh_sqrt5_joint_five_port_domain_20261017.py --verify
```
