# Joint port facets, cusp fibers and restricted dimension bounds

This successor preserves the [preceding checkpoint](sharp_finite_edge_flags_and_packing_20261017.md).
It establishes two restricted dimension bounds and a sharp mixed-triangle
inequality. Joint simplex/prism constraints strengthen the arithmetic
model, but its restricted obstruction remains above the target. The
construction search produces no improved chamber or new dimension.

The unrestricted finite-volume ceiling remains **163**, with the complete
internally checked proof described in the [status record](dimension163_proof_status_20261017.md).
The compact arithmetic ceiling21 is unchanged and is separate from the
already settled noncompact rational case. No lower unrestricted ceiling
is claimed here.

## Two restricted dimension bounds

The [affine-exclusion proof](fresh_simple_cusp_affine_exclusion_20261017.md)
and the [distance15 proof](fresh_simple_cusp_distance_bound_20261017.md) give:

| Hypothesis on a finite-sided, finite-volume hyperbolic Coxeter polytope | Conclusion |
|---|---:|
| Every ideal vertex is simple | Dimension at most61 |
| At least one facet is compact | Dimension at most127 |

The second assertion permits nonsimple cusps elsewhere and does not assume
that the compact facet has Coxeter angles. The first includes mixed
finite/ideal vertices; the stronger known compact and simple all-ideal
results remain separate. No publication priority or current-best claim is
made for these restricted bounds. The first strengthens the frozen
intermediate bound125; that earlier proof remains preserved.

If every ideal vertex is simple, every connected affine original-wall
subdiagram must be the entire diagram of an ideal vertex, with exactly n
nodes. Its positive kernel produces a feasible null ray, which is an
actual ideal vertex by finite volume and strict convexity of the Klein
ball. The simple Euclidean link is an irreducible affine simplex diagram.
For n>=10 this excludes every small affine subdiagram in a quasi-Lanner
core. The minimal indefinite cores used for polygon payments are therefore
Lanner, of order at most5. Actual simplex faces then have dimension at most4.

Triangles receive two corner paths of length at most3; quadrilaterals
receive one of length at most7. All-rank finite/affine graph formulas bound
the total distance7 corner cost by7n per vertex. The strict h-vector
inequality excludes n>=62, giving61. Generation and fresh replay pass;
the [independent review](fresh_simple_cusp_affine_exclusion_20261017_review.md)
checks the global affine argument, ideal quadrilaterals, all ranks and the
strict parity threshold. This conclusion does not extend to nonsimple cusps.

Closed opposite sides of an actual quadrilateral are disjoint even when
some vertices are ideal. Their original-wall blocks are indefinite. In a
minimal indefinite subdiagram, deleting either opposite black node gives
a positive-definite diagram, because the remaining walls support an actual
edge with hyperbolic interior. The complete frozen catalogue supplies384
marked pairs in26 higher-order quasi-Lannér types. Every node is within
distance7 of the marked pair. Smaller diagrams have radius at most3.

Two opposite-pair diagrams overlap or are joined, since orthogonal
indefinite subspaces cannot coexist in a Lorentzian space. This supplies
a path of length at most15 between consecutive sides. Triangles have two
corner paths of length at most8. Weighting corner pairs at distance at
most15 pays every polygon's demand5−k. Explicit all-rank finite/affine
graph formulas give vertex budget15n. This was the intermediate bound125
for a simple chamber, now superseded above. For a compact facet of dimension
d=n−1, it instead excludes d>=127, giving n<=127. The previous
compact-facet branch143 is improved; the nonsimple-cusp branch163 remains.

Generation and fresh saved replay pass for both packages. Root and exact_hub independently
reviewed the geometry, all-rank coverage and strict parity thresholds;
the [separate review](fresh_simple_cusp_distance_bound_20261017_review.md)
states its scope. This is internal review, not external refereeing or
proof-assistant formalization.

## A sharp mixed-triangle bound at cusps

Write T_j for triangular two-faces with j ideal vertices. At an ideal
vertex v, let its link be a product of simplices of dimensions r_vi,
and let a_v count edges to finite vertices. The
[cusp-fiber theorem](fresh_cusp_finite_port_triangle_capacity_20261017.md)
proves

\[
T_1+2T_2\le\sum_v a_v\sum_i\min(r_{vi},5).
\]

With simple cusps this gives the sharp inequality `T1+2T2 <= 5 E_fi`.
Within one factor fiber, the cusp is simple in the corresponding actual
face. A finite-ended port has a simple finite other endpoint. Balinski
connectivity completes its triangular neighbors to an actual simplex.
Distinct blocks have disjoint finite-ended ports, while their ideal-ended
ports may overlap. The existing bounds on every actual simplex face give
payment at most5 per finite-ended port. The classical H9 simplex with
seven finite and three ideal vertices attains equality.

The argument retains the unexcluded one-finite-vertex six-simplex case.
It avoids the previous T2 counterexample because its two triangular link
edges occupy different factor fibers. The factor sum cannot be dropped
at nonsimple cusps. The combined residual budget still requires control
of all-ideal triangles and quadrilaterals before a lower unrestricted
ceiling follows.

## Complete local simplex packing and new joint constraints

The [restricted local-cone calculation](fresh_sqrt5_restricted_flag_local_cones_20261017.md)
reconstructs all1,864 positive simplex/edge rows, including397 absent from
the earlier export. All593 contributing forests and3,485 block orbits are
covered. Explicit rational mixtures certify every nonnegative local
simplex-block charge:524 values are exact optima and69 have exact brackets.
The [independent cover checker](review_sqrt5_restricted_flag_local_cover_20261017.md)
verifies2,830 mixture terms and12,278 literal disjoint blocks. The
restricted witness remains at12.36702231545>12.

The successive joint scans retain the same2,219 positive face masses,
including355 actual selected dotted-prism provenances:

| Added family | Exact inequalities | Violated columns | Repaired ratio |
|---|---:|---:|---:|
| [Three-port capacity](fresh_sqrt5_joint_triangle_square_ports_20261017.md) |44,260|24|12.35526807578|
| [Four-port facet](fresh_sqrt5_joint_four_port_diagnostic_20261017.md) |90,619|47|12.35033157101|
| [Fifteen five-port facets](fresh_sqrt5_joint_five_port_diagnostic_20261017.md) |2,153,595|115|12.34494473291|

The [four-port theorem](fresh_sqrt5_joint_four_port_domain_20261017.md)
couples simplex-block rank and prism stars. The
[five-port theorem](fresh_sqrt5_joint_five_port_domain_20261017.md) uses the
fact that selected square components between any two triangle blocks
must be disjoint stars. Its complete necessary domain has11,985 labelled
patterns,10,445 accepted patterns and74 invariant vectors. A separate
recursion reproduces these vectors. Fifteen inequalities are verified
exactly as facets; geometric realizability and completeness of the facet
list are not claimed. An exact local witness separates a new inequality
from every previous pair, simplex, three-port and four-port capacity.
The [independent domain review](fresh_sqrt5_joint_five_port_domain_20261017_review.md)
checks the scope and enumeration logic.

Every scan averages the whole inequality over forest automorphisms and
uses common actual provenance. Exact-size histograms handle nonmonotone
five-port simplex coefficients. Only vertex masses increase, preserving
all signed flags and complete simplex-packing covers. The
[standalone Fraction union review](review_sqrt5_joint_port_repairs_20261017.md)
recounts every saved load, per-vertex maximum and all186 critical columns
without producer imports or an optimizer. The final exact obstruction is

\[
\frac{66908102699388132000000000000000}
     {5419878674792788039126804769429}
=12.344944732908832\ldots>12.
\]

These are fixed-witness diagnostics, not new optima or feasible primal
rules. No new global solve was run. The preceding partial-model optimum
bracket remains12.793956036279 to12.793969847724.

The next proposed family couples different vertices of a maximal simplex:
for fixed S and immediate coface U, selected prism bases have disjoint
vertex sets. Pointwise vertex–coface capacities do not exhaust their
integer packing constraints. Three edge bases in a triangular S cannot
coexist in one U, although three half-mass bases satisfy all individual
vertex capacities. Root, construction_audit and exact_hub checked this
scope. Labelled use requires simultaneous coface/subset automorphism
orbits and full spectator multiplicities. Any changed-alpha repair will
require new cover and joint-charge validation; existing covers cannot
silently be reused.

The [first exact integer-base diagnostic](fresh_sqrt5_coface_integer_base_scout_20261017.md)
has now checked the inexpensive coface-summed floor(k/t) family. All7,504
columns pass at the unchanged face masses;1,604 have positive loads and
there are no violations. It reconstructs all6,016 bases of the1,858 saved
actual provenances. The ratio therefore remains12.34494473291. The
simultaneous coface/subset-labelled family and full packing cone remain
untested; repeating the unlabelled scout is unnecessary.

## Construction: changed geometry, unchanged open boundary

The control remains the exact64-wall H9 chamber with198 finite vertices,
15 target-face open ends and certified infinite volume. The
[release47 four-cycle test](fresh_nonarithmetic_release47_cycle_patch_20261017.md)
loses two old finite-vertex certificates and the12 quadrilateral. All8,113
numerically coherent four-cap collections retain15 open ends.

The [changed12 family](fresh_nonarithmetic_release47_replace12_frontier_20261017.md)
replaces the rigid12–52 angle pi/6 by pi/3,pi/4 or pi/5. The sole new
native blocker is53, whose exact timelike projection is redundant on the
target face and whose old finite-vertex incidence is empty. Deleting47
and53 and replacing12 touches four old finite-vertex certificates and
the8 heptagon and12 quadrilateral. These costs remain in the ledger.
The38 seeds give1,138 numerically coherent three-cap collections, all
with15 open ends.

The [nonzero-A2 fourth-cap test](fresh_nonarithmetic_release47_nonzero_a2_bridge_20261017.md)
covers27 new frames. All1,553 collections passing its numerical screen
also retain15 open ends. One approximate near-pi/391 contact is explicitly
not an exact angle certificate and fails the required finite-corner test.
No candidate improves the complete boundary. New geometry remains
numerical discovery data; inherited deletion/projection bookkeeping is
exact. Other coordinated wall families remain open.

## Reproduction and archive

Named producer packages and reviews are frozen after generation and fresh
saved replay, with the scope of each review stated separately. The
[successor auditor](../scripts/audit_joint_port_facets_and_cusp_fibers_20261017.py)
and [saved audit](joint_port_facets_and_cusp_fibers_20261017_audit.json)
preserve the preceding archive's frozen sources, reports and top artifacts,
check recursive manifests including prepared model sources, and check
local links. Hash caching avoids rereading large inputs for each reference.
Archive integrity does not replace the mathematical proofs, exact
certificates or bounded numerical-search scopes.

```sh
python3 scripts/audit_joint_port_facets_and_cusp_fibers_20261017.py --verify
```
