# Raw357 quadrilaterals share a typed simplex-coface capacity

For every one of the47 raw357 source classes, assign a quadrilateral to
its canonically forced compact4-simplex S. If the full wall star of S
has b whole A2 components outside its five-node Lorentz core, then

    number of assigned raw357 quadrilaterals <= 2b * number of such S.

The coefficient is4,6 or8 on this complete family. There is no native
orientation multiplier and no factor for the two ends of the base edge.
This gives a zero-vertex-cost transfer between face counts, rather than
another allowance for markers at a vertex.

The [producer](../scripts/fresh_raw357_simplex_coface_transfer_20261018.py)
exports a [complete typed incidence certificate](fresh_raw357_simplex_coface_transfer_20261018.json).
Sixteen target simplex rows already exist in the current finite model;
six of the resulting signed inequalities strictly violate its current
exact dual. The other31 target simplex types are listed but not installed.
No matrix is changed, no repaired certificate is asserted, and no
optimization or dimension conclusion follows from this packet.

## Intrinsic simplex, edge and outgoing facet

The frozen [raw357 incidence proof](fresh_raw357_extra_vertex_incidence_pilot_20261018.md)
establishes actual finite vertices, compact-simplex equality and recovery
of the source quad. Its canonical ten-wall core has a unique quad with
opposite pairs {1,3} and {0,2}. There is no nontrivial native automorphism.
After adding the complete rank-five spectator, this is the entire
15-wall star of that actual quadrilateral.

Delete wall2. The canonically assigned S has fixed walls
{0,7,8,9} plus the spectator and side walls {1,3,4,5,6}. Its five vertices
are the positive-definite deletions of one of these side walls. The
strict-center projection and compact convex-hull argument in the frozen
proof make it an actual compact4-face. A finite vertex of the Coxeter
polytope is simple. Consequently the nine containing facets and the five
side facets exhaust its full wall star: an additional facet meeting S
would meet a vertex already incident to exactly thirteen facets.

The full star splits into the five-node Lorentz core1, with path labels
5,3,3,3, and a finite passive forest of rank9. The producer and independent
checker enumerate the120 permutations of that core and find only the
identity. Thus its vertices and the following edge are intrinsic. The
base edge e=Q intersect S has endpoints obtained by omitting raw1 and raw3
from the simplex sides. In core1 these are vertices4 and0. Choose vertex0
as the canonical endpoint v. This choice avoids counting the other end
of the same edge a second time.

The old prism/base template independently records deletion2, opposite
base facet0, this same core isomorphism and a single base occurrence.
In the wall star of S, {0,9} and {7,8} are two whole A2 components,
orthogonal to the Lorentz core and to each other. The rank-five spectator
may supply further A2 components. The facet q=raw0 contains S and is one
of the two endpoints of a whole A2 component. This is an actual facet,
not an abstract coface label.

## Injection into actual pairs (S,q)

For a fixed actual S, its full wall star fixes the distinguished edge e
and endpoint v. At the simple finite vertex v, dropping the containing
facet q gives one outgoing edge direction outside S. The direction along
e retains q, so the two directions are distinct. Two edge directions at
a simple vertex determine a unique actual2-face. A raw357 quadrilateral
assigned to (S,q) must be this face. Therefore at most one such quad can
be assigned to (S,q).

This argument does not choose a second A2 component or its orientation.
Different such native presentations with the same (S,q) cannot create
different actual quadrilaterals. Conversely, the unique source quad and
trivial source-core automorphism give one canonical assigned S and q
for each source occurrence. The full wall star of a face is intrinsic,
so the argument identifies the same S across vertex presentations.

There are exactly2b candidate facets q among the b whole A2 components
of this fixed normal diagram. This proves the stated upper bound. It
does not claim that every candidate can occur, that the bound is sharp,
or that any formal simultaneous attachment is realizable. It uses no
pairwise unknown as a geometric exclusion.

## Complete family and the signed transfer

Adding two A2 components to the source spectator gives the complete
passive forest of S. This operation is injective on the47 source types.
The producer compares every source context with the frozen complete
prism/base catalogue: each has one base of this type, with the same
core and passive forest. The cap2b belongs to this actual simplex type;
it is not borrowed from a finite-vertex presentation.

For source face Q and represented simplex row S, the exported signed
column has coefficient+1 on Q, coefficient-2b on S, and no vertex-budget
coefficient. Its total contribution over actual faces is nonpositive by
the incidence inequality, so it is a valid transfer. A candidate model
dual with face multipliers alpha must satisfy

    alpha_Q - 2b * alpha_S <= 0.

The simplex aggregate has required payment20, accounting for its ten
triangular2-faces. This demand is unchanged for the sixteen existing
rows. For the31 absent types, a future model extension must construct
their full face rows and associated valid coefficients before installing
the corresponding signed columns. They are not silently assigned zero
mass or treated as current model rows here.

The inherited [shared prism capacity](fresh_prism_simplex_shared_capacity_20261017.md)
uses a single capacity45 for all prisms attached to a4-simplex. Its
raw357 contribution has coefficient2. The new restriction uses the
specific intrinsic base edge and only its whole-A2 outgoing facets.
The certificate compares every represented forced-base capacity45 row
that it touches. This is not a census of every historical signed column.
Instead, strict violation by the separately verified full-model dual
proves that the six violated inequalities are not already enforced by
the old model, even as nonnegative combinations of its constraints.

The six source rows are70404,70406,70408,70414,70420 and70422. Their
positive signed loads sum37155891/10000000000. The largest is face70408:
its new capacity4 gives load18431473/10000000000. These numbers measure
violations of this particular exact dual. They are not a decrease in
an optimized bound, a repaired dual or a target12 certificate.

This mechanism escapes the qualified inherited-vertex-charge barrier:
it has zero vertex cost and changes a relation between face masses.
Increasing vertex allowances cannot repair such a violation while
holding the face masses fixed. A later exact model extension and solve
would still need to establish what new feasible bounds result.

## Verification scope

Final generation2329/7de300 terminated exit0/PASS after normalizing the
saved core-isomorphism pairs to JSON lists. An earlier generation30e984
computed the same mathematical result, but its replay0ed0d3 detected the
tuple/list schema difference. Separate final replay17273/3d73a3 and the
[independent checker](../scripts/check_raw357_simplex_coface_saved_20261018.py)
25100/56b910 also terminated exit0/PASS.

Root and the construction reviewer each completed a full mathematical
source/report/standalone-checker read, including the inherited actual
finite-vertex and compact-simplex argument. Both reviews passed. They
checked the intrinsic edge and canonical endpoint, actual containing
A2 facet, injectivity through the two incident edge directions, complete
47-type correspondence, signed-column direction, unchanged demand20,
and exact current-dual violation scope. Neither reviewer executed the
producer/helper or repeated the native/support/catalogue census or
current full-matrix verification. This report is final and frozen.

The checker imports no producer. It checks the full native core symmetry,
whole-A2 splitting, distinguished roles, exact47-type correspondence,
16 existing columns and31 missing target rows, capacities and current
dual-load arithmetic. The frozen actual-vertex/simplex proof, catalogue
completeness, current full-matrix feasibility and the geometric
(S,q)-injectivity argument remain separately stated inputs. No native
pair scan, clique census or old-matrix replay is repeated.

    python3 scripts/fresh_raw357_simplex_coface_transfer_20261018.py --verify
    python3 scripts/check_raw357_simplex_coface_saved_20261018.py
