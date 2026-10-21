# Coupled projected repairs of the H18-to-H20 lift

The new construction retains the exact 60-wall H20 lift of Vinberg's nonarithmetic H18 hybrid and adds **32 coherent coupled caps**. It has **92 actual facets**, creates **16 actual finite vertices**, and cuts all four height-sign images of the previous strict escape. The three previously completed cusp links stay bounded. A new exact strict spacelike point nevertheless proves **infinite volume**.

This changes the projected wall directions. It is not another assignment of the old roots to transverse axes. Before constructing it, an all-real argument closes a broader raw-height family: arbitrary finite planar height sets over the old short-root directions, provided their convex hulls contain zero and all 26 old long walls remain unchanged. The obstruction below applies to that family only.

Artifacts: [exact producer](../scripts/fresh_h20_coupled_lift_20261014.py), [exact normals and geometry](fresh_h20_coupled_lift_20261014.json), [bounded numerical proposal script](../scripts/fresh_h20_coupled_lift_20261014_scout.py), and [proposal data](fresh_h20_coupled_lift_20261014_scout.json). The frozen [60-wall predecessor](fresh_dimension_construction_routes_20261014.md) supplies the original Gram form, old-old Coxeter gates, signature and three bounded-link proofs. Its source manifest is checked before reconstructing the new data.

The [separate geometry and source review](fresh_h20_coupled_lift_20261014_review.md)
passes. Its Fraction implementation imports no producer arithmetic and
checks all norms, centers, escape signs, finite-vertex inequalities,
support ranks, sign transports and retained cusp signs. The interval
wall-pair gates receive source/proof review and a fresh producer replay;
the review does not claim a second independent all-pair implementation.

```sh
python scripts/fresh_h20_coupled_lift_20261014.py --verify
```

## Why arbitrary height polygons do not fix the raw lift

Keep every old long root unchanged. Over each old short root `r_i`, of square one, permit finitely many normals `r_i+u`, with `u` in the added Euclidean plane. Assume zero lies in the convex hull of that root's heights; then the old projected inequality remains enforced. Require the complete acute Coxeter pair gate.

Each of `L18,L25,L26` has a retained norm-two long neighbor with pairing `-1`. A lift has norm `q=1+|u|^2`, and this unchanged pairing has normalized cosine `1/sqrt(2q)`. Since `q>=1`, the only nonzero finite Coxeter possibilities are `q=1` (label four) and `q=2` (label three). Indeed values with `1<q<2` lie strictly between the consecutive label-three and label-four cosines, and values with `q>2` lie strictly between zero and the label-three cosine. The pair cannot be parallel or dotted.

A zero height cannot coexist with another height over the same old root, because their pairing would be positive. Otherwise every height has length one, and acuteness between two lifts requires their Euclidean pairing to be at most `-1`. Thus the height set is exactly an antipodal unit pair. Consequently each of the three distinguished roots has either a singleton zero height or an antipodal unit pair.

These three old roots are pairwise parallel. For two nonzero height pairs, an angle with `0<|u·v|<1` would create a normalized cross cosine `(1-|u·v|)/2` strictly between zero and one half. Therefore their axes must be equal or perpendicular. In a two-dimensional height space, some pair of the three roots has combined height span of dimension at most one; this is also immediate if a zero-height root occurs.

For that pair put `delta=r_i+r_j`. The certificate checks that it is an old feasible null direction and that every other short root pairs **strictly** negatively with it. Choose a unit height vector `w` perpendicular to the chosen pair's height span. For sufficiently large `T`, the point `T*delta+w` satisfies all lifted inequalities, strictly except at the old incident long walls and the chosen pair. Its square is one. A sufficiently small perturbation toward the old strict timelike center makes every remaining inequality strict while keeping the square positive and the time orientation future. Thus every chamber in this raw-height family has infinite volume, with no bound on the other eleven roots' heights or on their finite angle labels.

The suggested equilateral triple with height square two also fails directly. Its lifted roots have square three. Every short source root has a long pairing `-1` or `-2`, giving squared cosine `1/6` or `2/3`. In either case `2cos(2theta)=4cos^2(theta)-2` would be a rational noninteger, impossible for a finite Coxeter angle. This direct gate is consistent with, but weaker than, the distinguished-triangle obstruction.

These statements allow neither replacement of an old long wall nor additional projected directions. The construction below deliberately leaves this excluded family.

## Four support faces and coupled roots

Work in the predecessor's rational 21-dimensional form of signature `(20,1)`. Replace `sqrt(2)` by one in its saved strict center to obtain a particularly convenient rational point

```text
H=(214,468,728,627,530,438,352,267,186,111,40,-30,-94,
   -153,-208,-262,-312,-356,362,0,0),
H^2=-4296.
```

The producer checks every old inequality and its time orientation against the original center. Let `v` be the predecessor's literal boundary escape of square 25. Its 23 active walls contain four explicitly saved positive-definite rank-18 supports. For each support `S`, project `H` orthogonally to its common three-dimensional Lorentz complement, obtaining `h_S`. Choose the one-dimensional positive vector `w_S` orthogonal to `S`, `v` and `h_S`, using a literal rational nullspace normalization.

For a specified positive rational `t`, define

```text
m = v+t*h_S,
q = m^2 > 0,
d = q/w_S^2 > 0,
r_plus = m+sqrt(d)*w_S,
r_minus = m-sqrt(d)*w_S.
```

The two roots have square `2q`, pair to zero, and are orthogonal to all eighteen roots of `S`. They point toward the same strict center. The four recipes use only two dyadic parameters:

| Supports | `t` | New root norm squared |
|---|---|---|
| 0 and 2 | `75/4096` | `248123425/234881024` |
| 1 and 3 | `107/8192` | `6114849/33554432` |

Their radicands are rational and positive; all coordinates lie in the compositum of `Q(sqrt(2))` with two explicit real quadratic extensions. The artifact stores every coefficient and radicand, rather than rounded vectors.

The original 60-wall chamber is invariant under both transverse-coordinate sign flips. The producer explicitly reconstructs the corresponding old-wall permutations and transports all eight base caps through the four signs. Direct pair checks, rather than an assumption about new-root symmetry, certify the resulting 32 caps together.

## Exact compatibility and geometry

All **2,416** pairs involving a new wall pass: **592 right and 1,824 strictly dotted**. Every right pairing is an exact coefficient identity. Every dotted pairing is proved by outward integer enclosures at scale `2^100`, checking a negative upper bound whose square strictly exceeds the exact norm product. No floating-point predicate accepts any normal, pair, vertex or escape. The numerical scout supplies the dyadic proposals and an initial escape coordinate only.

The common strict center and all nonpositive cross pairings prove that every old and new wall is an actual facet. For each coupled pair there is an explicit finite vertex

```text
p = h_S - (h_S·m/q)*m.
```

Its norm is negative, it is orthogonal exactly to `S` and the two new caps, and every other one of the 92 inequalities is strict. Its incident Gram is the positive-definite rank-18 support block plus two positive orthogonal one-node blocks. The two resulting diagram types are `A11+A6+3A1` and `A10+D7+3A1`; no new `H3`, `H4` or `F4` component occurs at these corners. The four sign transports therefore give **16 distinct actual finite vertices**, with distinct verified incident sets. This is not an exhaustive census of all vertices of the chamber.

Every new cap is strictly negative at all three previously completed cusp rays. Thus their old positive affine relations and quotient rank-19 link proofs remain intact; the extra walls are absent sufficiently high in each corresponding horoball. The construction does not assert a complete inherited or new cusp census.

All eight corresponding caps cut each of the four sign images of the previous strict escape. This is checked by exact positive pairing bounds, not inferred from proximity to the original boundary point.

## Global volume and the next necessary change

The full 92-wall chamber has an exact rational strict future spacelike escape of square

```text
1411253767129583933/200000000000000000 > 0.
```

Its literal 21 coordinates, negative center pairing and all strict inequalities are stored and replayed. Thus the chamber still has infinite volume.

There is also a precise local limitation. Each new cap has **exactly nineteen finite neighbors**: its eighteen old support walls and its paired cap. All other walls are strictly dotted to it. At each of the sixteen new simple finite vertices, every one of the twenty outgoing edge germs retains at least one new cap. Every wall outside that starting vertex is dotted to a retained cap, so it cannot meet that edge at any finite or ideal point. The edge therefore reaches the light cone without another wall becoming active. It has an ideal-boundary limit, but no second finite vertex or ideal vertex with bounded link. This proves that **all 320 specified outgoing germs remain open**; it is not a claim about every edge of the whole chamber.

A future completion must introduce walls meeting these new caps and changing their contact supports, or replace some of them. Further banks whose only cross relations to these components are dotted cannot close these edges. This supplies a concrete frontier for subsequent work, rather than treating the extra facet count as volume progress.

## Arithmetic gate and scope

The retained norm-two pair `L1,R10` still gives the Cartan cycle `12+8sqrt(2)`, while the retained rational pair `R23,R24` has Gram determinant `-32`. The [predecessor's conditional non-quasi-arithmetic argument](fresh_dimension_construction_routes_20261014.md#arithmetic-gate-before-further-work) therefore applies unchanged: any finite-volume completion retaining these facets would be non-quasi-arithmetic, hence nonarithmetic. The present infinite-covolume group is not being called a nonarithmetic lattice. No prohibited noncompact arithmetic H20 lattice search is involved.

The exact construction adds only right and dotted relations; it introduces no noncrystallographic finite ridge angle. A separate [bounded frontier scout](../scripts/fresh_h20_coupled_lift_20261014_frontier_scout.py), with [saved numerical proposals](fresh_h20_coupled_lift_20261014_frontier_scout.json), tested labels two through thirty plus a parallel event on the four support faces. Its 6,822 candidate events yielded 112 numerically passing events (not a distinct-cap census), including an additional old label-six contact. It uses the original irrational source center, unlike the exact rational-center construction above. None of these finite-contact proposals has been accepted exactly or used to claim closure. The checkpoint stops at this exact 92-wall chamber and its explicit open frontier.
