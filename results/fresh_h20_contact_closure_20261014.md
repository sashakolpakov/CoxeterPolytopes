# Exact contact closure and mate replacement in the H20 hybrid lift

The final branch has **76 actual Coxeter facets, eight certified finite edges joining `D13 + D7` vertices to `A12 + E8` vertices, and an exact strict spacelike escape**. It therefore still has infinite volume. The construction changes which facets are retained: it removes the eight large orthogonal mates from the intermediate branch and replaces them by eight caps meeting different old walls. It keeps all 60 walls of the original lifted hybrid.

The finite angles remain among `pi/2`, `pi/3`, and `pi/4`. The unchanged old facets retain the conditional non-quasi-arithmetic obstruction from the [source lift](fresh_dimension_construction_routes_20261014.md). Thus a hypothetical finite-volume completion retaining them would be nonarithmetic; the current infinite-covolume chamber is not being called a nonarithmetic lattice.

## Exact checkpoints

| Branch | Retained/replaced walls | Certified finite geometry | Global volume |
|---|---|---|---|
| Initial contact bridge | Old 60 plus 16 rational caps; replaces the optional 32 caps of the frozen 92-wall branch | 16 specified finite vertices and eight finite edges | Saved 92-wall strict escape survives |
| Final changed-support replacement | Keeps old 60 and eight rational norm-two caps; replaces all eight norm-418 mates by eight norm-two caps over `Q(sqrt(79))` | Eight retained `D13 + D7` vertices, eight new `A12 + E8` vertices, eight finite edges | Previous escape is cut; a new exact strict escape remains |

The first checkpoint is in the [primary checker](../scripts/fresh_h20_contact_closure_20261014.py) and [certificate](fresh_h20_contact_closure_20261014.json). The final checkpoint is in the [replacement checker](../scripts/fresh_h20_contact_closure_20261014_replacement.py) and [certificate](fresh_h20_contact_closure_20261014_replacement.json). Both use exact acceptance predicates. The earlier [92-wall source](fresh_h20_coupled_lift_20261014.md) remains unchanged.

## Why the saved pi/3 event is useful

The old rank-18 support

`S = [0,2,3,4,5,6,7,8,9,10,11,13,14,15,16,19,27,35]`

has a three-dimensional Lorentzian orthogonal complement. In the rational basis with free coordinates numbered `11,19,20`, its Gram matrix and old wall 1 covector are

```text
Q3 = [143/2  -99/2   165]
     [-99/2   71/2  -115]
     [165    -115    381]
a  = [13/2, -9/2, 15].
```

The equations `r^2 = 2` and `(r,old1) = -1` reduce to

```text
16 y^2 - 20 yz + 3 z^2 = 4,
x = (9y - 30z - 2)/13.
```

They have the simple rational solution `r = (-19/2,13/2,6)` in this basis. The vector `s = (-191,99,110)` is orthogonal to it and has norm 418. Both are right to all 18 support walls. The first has a pi/3 contact with old wall 1, while the mate is dotted to that wall. All other required relations are checked exactly.

The anchor's Schur residual after the old 18 walls is `13/22`. Consequently the rank-20 Gram on `S + old1 + r` is positive definite: its last residual is `2 - 22/13 = 4/13`. This is the finite endpoint that the original all-dotted construction lacked. For comparison, prescribing pi/4 or pi/6 there does not give a positive-definite rank-20 extension. The saved pi/6 proposal on the other support type also fails this particular endpoint test.

Transporting this rational pair to the second specified support and taking both independent height sign changes gives 16 new walls. All 1,080 new pairs are checked: **296 right, eight pi/3, and 776 strictly dotted**. Projection of the common strict center certifies every actual facet and both endpoints of each of the eight edges. Their common 19 incident walls have a positive-definite Gram; exact feasibility and strictness of all nonincident inequalities prove the segment is the actual edge.

The separate [local frontier checker](../scripts/fresh_h20_contact_closure_20261014_frontier.py) and [certificate](fresh_h20_contact_closure_20261014_frontier.json) test all 320 outgoing germs at these 16 vertices. They give eight finite edges and **304 open germs**, with 88 exact Schur tests after dotted-pair exclusions. An open germ reaches the ideal boundary without another wall becoming active; this does not deny that ideal-boundary limit. This is a local census, not a census of every component of the chamber.

## The changed-support caps

A bounded search on those new corners found no cutter of the remaining escape compatible with all 76 walls. Allowing the norm-418 mates to be replaced produced useful candidates. The final selection uses the stronger candidate in each of four symmetry-related pairs and then completes the remaining height reflection.

For the first new cap, impose 19 prescribed pairings on

```text
[1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,19,27,34,68]
```

in the intermediate 76-wall numbering: pairing `-1` with wall 34, and zero with the other 18. Impose one further pairing `-1` with old wall 20. Solving the 20 linear equations gives `u + t w`, with rational vectors `u,w` saved explicitly. The norm-two equation is

```text
54 t^2 - 214 t + 35 = 0,
t = (107 + 11 sqrt(79))/54.
```

The other caps use the corresponding transported recipes and sign reflections. The final coefficient field is contained in `Q(sqrt(2),sqrt(79))`; the checker uses the explicit basis `1,sqrt(2),sqrt(79),sqrt(158)` with Fraction coefficients and rigorous integer square-root bounds for signs. It does not accept a numerical cosine fit.

The eight replacements have **572 new pairs: 144 right, 16 pi/3, and 412 strictly dotted**. They retain the eight old `D13 + D7` vertices and create eight `A12 + E8` vertices. The Dynkin types are checked from the literal rational Gram trees. Each new cap has exactly 20 finite neighbors. Exhausting its 19-neighbor subsets proves that it has exactly one finite vertex in this chamber. This last census is only for vertices incident to a replacement cap; it is not a global corner or cusp census.

All 76 facets have exact strict facet points. All three previously bounded cusp links are preserved because every added cap is strictly negative on their null rays. The new caps cut both distinct height-sign images of the prior 92-wall escape (four sign choices; the old escape's second height coordinate is zero).

The global LP is used only to propose a new witness. Rounding its coordinates and moving toward the strict center gives an exactly verified strict future spacelike point with

```text
E^2 = 269344645863312077/40000000000000000 > 0.
```

All 76 inequalities are strict at this point. The final chamber therefore has infinite volume, independently of the incomplete global vertex census.

## A support obstruction and its precise limit

The first exact checker also certifies eight positive-dual identities. On each of the four original support complements, and its second-height reflection, the saved 92-wall escape covector is a positive combination of three old wall covectors. The 18 support walls plus these three witnesses span all 21 ambient coordinates. Hence **every nonzero acute normal perpendicular to that support pairs strictly negatively with the escape**, with no field, norm, label, or parameter restriction.

This claim covers only the eight records explicitly saved in `positive_dual_support_obstructions`. It does **not** cover their first-height reflections: the numerical dual test fails there, and the separate [fixed-curve certificate](fresh_h20_structured_extension_20261014.md) exhibits acute but non-Coxeter cutters. The final replacement caps change the old contact support and are outside the asserted obstruction.

## Numerical search scope and next geometry

The [initial frontier scout](../scripts/fresh_h20_contact_closure_20261014_frontier_scout.py) examines 312 distinct outgoing supports at the 16 specified vertices. It uses 7,760 profiles with right contacts except at most three pi/3 contacts, positive Schur residual, and one additional anchor pairing zero or minus one. It tests 1,534,800 quadratic branches. No cutter survives while all intermediate 76 walls are mandatory.

The [backtracking scout](../scripts/fresh_h20_contact_closure_20261014_backtracking_scout.py) runs the same bounded proposal family retaining only the old 60 walls as mandatory, and returns eight numerically distinct cutters. These counts and exclusions are numerical and bounded; they are not all-label, all-parameter endpoint classifications. The selected replacement normals, their entire pair table, and their actual geometry are subsequently accepted by the exact replacement checker. The [whole-chamber replacement scout](../scripts/fresh_h20_contact_closure_20261014_replacement_scout.py) is likewise only a witness proposer.

The actionable new endpoint is `A12 + E8`: it is reached by changing a contact in an existing finite-edge component, rather than adding another isolated dotted component. The separate [E8-face reduction](fresh_h20_contact_closure_20261014_e8_face.md) turns it into an exact H12 completion problem. A [77-wall changed core](fresh_h20_contact_closure_20261014_h12_replacement.md) completes that face using the standard I12,1 source, and a [79-wall follow-up](fresh_h20_e8_cusp_completion_20261014.md) repairs its two ambient E8 cusp links. These are separately frozen certificates; the global escape survives throughout.

## Replays and review scope

```sh
python scripts/fresh_h20_contact_closure_20261014.py --verify
python scripts/fresh_h20_contact_closure_20261014_frontier.py --verify
python scripts/fresh_h20_contact_closure_20261014_replacement.py --verify
```

The [initial bridge review](fresh_h20_contact_closure_20261014_review.md) combines a fresh primary replay, source/proof review, and a separate small Fraction frame calculation. The [replacement review](fresh_h20_contact_closure_20261014_replacement_review.md) records its own precise scope. No claim of arithmetic independence is inferred from rerunning a primary checker.
