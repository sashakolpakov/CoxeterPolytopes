# Finite corners and continuous escape-cutter obstructions

The H20 construction now has an independently reconstructed **145 actual
facets, four proper cusps, and at least 171 finite vertices**, including one
new finite corner. The filtered H9 construction has an independently
reconstructed complete 52-wall cone. Both have exact infinite-volume
witnesses. No new finite-volume polytope or construction dimension has been
found. The arithmetic bound-21 and unrestricted finite-volume bound-163
proof candidates are unchanged; external validation and priority remain open.

## H20: a finite corner and a necessary contact for the next cutter

The [mixed isolated-support cap](fresh_h20_mixed_isolated_cap_20261016.md)
extends the reviewed 141-wall chamber to 142 walls and cuts E2. Three
additional rational residual normals then give the
[145-wall finite-corner cluster](fresh_h20_residual_h4_corner_20261016.md).
The four normals are pairwise orthogonal in the Lorentz five-space
perpendicular to sixteen retained support roots. Together with those roots,
they form a positive-definite rank-20 corner Gram of type
A10+A4+A2+4A1.

The [independent review](fresh_h20_residual_h4_corner_20261016_review.md)
reconstructs the old physical geometry and all four residual normals without
producer arithmetic imports. It checks 570 new pairs beyond wall 140:
71 right and 499 dotted. A separate rational linear solve constructs the
finite corner. Complete finite-neighbor exhaustion proves that it is the
only new finite vertex. All old certified finite vertices and all four
proper cusps persist, and all 145 facets are actual. Generation and fresh
saved replay pass.

The new corner does not close the chamber. An entire rational rectangle
of strict spacelike escapes survives, with exact minimum square
15000061/5000000000>0. The displayed E4 is one point of that rectangle.
The independent review certifies the whole rectangle by affine wall
inequalities and an exact norm minimum. There is no complete old ambient
cone census, so 171 is a lower bound on finite vertices.

The [four orthogonal edge-family obstruction](fresh_h20_h4_edge_escape_gate_20261016.md)
has an [independent standard-library review](fresh_h20_h4_edge_escape_gate_20261016_review.md).
Every new wall right to the sixteen support roots and three of the four
residual caps lies in an explicit Lorentz two-plane. Three entire families
are prevented from cutting E4 by old-wall acuteness. The fourth is prevented
by a normalized Coxeter cyclic invariant strictly between zero and one,
with a separate wall excluding its right-contact endpoint. These are
all-parameter conclusions, not finite angle grids.

A stronger [continuous residual-cutter theorem](fresh_h20_residual_cutter_domain_20261016.md)
now identifies which contact must be dotted. Every acute spacelike E4 cutter
right to the sixteen-root support is strictly dotted to wall 144. Normalize
its residual coordinates as r=C-sum(u_i*n_i), and put
z_i=u_i*sqrt(n_i^2/(-C^2)). Old walls 133 and 140 and the E4 halfspace give
three linear inequalities. Positive rational multipliers

    529/100, 268/25, 759/50

prove z0+z1+z2<=83/100. Hence projection into wall144-perpendicular has
square at most -(3111/10000)*(-C^2)<0. Acuteness then forces the strict
dotted contact. The exact verifier uses independently reconstructed geometry
and has generation and saved replay; the exploratory numerical residual
polytope is not a premise. This theorem has no norm, field, label, height,
or retained-center cutoff. It still permits cutters dotted to wall 144 and
walls outside this residual space. Those are the active search directions.

## H9: independent complete cone and a three-cap obstruction

The [independent filtered52 review](fresh_nonarithmetic_affine_filtered_20261016_review.md)
now reconstructs both the local three-cap formulas and the complete cone.
The final census is 178 finite, zero ideal, and 2,279 exterior extreme rays.
There are 127,764 final ray/inequality checks. All three new facets have
no finite vertex, and the renewed strict escape proves infinite volume.
Generation and fresh saved replays of both reviewers pass.

The [primary cutter-domain gate](fresh_nonarithmetic_escape_cutter_domain_20261016.md)
and [independent review](fresh_nonarithmetic_escape_cutter_domain_20261016_review.md)
clip this cone by r.E>=0. The independent implementation uses rank-eight
common-covector minors, yielding 146 retained and 196 new crossing
generators. Every one of the 342 generators projects strictly future
timelike into each of the three new cap hyperplanes. Convexity extends
this to every nonzero conic combination.

It follows that **every acute spacelike cutter is strictly dotted to all
three new caps**, over all real normals and labels. Thus a single wall
cannot both cut E and meet one of these caps. Multiple additions performing
different tasks remain possible. This explains the four
[bounded bridge pilots](fresh_nonarithmetic_ridge_bridge_20261016.md),
whose zero outcomes alone would not have justified a general exclusion.
The search is continuing with separate ridge closures and compatible
nonpolar cutters. The earlier first52 branch's persistent-cusp exclusion
remains frozen; it has not been confused with this filtered52 branch.

## Compact arithmetic dimension 21: freeing all old bonuses

Dimension 21 over Q(sqrt5) remains open. The
[restricted reallocation diagnostic](fresh_sqrt5_marked_reallocation_proposal_20261016.md)
still fails, but its fixed old coefficients were an artificial restriction.
The [full-free model](fresh_sqrt5_marked_free_bonus_proposal_20261016.md)
therefore makes all 30,169 inherited orbit/context bonuses independent
nonnegative variables, with no old primal coefficient fixed. Its exact
rational dual proves a model lower bound

    425633390406733/12138000036414 = approximately 35.06618793 > 35.

The [H4-count enlargement](fresh_sqrt5_marked_h4_count_free_bonus_20261016.md)
separates compensation by actual H4 count and yields the exact model lower
bound

    78271685370658/2232999997767 = approximately 35.05225502 > 35.

The subsequent [facet-mark/A1-context enlargement](fresh_sqrt5_marked_h4_facet_context_free_bonus_20261016.md)
also has a frozen exact obstruction,
508953802651448/14513000101591, approximately35.06882099. These models were
separated against different finite subsets of vertex inequalities, so their
reported lower bounds must not be compared as fully converged optima; the
larger numerical value does not mean enlarging the model worsened its true
optimum. The newest dual instead highlights the D5 and A1+D4 remainders
beside four H4 components, motivating a separate active refinement.

All three packages have generation, saved replay and literal root-pair recounts
of their dual vertices. They retain complete marked-forest separation and
the old quadrilateral obligations. They impose only selected triangle
families, and newly affected face coverage is not complete. These are
primary exact restricted-model obstructions, not independently reviewed
dimension proofs, fully converged global LP optima, or geometric existence
claims. The numerical lower bounds are also not feasible upper bounds.

The discount caps are inactive in the obstructing numerical relaxations.
The dominant dual vertices and faces identify missing H4 remainder-component
distinctions after the completed facet-mark split. Further enlargements have
separate namespaces and are not included here until frozen. No dimension bound is
promoted from these diagnostics.

## Preservation and scope

The [archive checker](../scripts/audit_finite_corners_continuous_cutters_checkpoint_20261016.py)
and [integrity certificate](finite_corners_and_continuous_cutter_obstructions_2026_10_16_audit.json)
preserve the entire preceding frozen archive and bind these completed
packages. Source hashes and link checks prove integrity, not mathematical
validity. Independent reviews, primary exact certificates, bounded numerical
scouts and active searches have the different scopes stated above. New
candidate walls under active exactification are outside this checkpoint.
