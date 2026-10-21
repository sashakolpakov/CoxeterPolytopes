# Cusp-link curvature, arithmetic field gates and boundary replacements

This checkpoint adds an internally checked **dimension bound85** for
finite-volume hyperbolic Coxeter polytopes whose every irreducible affine
cusp factor has rank at least9. It allows any number of factors per cusp
and imposes no arithmetic hypothesis. The unrestricted163, all-simple29
and compact-arithmetic21 ceilings remain unchanged. No new finite-volume
polytope or construction dimension has been established.

## Exact cusp-link shape and the dimension85 theorem

The [new counting inequality](fresh_cusp_link_shape_curvature_20261017.md)
retains the full product of simplices at each cusp. If d=n−1,
H(t)=product(1+t+...+t^r_i)=sum h_j t^j, V=H(1), and
G=sum h_j min(j,d−j), the local curvature is kappa=(5/2)G−dV/2.
The original polytope satisfies

` D5 >= beta_n v_f + sum_ideal_v kappa(L_v) + 5(n−1). `

This inequality has no rank or factor-count restriction. Its
[independent internal proof review](fresh_cusp_link_shape_curvature_20261017_root_review.md)
checks the actual vertex-figure edge counts and the quantitative slack
from the global minimum and maximum of a generic linear function.

Under the minimum-rank9 hypothesis, the
[distance-seven payment framework](fresh_cusp_factor_payment_framework_20261017.md)
gives finite-vertex cost at most7n and cusp cost at most
V sum_i min(7,r_i/2). The [dimension85 proof](fresh_cusp_minimum_factor_rank9_dimension85_20261017.md)
checks all50,330,195 rank partitions in ambient dimensions86 through163
by exact integer convolution, with a strict margin greater than1/5 per
link vertex. The earlier unrestricted163 proof covers every larger
dimension. A separate [review and partition counter](fresh_cusp_minimum_factor_rank9_dimension85_20261017_review.md)
checks enumeration completeness, integer bounds and all154 per-rank
counts through a different recurrence. Full generation and saved replay
both pass. This is a restricted theorem; a cusp with a factor of rank
at most8 remains outside it.

The framework also proves a positive-semidefinite spectator identity:
adding elliptic supporting roots orthogonal to an old affine factor
can only increase the relevant normalized four-root determinant. It
transfers strict scalar contradictions. It does not automatically
transfer every finer corner payment to nonsimple cusps. The new
[role-specific E7/D8 lemmas](fresh_cusp_spectator_role_payments_20261017.md)
close two exceptional spectator extensions and bound the proposed affine
weight by3 times factor order, but full face payment remains open.

## Fixed-field arithmetic constraints

The [two-spectator cap family](fresh_sqrt5_multispectator_cap_capacity_20261017.md)
adds555 actual cap columns with complete colored-core orientation and
support-intersection accounting. Its
[independent pair/capacity review](fresh_sqrt5_multispectator_cap_capacity_20261017_review.md)
checks exact signs, algebraic integrality and the combinatorial capacities.
The [three-spectator family](fresh_sqrt5_triple_a1_cap_capacity_20261017.md)
adds93 further actual columns.

For the fixed ground field K=Q(sqrt(5)), Vinberg cyclic products must both
be integral and belong to K. The
[fresh ground-field tightening](fresh_sqrt5_multispectator_groundfield_cap_capacity_20261017.md)
lowers187 of the555 capacities. The
[independent79-test review](fresh_sqrt5_spectator_groundfield_gate_20261017.md)
identifies four integral degree-four squared products outside K, two in
each spectator family. No exclusion over a larger ground field follows.
The earlier integrality-only bounds were valid weaker constraints.

The [combined local repair](fresh_sqrt5_a3_d5_cap_local_cover_update_20261017.md)
and [complete joint scan](fresh_sqrt5_a3_d5_cap_joint_ports_20261017.md)
give a common exact witness with ratio **12.12870037679824>12**, improved
from12.25714671504 at the preceding checkpoint. They retain2,188 simplex
rows and336 actual quadrilateral provenances on648 forests. All3,712
literal packing terms,2,056 cap checks (including preserved weaker
columns),94,722 inherited signed columns and2,711,432 joint inequalities
pass. The1,501 selected cap columns comprise851 prior columns,555
strengthened columns,93 triple-contact columns and two further A3/D5
columns. No optimizer ran.

The two [A3/D5 cap witnesses](fresh_sqrt5_a3_d5_cap_witnesses_20261017.md)
are Q61927<=2S81801 and Q75841<=S78789. They invalidate the preceding
12.1296605 masses; the final local repair and joint replay include them.
The current obstruction remains above the target12. These are
partial-model results; even a successful weighting at12 would require
the remaining face, field and forest coverage before an ambient
compact-arithmetic theorem. The next geometric target is compatibility
of adjacent cofaces in the dominant raw-core23 five-prism family.

## Construction tests and current exact control

The [all-label released-wall49 gate](fresh_nonarithmetic_release49_endpoint_gate_20261017.md)
excludes both alternating three-pentagon disks while all other63 old
walls are retained. Its896 recipes reduce to18 finite-corner recipes
and36 real branches, each with a strict obtuse blocker. This is an exact
local obstruction, not a no-completion theorem.

The [continuous ring deformation](fresh_nonarithmetic_release8_ring_deformation_20261017.md)
reaches the missing right angle numerically, but the two valid wall
selections have18 or15 open ends after actual deletion costs. A
[narrow source review](fresh_nonarithmetic_release8_ring_deformation_20261017_review.md)
checks that the invalid eleven-end selection is not promoted. This
remains a numerical diagnostic, not an exact new chamber.

The [joint replacement of49 and62](fresh_nonarithmetic_release49_62_replacement_20261017.md)
now supplies the new exact construction control: **63 actual facets,
199 complete ambient finite vertices and10 S6 open ends**, improved from
64/198/15. The new wall is right to S6,7,58, meets33 at pi/4 and has
strictly dotted contacts with all other retained walls. Exactly one old
finite vertex is lost and two are gained; every old compact face of
dimension at least2 persists. A [separate root review](fresh_nonarithmetic_release49_62_replacement_20261017_root_review.md)
checks the actual-facet and complete finite-census arguments.

The [complete intrinsic boundary certificate](fresh_nonarithmetic_release49_62_boundary_20261017.md)
checks all816 triples and all48 projective edges. It gives26 finite and
six exterior rays, with34 compact edges, ten one-open edges and four
edges outside the hyperbolic ball. No ideal or nonsimple ray is omitted.
The old compact S6 side faces are preserved.

The new wall cuts both the previous control's saved escape and the
older escape that62 was introduced to cut. A renewed strict spacelike
escape is certified separately, so the changed chamber still has
**infinite volume**. This is progress in the exact construction control,
not a new finite-volume example or a complete ambient exterior-cone census.

## Archive and verification

The [archive checker](../scripts/audit_spectator_field_gates_and_boundary_replacements_20261017.py)
and [saved audit](spectator_field_gates_and_boundary_replacements_20261017_audit.json)
preserve the [previous checkpoint](simple_cusp_dimension29_and_central_truncation_20261017.md)
and all its frozen predecessor hashes. In addition to source and link
integrity, the new checker independently composes the saved cap and
packing certificates using Fraction arithmetic. Its scope does not
include formal verification of the written geometry or rerunning every
historical calculation. External mathematician refereeing, Lean
formalization and publication-priority claims remain absent.
