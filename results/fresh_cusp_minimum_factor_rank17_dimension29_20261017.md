# Dimension29 when every affine cusp factor has rank at least17

**Complete internally checked proof.** The geometric ingredients and
their composition have separate internal mathematical reviews. Exact
all-rank scalar identities have a producer certificate and a separately
implemented rational-arithmetic checker. These are internal validations,
not external refereeing or proof-assistant formalization.

Every finite-sided, finite-volume real
hyperbolic Coxeter polytope whose irreducible affine cusp factors all
have rank at least17 has dimension at most29. There is no arithmetic
hypothesis and no bound on the number of factors at a cusp. A rank-r
affine factor has r+1 nodes. Compact polytopes satisfy the cusp condition
vacuously. This is a restricted bound, not an improvement of the
unrestricted163 ceiling or the compact arithmetic21 ceiling.

## Actual cusp factors and the weighting

Suppose, for contradiction, that the dimension n is at least30. The
[factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md)
identifies every connected affine subdiagram with an entire irreducible
factor at its actual cusp. Consequently there is no connected affine
subdiagram of order at most17. All cusp factors are classical A/B/C/D.

Use the original finite Vinberg weights at finite vertices. Their total
is at most3n. At an ideal vertex, give a corner weight1 if the two marked
roots have graph distance at most3 within its varying affine factor, or
belong to a fork-invariant end segment of order at most6; otherwise give
weight0. The marked roots lie in one varying factor by the actual
factor-fiber description, even when the full cusp has several factors.

For a cusp with factor ranks r_i, set s=#factors and V=product(r_i+1).
The [all-rank affine cost](fresh_cusp_spectator_role_payments_20261017.md)
is at most3(r_i+1) per factor. Each node pair has V/(r_i+1) fibers in
the actual cusp link. Thus the total cusp cost is at most3sV. All weights
are nonnegative.

## Every triangle receives at least two

The small-affine exclusion rules out quasi-Lanner minimal indefinite
diagrams; the remaining minimal indefinite diagrams are Lanner of order
at most5. This is the same classification input used explicitly in the
distance-seven factor framework. A triangle star has such a diagram
containing all three black side roots.

Apply the literal cases of Vinberg's triangle proof as audited in the
[finite/simple-ideal transfer](fresh_simple_cusp_vinberg_weight_frontier_20261017.md).
An order-three core or a four-cycle supplies two adjacent marked pairs.
The remaining five-cycle cases supply two payments by the published
finite endpoint/special-edge lemmas; if a contributing corner is ideal,
the required pair is connected in the core with the third black root
deleted, using at most four nodes. Its graph distance is therefore at
most3 in the varying affine factor and its ideal weight is1.

The high-bond case gives a marked pair in a finite H3/H4 component,
with finite weight2. This corner cannot instead be ideal: both marked
roots at an ideal corner belong to its varying affine component, which
is classical and cannot contain a bond5. Extra cusp factors do not move
either marked root into a spectator component. Finite corners keep their
full original finite diagrams, so no finite-payment extension is being
assumed. Every triangle therefore has total weight at least2.

## Every quadrilateral receives at least one

If a quadrilateral has an ideal corner, its varying affine factor has
order at least18. The [large-factor proof](fresh_large_affine_factor_quad_payment_20261017.md)
applies: the other supporting factors are spherical spectators, and
the required small-affine prohibition is available. Its all-rank cycle
and tree arguments, spectator determinant identity, exceptional marked
roles and additional scalar cases give payment at least1.

Now suppose all four corners are finite and unpaid. Use the actual
private/shared decomposition H=U union V union W_shared from the
[finite-quadrilateral audit](fresh_simple_finite_quad_vinberg_transfer_20261017.md).
The four corner diagrams are positive definite, the opposite-side
diagrams are indefinite, and H has negative index one and corank one.
These are actual-face facts and require no simplicity assumption at
other ideal vertices.

Any connected affine subdiagram of H-U contains the opposite v-black
pair; all its other roots are white. It has order at least18, and is
classical. The [cycle-carrier exclusion](fresh_affine_cycle_carrier_exclusion_20261017.md)
and [tree-carrier exclusion](fresh_affine_tree_carrier_exclusion_20261017.md)
exclude respectively its A and B/C/D possibilities. The same holds in
H-V. Thus both entire complements are affine-free. Since U is contained
in H-V and V in H-U, the private diagrams are affine-free as well.
This is stronger than merely saying that U,V themselves are not affine.

The preliminary finite-diagram reduction in Vinberg's proof now has its
required hypotheses. All subsequent affine-sensitive finite carriers
in the audited five cases have at most17 nodes and are excluded by the
global cusp-rank hypothesis. Larger pieces orthogonal to indefinite
subdiagrams are positive definite by the actual-center argument, without
an order restriction. The actual-star corank and the audited nonzero
determinant coefficient remain unchanged. Consequently every case of
the original finite-quadrilateral proof gives a contradiction. These
steps use the published local classification and determinant cases
through the explicit audit, rather than asserting that the ambient
polytope is compact.

Thus every quadrilateral receives at least1. For faces with at least
five sides the required lower bound5-k is nonpositive and follows from
nonnegative weights.

## Global contradiction in every dimension at least30

Let D5=sum_faces(5-number_of_sides), let v_f count finite vertices, and
write kappa(L_v) for the exact cusp-link curvature. The proved payments
give

` D5 <= 3n v_f + 3 sum_ideal_v s_v V_v. `

The [exact shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives the opposite estimate

` D5 >= beta_n v_f + sum_ideal_v kappa(L_v) + 5(n-1), `

where beta_n=n(n-6)/8 for even n and (n-1)(n-5)/8 for odd n.
For n>=30, beta_n>=3n. The
[all-rank curvature comparison](fresh_cusp_link_rank_monotonicity_20261017.md)
proves kappa(L_v)>=3s_v V_v whenever every r_vi>=17. Hence the lower
bound is at least the upper bound plus5(n-1)>0, a contradiction.
Equality of a local coefficient at n=30 causes no difficulty because
the global extremum slack is strictly positive.

This argument has no dimension-enumeration cutoff and does not invoke
the unrestricted163 theorem to cover a tail. Its consequence
is equivalently that any finite-volume example of dimension at least30
must have a cusp factor of rank at most16.

## Verification and scope

The all-rank geometric lemmas and exact scalar identities are the
essential inputs. Finite graph probes corroborate their calculations;
they are not a proof by rank cutoff. The
[composition checker](../scripts/fresh_cusp_minimum_factor_rank17_dimension29_20261017.py)
and its [certificate](fresh_cusp_minimum_factor_rank17_dimension29_20261017.json)
bind the completed geometric packages and verify the all-rank finite
parity inequalities, rank-monotonicity cap, positive increment and
multiple-factor curvature comparison. A
[separate internal review](fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.md)
records the written mathematical audit and independently derives the scalar
identities with rational polynomial arithmetic, without importing the
producer or its symbolic algebra system. Generation and fresh saved replay
are required for both packages and recorded in their certificates.

The tree-carrier proof has root and construction-agent mathematical/source
reviews; the full composition has frontier and construction-agent reviews.
These reads do not constitute a formal verification of the geometric
classification inputs. No unresolved step is identified in this chain.
External refereeing, proof-assistant formalization and publication-priority
claims are absent. The unrestricted bound remains163 and the compact
arithmetic bound remains21.
