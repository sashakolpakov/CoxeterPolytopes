# Dimension29 when every irreducible cusp factor has rank at least9

**Complete internally checked proof.** Every finite-sided,
finite-volume real hyperbolic Coxeter polytope whose irreducible affine
cusp factors all have rank at least9 has dimension at most29. There is
no arithmetic hypothesis and no bound on the number of factors at a
cusp. Compact polytopes satisfy the condition vacuously.

Equivalently, a finite-volume example in dimension at least30 must have
an irreducible cusp factor of rank at most8. This strengthens the
[previous rank17 theorem](fresh_cusp_minimum_factor_rank17_dimension29_20261017.md).
It is a restricted dimension bound; it does not itself change the
unrestricted163 ceiling or the compact-arithmetic21 result. No new
example or publication-priority claim is made.

## The cusp hypothesis and weights

Assume for contradiction that n>=30. The
[actual factor-fiber lemma](fresh_cusp_factor_payment_framework_20261017.md)
places every connected affine subdiagram in its actual ideal vertex and
identifies it with a full irreducible affine factor there. A rank-r
factor has r+1 nodes. Therefore the hypothesis excludes every connected
affine subdiagram of order at most9, and all remaining affine factors
are classical A/B/C/D of order at least10.

At finite vertices use the original finite Vinberg weights, of total
cost at most3n. At an ideal corner assign weight1 if its two marked
roots have graph distance at most3 in their varying factor, or belong
to its fork-invariant end segment of order at most6; otherwise assign0.
The [all-rank affine cost](fresh_cusp_spectator_role_payments_20261017.md)
is at most3m for a factor of order m.

If the cusp-link factor ranks are r_1,...,r_s and V=product(r_i+1), each
marked pair from factor i occurs in exactly V/(r_i+1) link fibers.
Consequently total cost at that ideal vertex is at most3sV. This
count does not assume that the ideal vertex is simple. Each face
corner's marked roots belong to a single varying affine factor; other
cusp factors appear as spherical spectators in its face star.

## Triangle and quadrilateral payments

The order-nine prohibition makes every minimal indefinite subdiagram
Lanner: a quasi-Lanner diagram has order at most10 and contains a
proper affine subdiagram of order at most9. Hence a triangle star has
a Lanner core of order at most5 containing all three black side roots.

Apply the published triangle cases through the
[finite-weight transfer interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md).
The three- and four-node cases supply two adjacent marked pairs. In
the five-cycle cases, the finite endpoint/special-edge arguments supply
two payments. If a contributing corner is ideal, its two marks are
connected after deleting the third black root by a path on at most
four nodes. They lie in the varying affine factor and are paid by
distance at most3. The high-bond alternative contributes finite
weight2 from an H3/H4 marked component. It cannot be an ideal marked
component: the varying factor is classical and has no bond5.
Thus every triangle receives at least2. This is the same triangle
transfer as in the previous theorem, with exactly the same small-affine
prohibition; no new ideal-port transitivity is assumed.

If a quadrilateral has an ideal corner, its varying factor has order
at least10. The [new complete ideal-corner proof](fresh_affine_order10_quad_payment_20261017.md)
gives payment at least1 with arbitrary spherical spectators. The
additional shorter-order case is the affine-D11 configuration in which
both replacements have E8 roles. Its low opposite contacts force
affine E8/A7 diagrams of order9 or8; after excluding them, its exact
normalized determinant is strictly positive. The remaining scalar,
finite-role and classical extension arguments are uniform in rank.

For a quadrilateral with four finite vertices, suppose all four
payments vanish. Its actual star has negative index one and corank one,
and its four corner diagrams are positive definite. Use the
private/shared decomposition H=U union V union W_shared. The
[new carrier proof](fresh_affine_carriers_order10_20261017.md)
excludes every affine A/B/C/D subdiagram of H−U or H−V already from
order10 onward. Its new small-order ingredient is the affine-E7 leaf
witness with arms(1,3,3); the grounded-resistance and actual shared-path
arguments then give the same complete private/carrier exclusion.

Both complements, and therefore U,V, are affine-free. This supplies
the preliminary hypotheses of the published finite-quadrilateral
casework. The [refined source audit](fresh_vinberg_finite_quad_affine_cutoff9_20261017.md)
shows that its remaining affine exclusions need only orders through9.
The former17 bound enclosed entire configurations rather than minimal
affine obstructions. The published finite case tables and determinant
calculations remain explicit inputs. Larger orthogonal pieces are
positive definite by the actual-center argument, and the actual-star
determinant coefficient is nonzero as in the frozen transfer.
Every published case is contradictory, so this quadrilateral too
receives at least1.

All weights are nonnegative. For polygons with at least five sides
the lower bound5−number_of_sides is nonpositive and is automatic.

## Exact global curvature contradiction

Let D5=sum_faces(5−number_of_sides), let v_f be the number of finite
vertices, and let s_v,V_v describe each cusp link. The proved payments
give

` D5 <= 3n v_f + 3 sum_ideal_v s_v V_v. `

The [exact shape inequality](fresh_cusp_link_shape_curvature_20261017.md)
gives

` D5 >= beta_n v_f + sum_ideal_v kappa(L_v) + 5(n−1), `

where beta_n=n(n−6)/8 for even n and (n−1)(n−5)/8 for odd n.
For n>=30, beta_n>=3n. The
[new all-rank curvature comparison](fresh_cusp_rank7_fine_curvature_20261017.md)
proves kappa(L_v)>=3s_v V_v whenever every factor rank is at least7
and the total rank n−1 is at least29. The present rank9 hypothesis
satisfies those assumptions.

For completeness, that scalar comparison uses a positive1/18 increment
under increasing a rank already at least7, a variance estimate for
arbitrarily many factors, and21 exact minimal profiles when there are
two through six factors. This is a proved monotonicity reduction, not
a finite maximum-rank census. Its only equality case is one factor of
rank29. The positive global term5(n−1) still makes the lower bound
strictly exceed the upper bound, including n=30. This contradiction
proves n<=29.

The argument has no dimension-enumeration cutoff and uses neither the
unrestricted163 theorem nor an optimization to cover a tail.

## Verification and scope

The [composition checker](../scripts/fresh_cusp_minimum_factor_rank9_dimension29_20261017.py)
and [certificate](fresh_cusp_minimum_factor_rank9_dimension29_20261017.json)
bind24 source/data/report inputs, rederive the all-rank finite parity
coefficients using rational polynomial arithmetic, check the affine
cost identities, and recompute the21 reduced curvature values.
It does not replay the old large rank/label enumerations. The new
computational geometric packages have successful generation/fresh
replay; all new geometric proofs and the report-only source audit have
separate root mathematical/source reads with their scope recorded.

The finite-quadrilateral source audit includes a fresh reading of the
primary PDF text but not a fresh visual reconstruction of its tables;
the original published casework remains a proof dependency. These
internal checks are not external refereeing or proof-assistant
formalization. Generation, fresh saved replay, syntax, all24 direct
bindings and local links pass. A separate root mathematical/source
review of the complete composition passes the actual affine-factor
equivalence, triangle roles, ideal and finite quadrilateral transfers,
weight costs, rank7 curvature and positive global slack. That review
does not duplicate the published tables or producer arithmetic. No
unresolved step is identified in this proof chain. Source, certificate
and report are final and frozen.
