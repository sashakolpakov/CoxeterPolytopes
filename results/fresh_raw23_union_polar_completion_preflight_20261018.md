# Raw23 unions: a complete extreme-ray polar-method preflight

The three retained raw23 unions cannot be made compact by the tested sequential polar procedure. This is a **restricted method obstruction**, not an exclusion of arbitrary acute completions or of the original raw23 pairs. No face inequality, vertex capacity, model column, or arithmetic bound changes.

Each original twelve-wall cone has50 projective extreme rays:18 timelike and32 spacelike, with no null extreme ray. In the mixed-kind case, none of the32 spacelike polar directions passes the necessary native tests. In either same-kind case, exactly one does. Adding that candidate cuts every original spacelike extreme ray, but the new cone has52 timelike and30 spacelike extreme rays. None of those30 polar directions passes the necessary native tests. An exact780-candidate scalar check also excludes every nonzero passive lift of these156 extreme projected directions. Cutting the original outside rays is therefore insufficient for compactness.

The [producer](../scripts/fresh_raw23_union_polar_completion_preflight_20261018.py), [saved certificate](fresh_raw23_union_polar_completion_preflight_20261018.json), and [standalone checker](../scripts/check_raw23_union_polar_saved_20261018.py) are new. They bind the frozen [complete raw23 relative and union packet](fresh_raw23_same_family_relative_pilot_20261018.md). Its235 source classes,175 relative cases, three retained unknowns and actual supporting-wall interpretation remain unchanged.

## Precisely restricted candidate method

The known twelve roots span a ten-dimensional Lorentz space. Their four passive directions in the ambient rank13 finite frame are orthogonal and positive definite. This preflight requires the Lorentz projection of every proposed additional normal to be a positive multiple of a current spacelike projective extreme ray of the known-wall cone. Arbitrary passive components are allowed in the necessary-condition relaxation and then excluded by the exact scalar test below.

Thus it tests a finite set of projected directions, not every possible additional facet. A positive combination of cone rays can define a different projected direction; such non-extreme projections, with or without passive components, remain open. A normal with zero Lorentz projection is also outside this method, and cannot cut any of the displayed rays in the Lorentz span.

A passing candidate is only a normal satisfying the necessary old-wall native pair tests and preserving the currently finite vertices. Higher-cycle arithmetic, compatibility with unknown future walls, and actual occurrence in the original polytope are not asserted. Treating these necessary-test survivors as eligible makes the tested method more permissive; its failure remains a statement only about this finite polar procedure.

## Complete cone enumeration

Choose ten independent known normals as a Lorentz basis. In its exact Q(√5) coordinates, let G be the basis metric and let the rows of M be the products of each known wall with the basis. The outward-normal convention gives the projective cone

    C = {x : Mx <= 0}.

The normal rows span the ten-dimensional space, so C contains no line. The producer exhibits a strict timelike interior point by summing the normalized finite rays and checking every wall inequality strictly. It also checks that all finite rays belong to one time orientation. Hence C is a full-dimensional pointed cone with an interior point in the relevant hyperbolic sheet.

Every extreme ray of a pointed ten-dimensional polyhedral cone has at least nine linearly independent active wall inequalities. The producer therefore enumerates every nine-wall active set. It computes its exact kernel, rejects rank-deficient sets and kernels whose two orientations both violate a wall inequality, and deduplicates the feasible oriented rays by their complete product vectors. Normalization makes the first nonzero wall product−1. No numerical tolerance is used.

For each initial twelve-wall cone there are220 active sets. All50 rays are found, with signs18 negative and32 positive. The finite rays have exactly the18 supports from the frozen union packet. The original matrices have physical inertia(9,1,2); the ten-dimensional basis removes the two redundant wall relations without changing the cone.

A spacelike ray in such a cone witnesses noncompactness of the known-halfspace intersection: the segment from the strict timelike interior point toward that ray crosses the future light cone while remaining in the cone. This explains the need for further cutting walls. It does not imply that such walls cannot exist.

## Exact polar filters

Let x be a spacelike ray, with squared norm N=x·x>0. Its candidate physical outward normal is sqrt(2/N)*x. If an old root has squared scale s, its physical norm2 cross-product square is

    2*(x·old_scaled_root)^2 / (N*s).

All old-wall products are nonpositive because x lies in C. The candidate must also have strictly negative product with every finite ray. A positive product would remove a known finite vertex, and zero would add a new facet through a vertex already supplied with its full simple set of facets.

For a finite Coxeter angle over Q(√5), the normalized product square below4 must be one of

    0, 1, 2, (3+sqrt5)/2, 3, (5+sqrt5)/2,

corresponding to labels2,3,4,5,6,10. The finite list follows directly from the fixed field: for a finite label m≥3, the square is4cos²(π/m)=2+2cos(2π/m). Membership in Q(√5) forces the real cyclotomic degree φ(m)/2≤2, leaving m=3,4,5,6,8,10,12. The values for8 and12 involve√2 and√3 respectively, so are excluded. Label2 gives zero. The ground-field and arithmetic hypotheses themselves remain inherited from the frozen model.

Every squared product must also be an algebraic integer in Q(√5): its trace and norm must be integers. Its nontrivial conjugate must lie in[0,4], as required by the conjugate two-root positive-semidefinite Gram matrix. Values at least4 are retained when these necessary tests pass; the parallel boundary value4 is not rejected merely for being parallel. No higher-cycle sufficiency is claimed.

All32 initial rays in every case preserve the18 finite vertices. The mixed-kind case has no candidate passing the finite-angle test. Each same-kind case has one candidate; its three nonzero physical contacts have squares1,1,(3+√5)/2. The two descriptions correspond under the frozen shared-edge endpoint interchange. The other31 poles in each same-kind case fail a native finite-angle condition.

## Arbitrary passive lifts along each extreme projection

Write a proposed full physical norm2 normal as

    n = sqrt(t) * n0 + z,    n0 = sqrt(2/N) * x,

where z is orthogonal to the known Lorentz span. Its passive metric is positive definite, so t>0 and the norm equation gives t≤1 and z²=2(1−t). If s_j is a zero-passive physical contact square computed above, the corresponding full contact square is exactly t*s_j: no unknown passive contact has been set to zero.

For every tested direction the certificate identifies a known wall with0<s_j<4. Then0<t*s_j<4. This contact can be neither right, parallel nor dotted, and must equal one of the five positive finite squares. Consequently t belongs to the exact five-element candidate list

    {1/s_j, 2/s_j, ((3+sqrt5)/2)/s_j, 3/s_j, ((5+sqrt5)/2)/s_j}.

Each candidate first passes the physical bound0<t≤1. Every known-wall scaled contact is then tested for a finite Coxeter angle when below4, algebraic integrality, and the conjugate two-root bound. The conjugate known-span Gram matrix is also checked positive definite. Since x has coordinates in Q(√5), the conjugate squared norm of its normalized projection is2*t′. A real conjugate extension requires0<t′≤1: positivity makes the projection real, and the bordered-Gram Schur complement is2(1−t′)≥0. This extra known-span condition is checked independently. Products with unknown passive walls, their higher cycles and their further conjugate restrictions are left unrestricted; omitting them enlarges the candidate set and cannot invalidate an exclusion.

Among all780 scalar candidates, the only survivors are t=1 for the same two initial same-kind poles already found without passive components. The norm equation then forces z=0. No subsequent extreme projected direction admits a lift. Thus allowing arbitrary passive components along these extreme projections does not rescue the sequential procedure. This does not test a non-extreme projected direction.

## Why one apparently successful cut still fails

The sole passing same-kind normal has positive product with all32 original spacelike rays. Its halfspace removes every one of those directions. It also has strictly negative product with all18 old finite rays, so none of those vertices is clipped.

Nevertheless, the new thirteen-wall cone has new extreme rays on previously intersected faces. The producer enumerates all715 nine-wall active sets for each same-kind case, rather than inferring compactness from the original rays. The exact result is82 rays:52 timelike and30 spacelike, with no null ray. Thus the new known-wall intersection remains noncompact.

All30 new spacelike directions are again tested as polar normals against all13 current walls and all52 finite rays. Every one fails either a finite-angle or an algebraic-integer two-cycle requirement. There is no second eligible extreme-ray polar truncation. Since the first same-kind candidate was unique and the mixed-kind case had none, this exhausts the stated sequential method.

Across the three initial cases and two one-pole successors, the packet checks2,090 active sets,156 spacelike projected directions and780 passive-lift scalar candidates. It introduces no assumption that the32 or30 outside rays are actual vertices of a compact Coxeter polytope; they are projective rays of the explicitly constructed partial halfspace system.

## Verification and unresolved geometry

The final passive-lift extension passed generation24674/f98434, separate producer replay53717/8f7c0d and independent saved checker59090/dfca87, all exit0/PASS. Generation also compared every cone/ray record and original passing-pole list with the initial zero-passive version and found them unchanged. The earlier version's generation22785/ff9e03, replay9408/08c7f7 and checker68063/aae946 are historical runs only. Metadata preflighte2b13c checked all6 direct source bindings, both Python ASTs and all4 local links.

Root completed a full mathematical source/report/standalone-helper read (fe5d03/68c5d0/f84d48): PASS. Its scope covered the pointed rank10 cone and complete rank9 active-set domain, normalization/contact scaling, finite-label list, five-scalar passive-lift completeness, conjugate known-span Schur condition, preserved-vertex restriction and successor recomputation. The distinct producer/helper eliminations were checked by reading. Root did not execute the scripts or independently reprove the frozen actual-vertex, catalogue or field hypotheses. The report read preceded only terminal-run bookkeeping.

Construction independently completed the full producer/report/standalone-helper mathematical read, including the final passive-lift extension: PASS. The read covered the same complete active-set and five-scale domains, exact contact and conjugate conditions, finite-vertex retention and both successor cones, as well as the helper's separate field arithmetic and reverse-column elimination. It did not execute the scripts, replay historical cones or reprove the original union, catalogue, field or actual-vertex inputs. Both reviews preserve the restriction to sequential extreme projected directions. Source, certificate, report and helper are final and frozen after these completed scopes.

The standalone checker imports no producer. It uses its own rational Q(√5) arithmetic and a reverse-column kernel elimination, reconstructs every active-set decision and every oriented ray, verifies strict interior points and signs, and rechecks all normalized contact tests, the five-scale passive-lift domains, the conjugate known-span bounds and the unique-pole successors. The original native union geometry, actual-vertex theorem, source catalogue and fixed-field arithmetic setting remain frozen inputs. This is not a historical catalogue replay, an optimization, or a full chamber verification.

The result narrows the next search: a successful completion must depart from this restricted sequence, for example by using a normal whose Lorentz projection is a non-extreme direction of the acute-extension cone. Proving a genuine pair exclusion would require controlling such directions and their passive lifts. No such conclusion, no new occurrence quota, and no improvement to the surviving full-model witness follows from this packet.
