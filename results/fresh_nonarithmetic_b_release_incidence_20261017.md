# A three-cap incidence restriction after releasing52

Releasing52 is the cheapest of the three examined changes in [compact-face deletion cost](fresh_nonarithmetic_b_release_costs_20261017.md). Nevertheless, a repair by at most three new actual S6 facets, each meeting at most two retained S6 facets along geodesic edges, cannot improve on the control's15 open ends **if the resulting S6 chamber has no ideal vertices**. The statement is independent of labels, coefficients and the seven possible A2 contact profiles.

The hypotheses retain every old wall except52, require all normals mutually acute and a common strict future timelike interior point, and count only actual facets of the S6 face. The no-ideal hypothesis is essential to interpreting the count as unclosed edge ends. This is not an unrestricted finite-volume obstruction.

The frozen complete boundary has45 geodesic edges and25 finite vertices. Exactly six geodesic edges and five finite vertices involve52. After its deletion the remaining counts are therefore E=39 and V_f=20, so I=2E−3V_f=18. Exterior-only edges of the auxiliary projective cone are excluded from E.

These old incidences persist under the stipulated additions. For any positive-definite old support, write the projection of a common interior point X as X−Σc_i n_i, where c=G⁻¹(⟨X,n_i⟩). The Stieltjes inverse is nonnegative, hence c_i≤0. For an acute new normal r,

`⟨r,X−Σc_i n_i⟩ = ⟨r,X⟩−Σc_i⟨r,n_i⟩ ≤ ⟨r,X⟩ < 0`.

For a finite-vertex support this proves strict retention of that vertex. For the white support plus a retained positive-definite two-side block, it proves the old geodesic edge remains actual. No whole old edge or finite vertex is silently removed. Old-only positive-definite triples were already actual in the acute control, so deleting52 introduces no additional old-only finite vertices.

Let a count old–new geodesic edges, b new–new geodesic edges, and x,y,z new finite vertices incident to respectively one,two,three new facets. Finite vertices are simple: otherwise a dependent positive-semidefinite acute normal block would have a nonnegative kernel vector, contradicting the common strict interior point. Thus

`ΔI = 2(a+b)−3(x+y+z)`.

At most three new facets with at most two old neighbors give a≤6 and b≤3. Each cap supports at most one finite vertex using two old neighbors, so x≤⌊a/2⌋. Counting finite endpoints of old–new edges gives x+y≤a. Counting those of new–new edges gives y+3z≤2b. Finally z≤1, and z>0 requires b=3.

The [small exact checker](../scripts/fresh_nonarithmetic_b_release_incidence_20261017.py) exhausts all238 integer tuples satisfying these necessary conditions. It verifies ΔI≥−3, with equality only at `(a,b,x,y,z)=(6,3,3,3,1)`. Hence I remains at least15. The extremal tuple is a relaxed incidence possibility, not a constructed configuration.

A separate analytic check gives the same minimum. If z=1, then b=3 and x+y+z≤a+1, so ΔI≥3−a≥−3. Equality forces a=6 and x+y=6; the bounds x≤3 and y≤3 then force x=y=3. If z=0 and b=0, then y=0 and ΔI≥0. For b=1, use x+y≤min(a,⌊a/2⌋+2): when a≤4, ΔI≥2−a≥−2; when a≥4, ΔI≥a/2−4≥−2. For b=2, ΔI≥4−a≥−2; for b=3, ΔI≥6−a≥0. This proves the minimum and its unique equality tuple without relying on the enumeration.

When ideal vertices are present, I counts their incident ends too. The number of genuinely unclosed ends is `I−Σ deg(v_ideal)`. Their valences may exceed3. The checker therefore explicitly assumes no ideals for its open-end conclusion; it does not discard that correction or invoke a cusp theorem to extend its scope.

A sharper subcase explains the continuous scout. Two new walls with exactly the old finite contacts `(8,12)` and `(12,47)`, a finite mutual contact, all other old native contacts dotted, and all three designated endpoint triples positive definite have `(a,b,x,y,z)=(4,1,2,1,0)`. They close the12-face but increase I by1, giving19. This conclusion does not depend on the sampled conic parameters.

Thus a successful no-ideal release52 repair with at most three caps needs a cap with at least three retained intrinsic geodesic neighbors. Alternatively at least four new actual S6 facets are needed. Neither option is excluded here.

Generation and fresh saved replay passed for the [exact result](fresh_nonarithmetic_b_release_incidence_20261017.json). The parent agent separately checked the geometric retention argument and supplied the analytic minimization above; the238-case arithmetic is one implementation. No new wall or finite-volume example is claimed.
