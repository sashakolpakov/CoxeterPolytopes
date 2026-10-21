# A/P-circle diagnosis of the fourth-cap angle mismatch

The continuous necessary domain admits finite K–R contacts. It does not
explain the preceding failure by a uniform dotted-separation inequality.
Instead, the1,016 returned finite K–R angle events all give a non-Coxeter
K–L product. In263 events K–L is the only native blocker. No chamber is
promoted, and the exact63-wall/199-finite-vertex/ten-open-end infinite-volume
control remains unchanged.

The [source](../scripts/fresh_cu_release47_ap_circle_20261017.py) and
[saved result](fresh_cu_release47_ap_circle_20261017.json) diagnose the same
83 selected event frames and516 saved R caps as the frozen
[fourth-cap subdivision](fresh_cu_release47_interior_cap_20261017.md).
There is no new event grid. K–L is allowed to vary continuously; symbolic
high K–A/K–P bonds are reduced to finite frame equations. These are analytic
domain reductions evaluated numerically, not an interval-certified
all-label obstruction for the frozen numerical frames.

## Whole-circle necessary domain

Use intrinsic normals in the four-dimensional Lorentz complement of S6.
The span of A and P is Lorentzian at every selected frame. If
\(G=\operatorname{Gram}(A,P)\), and the prescribed native K contacts are
\(-\theta_{m_A},-\theta_{m_P}\), let \(r\) be these two products after
subtracting the A2-profile offsets. The projection and residual radius are

\[
 k_0=r^T G^{-1}(A,P),\qquad
 \rho^2=q_K-k_0^2>0.
\]

The orthogonal complement of A/P is a positive two-plane. For an
orthonormal basis \(e_1,e_2\), **every** norm solution with those two
contacts is

\[
 K(t)=k_0+\rho(\cos t\,e_1+\sin t\,e_2),\quad0\le t\le2\pi.
\]

Thus both physical sides are included without choosing a single norm
branch. The source computes and checks the positive complement and
orthogonality before using this parameterization.

Every old/native acuteness condition, the common-center condition and
each intrinsic nonincidence bound is linear in \(\cos t,\sin t\).
The stipulated K face has only A,L,P,R as sides. Accordingly, for every
retained intrinsic old side and C it imposes the necessary bound
\(K\cdot n_i\le-\sqrt{q_Kq_i}\). A positive-definite pair of actual
acute side normals would give an actual edge by strict-center projection,
contradicting this incidence. Causal old projections are included in the
native checks but are not incorrectly treated as intrinsic side facets.

For either required finite corner A/L/K or P/L/K, the Schur complement
gives a closed necessary interval. With the first root of square \(q_0\),
L of square \(q_L\), and fixed products \(g_{0L},g_{0K}\), it is

\[
 \frac{g_{0L}g_{0K}-
 \sqrt{(q_0q_L-g_{0L}^2)(q_0q_K-g_{0K}^2)}}{q_0}
 \le g_{LK}\le
 \frac{g_{0L}g_{0K}+
 \sqrt{(q_0q_L-g_{0L}^2)(q_0q_K-g_{0K}^2)}}{q_0}.
\]

The source intersects both corner intervals with the native finite-contact
range \([-2,0]\), retaining its parallel endpoint as a relaxation.
All these inequalities cut closed arcs from the full circle. Numerical
half-planes are slightly enlarged; interval endpoints and isolated
boundary points are retained. This enlargement and the numerical endpoint
calculation do not replace rigorous interval enclosure.

For each R, the required finite A/R/K and P/R/K classifications supply the
permitted K–R labels2 through5 and any symbolic label-at-least6 family.
The same Schur formula bounds the projected K–R product. Its range on each
arc is obtained from endpoints and stationary points. Every allowed low
label is then solved on the circle. A high K–R family is relaxed to its
whole product interval \([-2,-\sqrt3]\), rather than sampled at a cutoff.
No such high interval is reached in this run.

## Symbolic high A/P contacts

The inherited finite-corner classification has33 high K–A/K–P families
at these frames. Every one has zero K A2-profile and exactly one high bond;
the other A/P contact has label at most5.

If K–A has label at least6, finiteness of A/L/K forces A–L and K–L right,
and finiteness of A/R/K forces A–R and K–R right. The analogous conclusion
holds for a high K–P bond. This follows from the norm2 three-root
determinant: a bond of magnitude at least \(\sqrt3\) sharing a node with
any second nonright bond of magnitude at least1 makes that determinant
nonpositive. The equality case is not finite.

Consequently every surviving high family has K right to L and R. Its
remaining low A/P contact defines the three-frame problem
\((L,R,\text{other A/P root})\). Both norm branches are solved, and the
resulting high product is tested against the entire closed interval
\([-2,-\sqrt3]\) and then the native Coxeter angles. This removes the
high-label cutoff for this necessary test. Singular unresolved frames
would be recorded; none occurs here.

## Saved diagnosis

There are837 distinct low K–A/K–P recipes. Of their relaxed circles,625
are nonempty and212 empty. The first empty constraint is old7 native
acuteness in30 cases, old7 intrinsic nonincidence in139, old53 native
acuteness in10, and L native acuteness in33. These first-blocker counts
depend on the checking order and are not unique obstruction assignments.
The circle reconstruction also recovers all454 predecessor K caps whose
two A/P labels are at most5.

The low domain gives4,895 R/circle tests and1,016 finite K–R events. Every
K–L product is in a native angle gap; the range is approximately
−1.49725924 to−0.00007904819. The full native blocker sets are:

| Blockers | Events |
|---|---:|
| K–L alone |263|
| old7 and K–L |620|
| old7, old30 and K–L |124|
| old7, old10, old30, old52 and K–L |9|

The closest K–L value to an allowed angle occurs at selected event1352,
R4430: K is right to A,P,R but has K–L approximately−0.00007904819,
close to, but distinct from, a right angle. This is an oldP3 control event,
not a validated deformation or a new candidate.

For the high A/P families,247 R/family combinations fail a required
finite-corner classification. The121 surviving frame equations give242
norm branches. Only three branches lie in the high product interval.
All three have product approximately−1.8327004801, strictly numerically
between \(-2\cos(\pi/8)\) and \(-2\cos(\pi/7)\); the high bond itself
is their sole native blocker. None yields a candidate. There are no saved
unresolved records.

The diagnosis therefore concerns simultaneous native-angle constraints.
The continuous A/P circle can reach a finite K–R angle, while the required
K–L angle misses the Coxeter set. Increasing the sampled K–L label cutoff
does not address the low-A/P events found here: all their K–L magnitudes
are below1.498. A genuinely new family would have to change the coupled
A/P/L/R geometry or the required adjacency, while preserving the old
native constraints; this report does not assert that such a change exists.

## Scope and verification

The fixed input frames are still numerical moving-wall events. The circle
and high-bond reductions are mathematical necessary conditions, but arc
intersections, eigenvectors, Schur intervals, level roots and native-angle
separation use floating point. Near-degenerate numerical behavior is not
independently interval-certified. No claim excludes all labels by exact
arithmetic or all moving parameters. The result does not cover the other
unselected events, different P/A families or different disk incidences.

All62 retained old walls are checked in the necessary domain; the inherited
exact deletion debit remains two old finite supports and the old12
quadrilateral. Every returned finite K–R point rechecks native products,
the strict center, norm/contact residuals and intrinsic nonincidences.
There is no full native candidate, so the later complete-cone/PD scoring
branch is not reached. The previous frozen package remains the source
of the4,128 fully scored K/R cones, all with ten open ends. This diagnosis
does not claim a new ambient finite-vertex census or finite-volume chamber.

**FINAL/FROZEN.** Generation and a separate fresh saved replay both
terminated exit0/PASS. Eleven direct source/data bindings, Python syntax
and three report links pass. A standard-library saved-data check importing
no producer verifies all83 domains,1,016 low angle events,121 high frame
records and242 branches, the native gaps, counts and residual/center
guards. It does not reconstruct physical projections, construct circle
arcs or solve the level equations.

Root's independent full source/report mathematical read passed within the
numerical scope: the Lorentz A/P projection and full positive circle,
enlarged half-plane arcs, Schur intervals, low K–R roots and high necessary
interval, high A/P isolation and finite-frame reduction, and the explicit
residual/native/center/nonincidence rechecks. It did not independently
recompute the numerical arcs or physical coordinates. No frozen
predecessor is modified.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_ap_circle_20261017.py
OPENBLAS_NUM_THREADS=1 python scripts/fresh_cu_release47_ap_circle_20261017.py --verify
```
