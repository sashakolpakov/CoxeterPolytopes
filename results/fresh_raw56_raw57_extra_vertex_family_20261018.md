# Complete raw56/raw57 extra-vertex families

This packet proves new coupled inequalities for the current **compact
arithmetic $H^{13}/\mathbb Q(\sqrt5)$** model. It covers all 704 spectator
classes of raw56 and all 704 of raw57, including zero-mass classes. It does
not change a historical ordinary-corner graph or run an optimizer.

The resulting 1,212 one-vertex columns have 60 strict violations at the
unchanged current face masses. Their least monotone exported-column repair
has the exact ratio

\[
\frac{14671331978990040}{1220681360116301}
=12.01896945292277\ldots>12.
\]

The predecessor ratio is $12.273467266364642\ldots$. These are exact
fixed-mass certificates. The new number is neither a complete enlarged-cone
optimum nor a global optimum, and it gives no new dimension bound. The
completed unrestricted bound remains 59. The predecessor's sole global
optimization attempt timed out; this packet does not retry it or infer new
face masses.

The source, saved data and independent saved-data checker are
[producer](../scripts/fresh_raw56_raw57_extra_vertex_family_20261018.py),
[certificate](fresh_raw56_raw57_extra_vertex_family_20261018.json) and
[checker](../scripts/check_raw56_raw57_extra_vertex_saved_20261018.py).
The bound predecessor is the
[raw76 full-model interval](fresh_raw76_complete_family_full_model_20261018.md).

## Two actual geometric incidences

Raw56 and raw57 are canonical six-wall cores 45 and 46. Each has exactly
one quadrilateral partition, with opposite side pairs $\{1,3\}$ and
$\{0,2\}$, and fixed core walls 4 and 5. Their passive finite spectator
forest has rank nine. Both scaled six-wall Grams have inertia $(4,1,1)$.
The complete 15 possible four-root supports have precisely five positive
definite members, with omitted pairs

\[
\{0,1\},\{0,3\},\{1,2\},\{2,3\},\{2,4\}.
\]

The first four are the source quadrilateral's four finite vertices. The
last is an additional support $\{0,1,3,5\}$. Its whole active finite
forest is $2A_1+I_2(10)$ for raw56 and $A_2+I_2(10)$ for raw57. With
the spectator it has rank 13.

These are **actual vertices**, not merely positive principal submatrices.
For completeness, let $h$ be a strict interior future timelike vector of
a realizing Coxeter polytope. For selected actual facet normals $n_i$,
write $G=((n_i,n_j))>0$, $b_i=(h,n_i)<0$, and $c=G^{-1}b$. The
Gram is acute, so its inverse is entrywise nonnegative and $c_i\le0$.
Then

\[
x=h-\sum_i c_i n_i,
\qquad x^2=h^2-b^TG^{-1}b<0.
\]

It is orthogonal to the selected normals. For any other actual facet
normal $m$, acuteness gives

\[
(x,m)=(h,m)-\sum_i c_i(n_i,m)\le(h,m)<0.
\]

Also $(x,h)=x^2<0$, so the ray is future timelike. Thirteen independent
normals leave exactly this timelike line. Every additional facet is
strictly inactive there. This proves a genuine finite vertex even in the
presence of arbitrary extra facets. At the fifth vertex, wall 4 is
strictly inactive, so this vertex lies outside the source quadrilateral,
whose fixed support includes wall 4. This is the same strict-center
criterion proved in the frozen
[1325 incidence packet](fresh_forest1325_cross_vertex_incidence_20261017.md).

Both cores also force an actual compact triangle. Fix walls 0 and 5,
together with the spectator. The three remaining sides are 1, 3 and 4.
The saved Schur Gram has inertia $(2,1,0)$, strictly negative inverse,
and all three pairwise side intersections are among the five actual
finite vertices just proved. Its three sides therefore bound a compact
triangle. The polytope's face lies inside this triangle, while convexity
and the three actual vertices put its entire convex hull in the
polytope. Thus the face equals that triangle; further facets cannot
truncate it. At the fifth vertex its two ports are 1 and 3. An ordinary
quadrilateral occupying those same ports is impossible, since the two
ports of a simple finite vertex determine a unique actual two-face.
This is a compatibility exclusion only. No additional triangle demand
or triangle allowance is introduced.

## Recovery and multiplicity

A marker records the complete active-role map into the actual finite
13-wall frame, the full native contact covectors of omitted walls 2 and
4 (zero on every spectator coordinate), their norm-two normalization,
and their strictly positive timelike heights. The full finite Gram is
positive definite. A contact covector fixes a missing wall's spacelike
projection uniquely; norm two fixes the magnitude of its remaining
timelike coordinate. The common future-side sign fixes that coordinate
uniquely. Thus the marker recovers both actual missing walls, not merely
an abstract diagram congruent to them.

The recovered six core walls and the actual spectator walls determine
the source quadrilateral's fixed support and its four sides. Every
source quad has exactly one fifth support, and the saved marker keys are
duplicate-free unordered pairs of missing full-frame covectors. Complete
whole-component embeddings enumerate all possible active role maps. A
permutation of unused spectator components changes no covector or port
and is not another occurrence.

Each connected six-wall native diagram has only the identity
automorphism, and the two native diagrams are not isomorphic. Each has
only the one quadrilateral partition stated above. A compact actual quad
has a unique full supporting-wall set: at its simple finite vertices,
the fixed and side walls are exactly the incident facets. This is the
frozen [full-star aggregation interface](fresh_quad_product_aggregation_20261017.md).
Consequently two different descriptions of that same quad cannot create
a second marker, nor can the same quad count in both raw56 and raw57.
The connected indefinite six-wall core also distinguishes its walls from
the passive finite components. The incidence coefficient is therefore
**one marker per actual source quad**, independent of orientation and
component stabilizers.

The complete rank-nine spectator catalogue gives 704 classes in each
family. There are 1,212 distinct target forests because 196 targets admit
both families. The 21,584 marker orientations comprise 16,472 raw56 and
5,112 raw57 orientations. These orientation counts are not face-demand
multipliers. Only 45 raw56 and 30 raw57 classes have positive mass at the
current witness; all 1,408 classes are included in the inequalities.

## Exact $I_2(10)$ projections

The previous native projection helper handles active 3/5 components and
skips non-$K$ components when they are passive. Here $I_2(10)$ is
active, so this packet adds its exact inverse explicitly. With
$\tau=(5+\sqrt5)/2$, its norm-two Gram and inverse are

\[
\begin{pmatrix}2&-\sqrt\tau\\-\sqrt\tau&2\end{pmatrix},
\qquad
\frac1{4-\tau}
\begin{pmatrix}2&\sqrt\tau\\\sqrt\tau&2\end{pmatrix}.
\]

Every marker's self-projections, positive heights and mutual omitted-wall
product are reconstructed in a fresh immutable squareclass field. The
ordinary nodes in all affected current domains are explicitly checked
to have zero contacts on the $I_2(10)$ blocks. No cached field
coefficients migrate between wall pairs. The existing numeric native
decoder represents the extra magnitude $\sqrt\tau$ by its squared
$K$ value; this does not reclassify a Coxeter-10 edge as a dotted edge.

The saved checker takes a different route: it scales one root of each
I2(10) block by the positive square root of tau and inverts the resulting
K-valued contacted-block Gram by Gauss–Jordan elimination. Untouched
passive blocks with other labels are replaced by harmless diagonal
blocks; explicit zero-contact and whole-component guards ensure that
they contribute to no checked projection. The contact covectors are
scaled consistently, so the inverse is exact on all contacted coordinates.
The checker independently implements ordered
$K$ and multiquadratic arithmetic, signs, squareclass membership and
quadratic integrality.

## Complete marker compatibility from 31 cases

Write $a$ for the number of whole $A_1$ components and $h$ for the
number of whole $I_2(10)$ components at a target vertex.

A raw56 marker chooses an ordered pair of distinct $A_1$ roles and an
oriented $I_2(10)$ component. Two ordered $A_1$ pairs have seven
relative equality patterns. Their $I_2(10)$ roles are equal with the
same orientation, equal with opposite orientation, or different. After
removing the identical marker, there are 20 cases. The exact pair gates
retain precisely:

- the same ordered $A_1$ pair and different $I_2(10)$ components;
- disjoint $A_1$ pairs and the same oriented $I_2(10)$ component;
- disjoint $A_1$ pairs and different $I_2(10)$ components.

Thus distinct ordered $A_1$ pairs in a clique form a matching, and each
such pair occurs at most once for each $I_2(10)$ component. The marker
capacity is

\[
c_{56}=\lfloor a/2\rfloor h.
\]

A raw57 marker chooses an oriented whole $A_2$ and an oriented whole
$I_2(10)$. The eight nonidentical relative cases retain only the same
oriented $A_2$ with different $I_2(10)$ components. Hence

\[
c_{57}=h.
\]

For a mixed raw56/raw57 pair, the $A_1$ and $A_2$ roles are disjoint
whole components. The three relative $I_2(10)$ cases all exclude the
pair. The two marker families cannot coexist at one target vertex.
Three rank-13 host forests realize all 20+8+3 relative patterns; passive
orthogonal components cannot affect a product. Each family's capacity
also has an explicit attaining clique in the conservative marker graph.
This is graph attainment, not simultaneous geometric realizability.

The exact necessary gates are positive product, the forbidden open
norm-two product interval $(-1,0)$, and failure of the squared product
to be an algebraic integer of $K$. Every surviving gate outcome is
retained as an unknown.

## Ordinary-domain transport and one-budget columns

The current chronological 327-class ordinary domain has 600 target
forests. Ninety-one of the new targets contain ordinary nodes, with
2,883 such nodes in total. A full direct comparison would have 108,604
marker–ordinary positions and 577,768 marker–marker positions.

For each represented target and each available new family, the producer
tests one anchor marker against every ordinary node. Every other marker
is reached by an explicitly saved automorphism of the complete finite
frame. It permutes all native ordinary fields and preserves each
ordinary orbit's owner, hence every actual-origin coefficient. It also
permutes the complete marker set. This transports the anchor decisions
without requiring the historical ordinary–ordinary graph itself to be
automorphism invariant. All old ordinary edges, exclusions and unknowns
remain unchanged.

There are 4,115 anchor tests. Their 414 unknown outcomes transport to
13,104 allowed marker–ordinary positions. None is silently excluded.
All 31 marker-relative cases are tested separately from these anchors.
This is complete structural coverage, not a positive-mass sample or a
new ordinary-pair census.

At a target $v$, take any bound nonnegative old ordinary charge $w$
whose load is at most one on every saved complete ordinary count vector.
For the new charge only, set to zero every coordinate of an ordinary
orbit compatible with any marker. Call the result $w'$. All retained
positive ordinary coordinates are then incompatible with every new
marker. If $k_{56},k_{57}$ are actual marker counts, the preceding
capacity and mutual-exclusion theorem gives

\[
w'\cdot q+\frac{k_{56}}{c_{56}}+\frac{k_{57}}{c_{57}}\le1,
\]

omitting a term when that family is absent at $v$. When a marker is
present the ordinary term is zero and only one marker family can occur;
when none is present the old ordinary bound applies. Nonnegative weights
extend the inequality from maximal ordinary vectors to every subclique.

Ordinary coefficients count each actual corner origin once, and marker
coefficients count each actual source quad once at its unique fifth
vertex. No base-triangle count is added. The entire inequality has one
vertex cost, not separate allowances for the ordinary and marker parts.
Existing ordinary or raw76 columns remain independent valid inequalities;
their budgets are not removed or reduced. The export contains 1,212
distinct one-budget columns, one per target forest.

## Exact repair and limits

Keep every face mass fixed and raise each vertex multiplier to the
largest of its old value and all new column loads there. Exactly 60
targets rise. This is the least monotone repair for these exported
columns. The source explicitly checks all 242 inherited raw76 columns
and all 14 multi-owner columns with their full costs.

The remaining inherited matrix constraints are also preserved. Every
nonzero vertex-budget coefficient in the frozen model encoding is
positive. Holding face multipliers fixed preserves every zero-budget
left side exactly and increases every nonzero-budget right side, including
the constraints attached to signed variables. This is the inherited
`cost_rows`/`verify_dual` interface used by the bound full-model replay;
no sign of a face coefficient is being assumed. Thus monotone repair
preserves all 440,098 predecessor columns, including all historical
unknowns and the 242 raw76 marker columns.

The reported ratio remains above 12. It is a valid witness for the
new exported model, not proof that no stronger mixed inequality can
help. The columns are valid for arbitrary face masses, but no new global
optimization, complete enlarged local cone, or new physical chamber is
claimed. No ordinary class is declared fully represented merely because
its extra-vertex marker has been covered.

## Verification status

Final source-bound generation passed in session54545, terminal chunkd4768e.
Separate fresh reconstruction passed in session62410, terminal chunk3af139.
The independent saved checker passed in session32167, terminal chunk85b594:
31 relative patterns, 21,584 marker maps, 4,115 anchor tests with 6,515
wall expressions, 1,285 ordinary charge inequalities, 1,212 columns,
242 preserved raw76 loads and all 14 multi-owner loads. A docstring-only
scope correction also replayed successfully in session98914, terminal
chunk4d4179. The final helper with duplicate-target completeness guards
passed in session7408, terminal chunkd1e77d, with the same counts. All
processes exited zero; no job remains live. The final metadata check
passes all 55 direct bindings, both Python syntax trees and all six local
report links. Source, data, report and helper are final/frozen.

Root completed a full producer/report/helper mathematical read of the
actual fifth-vertex and triangle arguments, full normal recovery and
unique quad count, 31 patterns and matching capacities, owner-preserving
transport, unknown retention, one-budget charges, prior raw76/multi-owner
loads and monotone full-matrix preservation. The construction peer
completed a separate full producer/report/helper read and reread inherited
active embeddings, projection and pair-field gates, chronological current
domains/native transport, and full-model `cost_rows`/`verify_dual`.
Both reviews passed after the contacted-block scope correction and
duplicate-target helper guards. Neither reviewer executed or repeated
the physical, catalogue, ordinary clique, optimizer, matrix or saved
arithmetic computations. The root review does not claim an additional
inherited-helper reread. These are internal mathematical reviews, not
external refereeing or proof-assistant formalization.

The independent checker replays the two complete 15-support PD decisions,
triangle residual inverses, native nonisomorphism/automorphisms, all
complete active-role maps, future heights and omitted products, every
new anchor gate and all 31 relative cases, saved full-frame/owner
permutations, old count-vector charge inequalities, coefficients and
exact repair arithmetic. It inherits historical owner chronology,
ordinary clique completeness, complete spectator-catalogue provenance
and the old full-matrix theorem from the directly bound packets. It does
not claim another physical historical census, local optimization or
global matrix replay.
