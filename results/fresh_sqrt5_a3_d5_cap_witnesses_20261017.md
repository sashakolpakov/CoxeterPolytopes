# Two remaining finite-component cap obstructions

Two exact additional capacity constraints violate the [final12.12966051709919 witness](fresh_sqrt5_triple_groundfield_cap_joint_ports_20261017.md). The selected A1+A3 profile has capacity2 while its quadrilateral mass is three times the containing-simplex mass; the A1+D5 profile has capacity1 while its quadrilateral mass is twice the simplex mass. No repair or optimizer is performed here.

The [source](../scripts/fresh_sqrt5_a3_d5_cap_witnesses_20261017.py) and [certificate](fresh_sqrt5_a3_d5_cap_witnesses_20261017.json) reconstruct literal base4465/type110693 and base5885/type8351. Both are maximal-edge bases. Their full containing-simplex stars are nondegenerate of rank14. Hence a new norm2 root is uniquely determined by its complete contacts with that star. The source enumerates every image of the selected core contact, every allowed finite-component contact position, and every identical A1 spectator copy. Distinct placements represent distinct roots; an identical root is one incidence. A fixed extra root determines one full canonical quadrilateral star, so it is not counted as several independent caps at the same containing simplex.

These are constraints on two specified actual provenance rows. Each such quadrilateral occurrence supplies the displayed base; ignoring its other bases weakens the inequality. No merged provenance is substituted and no completeness assertion about all79 untouched A1+A3/D4/D5 profiles is made.

## A1+A3: capacity two

For base110693, the Lorentz core has four norm2 roots. Roots0,1 have product−(1+sqrt5); each has product−1 with roots2,3, which are mutually orthogonal. Its two spherical deletions are0,1. The cap has product−1 with exactly one of these roots, giving two core orientations. The normalized inverse entries are

\[
J_{00}=J_{11}=(2-\sqrt5)/4,\qquad J_{01}=-1/4.
\]

The cap also has product−sqrt(tau) with one of the three isolated A1 copies, where tau=(3+sqrt5)/2, and product−1 with an endpoint of the unique A3 component. The endpoint projection has norm3/4, and the two endpoint projections have product1/4.

Two distinct placements in the same core orientation have positive mutual product, whether their A1 copy agrees or differs. Across the two orientations, compatibility is possible only when the A1 copies differ and the A3 endpoints are opposite; the mutual product is then zero. Thus at most two placements coexist. The exact twelve-position compatibility graph and complete maximum-clique recurrence attain2. The resulting necessary inequality is

\[
\alpha_{Q61927}\le2\alpha_{S81801}.
\]

The saved masses give slack−1057063/30000000000.

## A1+D5: capacity one

For base8351, the Lorentz core has two norm2 roots with mutual product−sqrt(t), t=3+sqrt5. The cap has product−1 with one endpoint, so again there are two core orientations. Its D5 contact is the unique long-arm endpoint, whose projection has norm1 and is unchanged by exchanging the core orientations. There are two A1 copies, each with cap contact−sqrt(t).

In the same orientation, two distinct A1 placements have product

\[
1+\frac2{4-t}=(1-\sqrt5)/2\in(-1,0).
\]

This is forbidden for distinct Coxeter normals of norm2: the possible finite products begin at0 and−1, while parallel or divergent walls have product at most−2.

In opposite core orientations the product is

\[
A+\frac{\sqrt t}{4-t},\qquad
A=1\ \text{or}\ 1+t/2,
\]

according as the A1 copies differ or agree. Here sqrt(t) is outside K=Q(sqrt5), while both A and1/(4−t) are nonzero in K. Its square therefore lies outside K, violating the arithmetic two-cycle ground-field condition. All pairs of the four placements are forbidden, so capacity1 is exact. The resulting inequality

\[
\alpha_{Q75841}\le\alpha_{S78789}
\]

has slack−114553/5000000000.

Both literal finite Gram inverses, core contact maps, root automorphisms, complete position graphs, pair signs/field tests and clique recurrences are saved. Generation and a fresh reconstruction replay pass. A separate narrow conceptual check confirms the D5 product and field mechanism, conditional on the literal identifications which this source reconstructs. This is a local compatibility result in the fixed arithmetic field, not a global chamber realization or dimension theorem.
