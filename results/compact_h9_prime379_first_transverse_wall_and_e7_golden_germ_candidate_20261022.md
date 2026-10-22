# Norm379: the first transverse wall and its full E7 + I2(5) germ

22 October 2026. Root-authored mathematical/runtime composition, **pending
complete outside-author review**. The five profile units and separate combine
have actual bounded PASS records. This note identifies one new actual edge
endpoint of the literal lattice. It does not prove finite volume, reflectivity,
a chamber period, a new polytope example or a dimension-bound improvement.

## Statement and exact mathematical boundary

Let phi=(1+sqrt(5))/2, O=Z[phi], D=1+20phi and b=20+21phi. Use the ORIGINAL
module L=E8(O) perpendicular [[-2phi,1],[1,10]], in coordinates
(X,Y,e0,...,e7). Write

    q=(18+25phi)X-(12+17phi)Y,
    m=(2+5phi)X-(2+3phi)Y,
    H=(2,3,4,6,5,4,3,2) in the E8 summand.

The accepted first binary endpoint has full outward simple germ
(-e0,...,-e7,m). Release -e7. The first new simple wall on that actual edge is

    r=a q-(phi/2)m+z H,
    a=2572+(8323/2)phi, z=14856+24036phi.                 (1)

It has square2. Its collision parameter on q+tau H is

    tau=(9575+23982phi)/8376.                            (2)

The full new finite simple germ is

    (-e0,...,-e6,m,r), of type E7 + I2(5).              (3)

This statement uses the complete first-binary/folding theorem and its
explicit original-module, norm-alphabet, finite-local-chamber and arithmetic
premises. It neither assumes a folded submodule is the original lattice nor
transfers an outcome from the different norm239 lattice.

The exact source-bound mathematical domain is
`results/compact_h9_prime379_first_binary_edge_and_quarter_folded_transverse_domain_candidate_20261022.md`,
SHA-256 `7e7ffd616d0e025aeaff84abae1c9fcfb7d930c8a3879a7bb43734f308cacd3a`.
Its root acceptance is
`results/review_prime379_first_binary_and_quarter_folded_domain_root_acceptance_20261022.md`,
SHA-256 `0f79f68d291136cd7d1963c70dd3d92d04246fbc55e9ffbd4c80a48a40497a84`.
Its complete outside mathematical review is
`results/review_prime379_first_binary_and_quarter_folded_domain_full_scope_foundations_20261022.md`,
SHA-256 `28fc56c79de398092083bc42c5ab721fa218451762715632bf6d319e334dd7d5`.

## Algorithm, termination and completeness

The algorithm is specified in full in
`results/compact_h9_prime379_quarter_folded_transverse_checker_spec_20261022.md`,
SHA-256 `deaa343294d09e35ab89cc48d84eba1ddfb9ad5f337940c5c01fa49422d98c10`.
The implementation is
`scripts/check_compact_h9_prime379_quarter_folded_transverse_batches_20261022.py`,
SHA-256 `e6310c0a345d682b7769a4d6d7feb9dedd52b3229582c67858300900a426b7cd`.
Both were completely read and mathematically reviewed before execution;
the root source acceptance and complete outside code review are bound by the
runtime manifest below. The preparation-era source headers are preserved.

Here is the finite procedure whose output is used. The five possible first
simple-wall profiles (d,kappa,epsilon) are

    (1,0,0), (1,1,0), (1,phi,0), (1,0,1), (b,0,0).

The accepted domain proves every first simple wall lies before the actual
original integral mirror V^10 e7 and is covered by a nearest phase j in
{0,...,10} of the explicitly proved quarter-integral boost. The WHOLE folded
seed is quarter-integral; its time coefficient is not separately restricted.

For each short profile and phase, enumerate -30<=e<=30 and all integers c
in the closed conjugate strip |c+e phi'|<=4, with 4z=c+e phi. Retain the
proved physical/conjugate and folded-window bounds. For the long profile use
z=b g, 4g=c+e phi, -4<=e<=4 and |c+e phi'|<=8/5, then retain the corresponding
proved bounds. There are at most 549 short or36 long pairs per phase, hence at
most24,552 across all55 selectors. Solve each norm equation exactly in Q(phi),
including both time signs and all rational/pure-sqrt(5) branches. Each pair has
at most two roots. These finite bounds prove termination independently of
finding any wall.

For each root, reconstruct the entire original ten-coordinate vector, unfold
it by its phase and check original integrality, exact norm, the FULL integral
reflection matrix, outgoing orientation, positive collision not beyond the
original mirror upper, and acuteness against all thirteen known original
walls. The norm/alphabet premise supplies original primitivity. A profile unit
stores its minimum and EVERY distinct original tie; a tie or storage overflow
is a failure, never a truncated successful census.

A separate combination authenticates all five part hashes and all55 distinct
selectors, revalidates every saved tie in the original coordinates, and takes
the global exact minimum. It tests all candidate simple cones using the eight
retained walls plus a new wall: full Gram positive definiteness, endpoint
orthogonality, and nonnegative expansion of EVERY minimum tie checked in ALL
ten coordinates. A unique cone is necessary for the reported full-germ PASS.
The code review proves these operations implement the stated finite algorithm;
the actual execution records supply its output as a separate premise.

Completeness at the mathematical boundary has both directions. Every first
SIMPLE wall of the full original chamber is in the accepted five-profile
domain. Conversely an original arrangement mirror cannot cut the interior of
that chamber before its first boundary. Therefore the complete admitted
minimum identifies the first collision. The full local cone then identifies
the simple germ. Checking only thirteen walls, or merely finding an integral
vector of the right norm, would not establish that conclusion.

## Actual execution and exact output

The immutable runtime manifest is
`results/prime379_transverse_runtime_manifest_20261022.json`, SHA-256
`ae53f5fc76d9d62c28c6bc4a905ce094acb4c6165c483cea9a3c6f11357e8c20`.
It binds all five result files, all five actual wrapper ledgers/stdout files,
the separate combine result/ledger/stdout, current worker/runner/contract and
both full source reviews. Every recorded binding matches current bytes.

The five units completed all55 selectors and actually visited21,835 scalar
rows;20,559 were inside the strips and14,795 inside the folded window. There
were132 field time roots/quarter seeds,28 original-integral occurrences and12
admitted occurrences. These are occurrence counts across phases, not distinct
root counts. The combine has exactly ONE distinct global minimizer, in profile
`norm2_golden`, phase2. Its output is
`results/prime379_transverse_combine_20261022.json`, SHA-256
`861e3b5385339946ab1691e70ed3d9bfbf37932a1b814b6deb9ddbb6d958fce3`.
It reports `PASS_COMPLETE_DOMAIN`,55 selectors and one unique certified cone.

For direct coordinate comparison, the output wall has X,Y coefficients

    (150331+243241phi, -101608-164405phi),

and E8 part z H. These are exactly (1): multiplying a by q and subtracting
phi*m/2 gives those two binary coefficients, while zH has first coordinate
29712+48072phi, with all remaining coordinates in the displayed H ratios.
The combiner checks its original square, reflection matrix and collision
exactly; substituting (1) gives q.r=-2Da and H.r=2z, so the parameter is Da/z,
which simplifies to (2).

All six jobs ran serially through the unchanged shared-lock resource runner,
with59-second/1-GiB sampled process-tree limits and one thread. Each passed
40-GiB disk/6-GiB available-memory admission;8-GiB/3-GiB abort floors remained
active. The largest elapsed time was0.628s, and maximum sampled process-tree
RSS was59,801,600 bytes. No scientific worker was imported outside the runner,
and no successful unit was rerun. No build products were generated; parent
and child used Python -B. All result and stdout evidence is retained.

## Full-germ interpretation and limits

The eight retained normals are -e0,...,-e6,m. The first seven have the usual
E7 Gram and are orthogonal to H,q,m, hence to r. The binary pair has

    m²=r²=2, m.r=-phi.

Its angle is pi/5, and its Gram determinant4-phi² is positive. Consequently
the output's complete nine-normal Gram is precisely E7 plus I2(5), as in (3).
The actual negative future endpoint q+tau H is orthogonal to all nine normals.
The checked unique complete cone and the accepted firstness/local-recognition
premises make this the FULL finite germ, not just a finite subdiagram found
among known walls.

There is no conclusion about a finite-volume chamber: other outgoing edges
have not been completed. Nor has any future integral isometry or period of
the chamber been proved here. This gives one actual step of the literal
norm379 edge walk. It supplies no field-wide/form-class/lattice-family
exclusion, new example or improved dimension ceiling. Lean was not run.
