# Norm379: the released-m edge and its full E7 + 2A1 germ

22 October 2026. Root-authored mathematical/runtime composition, **pending
full outside-author review and root acceptance**. Five newly executed
profile units and a separate authenticated combination have actual bounded
PASS records. This identifies one more finite edge endpoint in one literal
reflection chamber. It proves no finite volume, reflectivity, chamber
period, new polytope example or dimension-bound improvement.

## 1. Statement and exact mathematical inputs

Let phi=(1+sqrt(5))/2, O=Z[phi], D=1+20phi and b=20+21phi. Work in the
ORIGINAL module

    L = E8(O) perpendicular [[-2phi,1],[1,10]],

with coordinates (X,Y,e0,...,e7). Put

    q=(18+25phi)X-(12+17phi)Y,
    m=(2+5phi)X-(2+3phi)Y,
    H=(2,3,4,6,5,4,3,2),
    r=a0 q-(phi/2)m+z0 H,
    a0=2572+(8323/2)phi, z0=14856+24036phi,
    tau0=(9575+23982phi)/8376.

The accepted starting future ray q0=q+tau0 H has FULL outward simple
germ (-e0,...,-e6,m,r), of type E7+I2(5). Define

    u=m+r, d=2-phi, F=s_r s_u,
    s_r(v)=v-(v.r)r, s_u(v)=v-[(v.u)/d]u.

The accepted domain theorem proves that F and its inverse preserve the
entire original module, Fm=r, Fr=-m-phi r, F fixes q0 and E7, and F has
order5. No chamber invariance of F is assumed. Release m at q0 and
parametrize the actual outgoing edge by x(tau)=F(q+tau H), tau>tau0.

**Theorem at the explicit inherited geometric and exact-output boundaries
below.** The first new simple wall of that edge is the original vector

    w=(28490145+46098027phi)X
      -(19256316+31157373phi)Y
      +(2815372+4555369phi)H.                         (1)

Its square is2b and its collision parameter is

    tau1=(1154937+1948905phi)/745919.                 (2)

The full new finite simple germ is

    (-e0,...,-e6,r,w), of type E7 + A1 + A1,          (3)

with r²=2, w²=2b and r.w=0. The two A1 factors have different original
root norms; the diagram type alone does not identify their lattice data.

The exact complete domain is
`results/compact_h9_prime379_third_edge_released_m_complete_folded_domain_candidate_20261022.md`,
SHA-256 `0ff2e60377f26528378c9803024765fc443c7503ae071cc797936d49db3f56ac`.
Its full outside mathematical review is
`results/review_prime379_third_edge_released_m_complete_domain_full_scope_compact_reviewer_20261022.md`,
SHA-256 `5e578817b4b08f3b99720f8587c78528d3b3ebcfefd5ff68ffc5c6d422a4e2c3`;
its root mathematical acceptance is
`results/review_prime379_third_edge_complete_domain_root_acceptance_20261022.md`,
SHA-256 `a1d2763bf66faa04d3f91a65633ee05a5708a0ca60936053f9196a3096ed67cb`.
That domain pins the accepted starting full germ, original glue, primitive
norm alphabet2/2b, finite/native diagram classification, local finiteness
of the full original integral-reflection arrangement, and finite local
chamber-recognition premises. These remain explicit mathematical inputs.

## 2. Algorithm and why it terminates and is complete

The NEW checker contract is
`results/compact_h9_prime379_third_edge_released_m_checker_spec_20261022.md`,
SHA-256 `6b43db1fcd1d38e440288c371a48877c4e41c38c8e018df32e5df666353f42f0`.
Its standalone implementation is
`scripts/check_compact_h9_prime379_third_edge_released_m_batches_20261022.py`,
SHA-256 `56d8294bc78fa037de7413ac076799bf572da09659bac1ae5ca86301abe2bd7d`.
Their entire source/contract received full other-author and root reviews
before execution. Those two review files are bound by the runtime
manifest in Section3. Preparation-era headers are preserved history.

Here is the mathematical algorithm whose actual output is used. The
domain supplies the original integral mirror upper F V^10 e7 at a
strictly later finite parameter tauU<sqrt(D), where V is the explicitly
proved quarter-integral boost. Every first simple wall w has inverse
image v=F^-1 w of one of the five profiles

    (d_w,kappa,epsilon) =
    (1,0,0), (1,1,0), (1,phi,0), (1,0,1), (b,0,0).

For each profile and j=0,...,10, enumerate its nearest-phase folded
seeds V^-j v. Let Z0=Omega6-3H/2, with
Omega6=(4,6,8,12,10,8,6,3). A seed has the unique expression

    a q-(kappa/2)m+z H+epsilon Z0,
    a²=[2z²+kappa²/2+(3/2)epsilon²-2d_w]/(2D).

For short profiles enumerate4z=c+e phi, -30<=e<=30, and every integer c
in the closed conjugate strip |c+e phi'|<=4. Retain the domain's physical,
conjugate and folded-window inequalities. For the long profile put z=b g,
enumerate4g=c+e phi, -4<=e<=4, |c+e phi'|<=8/5, and retain its bounds.
There are at most549 short or36 long scalar pairs per selector, hence
at most24,552 across all55 phase/profile selectors. Each scalar gives at
most two exact field roots for a. Both signs, zero and rational/pure-sqrt5
branches are retained. These bounds prove termination independently of
finding a successful candidate.

For every seed reconstruct ALL ten original coordinates, check whole
quarter membership, unfold and apply F. Check original integrality,
exact norm2 or2b, primitivity by the even original module and prime norm
alphabet, and the FULL integral reflection matrix. Require the outgoing
derivative positive and tau0<tau<=tauU. Test acuteness against all FOURTEEN
known original actual walls. A known wall is exempted as self only after
full ten-coordinate positive proportionality is proved. Neither the old
untransformed filter nor a zero Gram determinant is used for self recognition.

Each independent profile stores its exact minimum and EVERY distinct
minimum tie. A tie/storage overflow is a failure. The separate combination
authenticates all five part hashes, exactly55 distinct selectors and the
same current source/premise bindings; it revalidates every saved root,
inverse-F profile, inverse-V phase, original conditions and self exemptions.
It then takes the global minimum. It does not rerun the profile domains.

At that minimum the combination checks the original endpoint x(tau), its
negative square, future orientation against q0, all known-wall inequalities,
and orthogonality to the eight retained walls plus a proposed new simple
wall. Its nine-normal Gram must be positive definite. Every minimum tie
must reconstruct in ALL ten coordinates as a nonnegative combination of
that basis. Exactly one supported simple cone is required. An unresolved
or nonunique cone would give no full-germ conclusion.

Completeness uses the domain theorem, not the fact that a finite loop ran.
Every actual first simple wall has a listed profile, nearest phase and
original coordinate witness, and passes every distinct known-wall test.
Conversely any admitted root is a mirror of the full original reflection
arrangement and cannot cut the interior of its starting chamber. The
complete admitted minimum is therefore exactly the first boundary. Its
complete tie cone and the inherited finite local recognition identify the
FULL germ. A feasible subdiagram or a minimum among only some profiles
would not suffice.

## 3. Actual runtime output

The source/runtime manifest is
`results/prime379_third_released_m_runtime_manifest_20261022.json`,
SHA-256 `ec8eb9bb4558703db8a2d9253bd7c6c425f71c2d0820426e6e565494593b2824`.
It binds all five NEW outputs, actual wrapper ledgers and raw stdout,
the separate combine output/ledger/stdout, the worker/contract/domain,
all fixed premise bindings and both full worker reviews. Root streamed
and matched every binding when preparing it. This metadata validation
does not rerun the algorithm or replace geometric review.

All55 selectors completed. The actual occurrence counters are

| Counter | Occurrences across all phases |
|---|---:|
| Scalar pairs visited |21,835|
| Inside scalar bounds |20,559|
| Inside folded window |14,795|
| Field time roots |132|
| Whole-quarter seeds |132|
| Original-integral occurrences |28|
| Admitted occurrences |11|

These are phase occurrence counts, not distinct root counts. The
authenticated combination has ONE distinct global minimum, profile
`norm2b_orthogonal`, phase2, with no known-wall self exemption. Its
complete output is
`results/prime379_third_released_m_combine_20261022.json`, SHA-256
`6d7627bbf2c2b7224aaa539eff2678da30ed061d5890e01a8a0f0d1cdee4917e`.
It reports `PASS_COMPLETE_DOMAIN`,55 selectors and ONE supported full
simple cone. The stored root is exactly(1): its E8 coordinates are
(2815372+4555369phi) times the displayed H. Its encoded collision is
exactly(2).

All six jobs ran serially through the unchanged shared-lock research
runner, with59seconds/1GiB sampled process-tree RSS/one thread. Each
passed40GiB disk/6GiB available-RAM admission;8GiB/3GiB abort floors stayed
active. Maximum elapsed time was0.607seconds and maximum sampled
process-tree RSS was61,095,936 bytes. No worker was imported outside the
runner. No build products were generated; parent and child used -B.

The earlier resource refusal, before any scientific child started, and
its empty historical stdout are preserved. They earn no execution credit.
The later five successful profiles and separate combine above provide the
runtime evidence. No successful unit was rerun merely for a status change.

## 4. Full finite germ and limitations

F fixes E7 and takes m to r. The winning inverse-F profile is long and
orthogonal to E7+m. Therefore its original wall w is orthogonal to
E7+r. The first seven retained normals have the E7 Gram, while r²=2
and w²=2b>0. Consequently the authenticated nine-normal Gram is exactly
E7 perpendicular [2] perpendicular [2b], proving the diagram type(3).
The checked future timelike endpoint and unique complete cone, together
with Section2's firstness/recognition argument, prove the full finite
corner statement rather than only a diagram found among known walls.

This is one additional actual edge of the literal norm379 chamber.
Reflectivity, finite volume and other outgoing edges remain unresolved.
No chamber-preserving isometry or period has been established, and no
claim for other lattices, form classes or fields follows. No Lean job
was run for this construction. The six frozen papers and all global
dimension ceilings are unchanged.
